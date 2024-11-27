function ceq = seq2ceq(seqarg, varargin)
% function ceq = seq2ceq(seq)
%
% Convert a Pulseq file (http://pulseq.github.io/) to a PulCeq struct.
% See github/HarmonizedMRI/PulCeq/src/pulCeq.h
%
% Input
%   seqarg     a seq object, or name of a .seq file
%
% Input options with defaults
%   nMax                  [1]            Only parse the first nMax blocks in the .seq file  [all blocks]
%   ignoreSegmentLabels   true/false     Treat each block as a segment. Use with caution! [false]
%   verbose               true/false     Print some info to the terminal [false]
%
% Output
%   ceq        struct, based on github/HarmonizedMRI/PulCeq/src/pulCeq.h

GAM = 4257.6;   % Hz/Gauss

%% parse inputs

% defaults
arg.verbose = false;
arg.nMax    = [];
arg.ignoreSegmentLabels = false; % Don't define segments using the EXT column

% Substitute specified system values as appropriate (from MIRT toolbox)
arg = vararg_pair(arg, varargin);


%% Get seq object
if isa(seqarg, 'char')
    seq = mr.Sequence();
    seq.read(seqarg);
else
    if ~isa(seqarg, 'mr.Sequence')
        error('First argument is not an mr.Sequence object');
    end
    seq = seqarg;
end

nEvents = 7;   % Pulseq 1.4.0 
blockEvents = cell2mat(seq.blockEvents);
blockEvents = reshape(blockEvents, [nEvents, length(seq.blockEvents)]).'; 

% number of blocks (rows in .seq file) to step through
if isempty(arg.nMax)
    ceq.nMax = size(blockEvents, 1);
else
    ceq.nMax = arg.nMax;
end



%% Get parent blocks
% parent blocks = unique up to a scaling factor, or phase/frequency offsets.
% Contains waveforms with maximum amplitude across blocks.
% First find unique blocks, then determine and set max amplitudes.
% parentBlockIDs = [nMax], vector of parent block IDs for all blocks

parentBlockIndex = []; 

parentBlockIDs = [];

fprintf('seq2ceq: Getting block %d/%d', 1, ceq.nMax); prev_n = 1;
for n = 1:ceq.nMax
    if ~mod(n, 500) || n == ceq.nMax
        for ib = 1:strlength(sprintf('seq2ceq: Getting block %d/%d', prev_n, ceq.nMax))
            fprintf('\b');
        end
        prev_n = n;
        fprintf(sprintf('seq2ceq: Getting block %d/%d', n, ceq.nMax));
        if n == ceq.nMax, fprintf('\n'), end
    end

    b = seq.getBlock(n);

    % Skip blocks with zero duration (only contains label(s))
    if b.blockDuration < 2*eps
        parentBlockIDs(n) = -1;
        continue;
    end

    % Pure delay blocks are handled separately
    if isdelayblock(b)
        parentBlockIDs(n) = 0;
        continue;
    end

    if isempty(parentBlockIndex)
        parentBlockIndex(1) = n;
    end

    for p = 1:length(parentBlockIndex)
        n2 = parentBlockIndex(p);
        IsSame(p) = compareblocks(seq, blockEvents(n,:), blockEvents(n2,:), n, n2);
    end
    if sum(IsSame) == 0
        if arg.verbose
            fprintf('\nFound new parent block on line %d\n', n);
        end
        parentBlockIndex(p+1) = n;  % add new block to list
        parentBlockIDs(n) = p+1;
    else
        I = find(IsSame);
        parentBlockIDs(n) = I;
    end
end

ceq.nParentBlocks = length(parentBlockIndex);
for p = 1:length(parentBlockIndex)
    ceq.parentBlocks{p} = seq.getBlock(parentBlockIndex(p));
    ceq.parentBlocks{p}.ID = p;
end


%% Get segment (block group) definitions
% Segments are defined by their first occurrence in the .seq file
previouslyDefinedSegmentIDs = [];
segmentIDs = zeros(1,ceq.nMax);  % keep track of which segment each block belongs to
for n = 1:ceq.nMax
    b = seq.getBlock(n);

    if parentBlockIDs(n) == -1
        continue;    % skip
    end

    % Get segment ID label (TRID) if present
    if isfield(b, 'label') 
        hasTRIDlabel = false;
        for ii = 1:length(b.label)
            if strcmp(b.label(ii).label, 'TRID')
                hasTRIDlabel = true;
                break;
            end
        end
        if hasTRIDlabel  % marks start of segment
            activeSegmentID = b.label(ii).value;

            if ~any(activeSegmentID == previouslyDefinedSegmentIDs)
                % start new segment
                firstOccurrence = true;
                previouslyDefinedSegmentIDs = [previouslyDefinedSegmentIDs activeSegmentID];
                Segments{activeSegmentID} = [];
            else
                firstOccurrence = false;
            end
        end
    end

    if ~exist('firstOccurrence', 'var')
        %error('First block must contain a segment ID');
    end

    % add block to segment
    if firstOccurrence
        Segments{activeSegmentID} = [Segments{activeSegmentID} parentBlockIDs(n)];
    end

    segmentIDs(n) = activeSegmentID;
end

% In the above, the Segments array index equals the Segment ID specified in the .seq file,
% which is an arbitrary integer.
% Now we squash the Segments array and redefine the Segment IDs accordingly;
% this is needed since the interpreter assumes that segment ID = index into segment array.
if ~arg.ignoreSegmentLabels
    iSeg = 1;    % segment array index
    for segmentID = 1:length(Segments)
        if ~isempty(Segments{segmentID})
            segmentID2Ind(segmentID) = iSeg;
            ceq.segments(iSeg).segmentID = iSeg;
            ceq.segments(iSeg).nBlocksInSegment = length(Segments{segmentID});
            ceq.segments(iSeg).blockIDs = Segments{segmentID};
            iSeg = iSeg + 1;
        end
    end
end

ceq.nSegments = length(ceq.segments);


%% Get dynamic scan information, including cardiac trigger
ceq.loop = zeros(ceq.nMax, 14);
physioTrigger = false;
for n = 1:ceq.nMax
    b = seq.getBlock(n);
    if isfield(b, 'trig')
        if strcmp(b.trig.channel, 'physio1')
            physioTrigger = true;
        end
    end
    p = parentBlockIDs(n); 
    if p > -1 % if p == 0  % delay block
        ceq.loop(n,:) = getdynamics(b, segmentID2Ind(segmentIDs(n)), p, physioTrigger);
        physioTrigger = false;
    %end
    %if p > 0
    %    ceq.loop(n,:) = getdynamics(b, segmentID2Ind(segmentIDs(n)), p, ceq.parentBlocks{p});
    end
end


%% Remove zero-duration (label-only) blocks from ceq.loop
ceq.loop(parentBlockIDs == -1, :) = [];
ceq.nMax = size(ceq.loop,1);


%% No trigger for now. TODO
for p = 1:ceq.nParentBlocks
    ceq.parentBlocks{p}.trig.type = 0;
end


%% Check that the execution of blocks throughout the sequence
%% is consistent with the segment definitions
n = 1;
while n < ceq.nMax
    i = ceq.loop(n, 1);  % segment index

    if (n + ceq.segments(i).nBlocksInSegment) > ceq.nMax
        break;
    end

    % loop through blocks in segment
    for j = 1:ceq.segments(i).nBlocksInSegment

        % compare parent block id in ceq.loop against block id in ceq.segments(i)
        p = ceq.loop(n, 2);  % parent block id
        p_ij = ceq.segments(i).blockIDs(j);
        if p ~= p_ij
            warning(sprintf('Sequence contains inconsistent segment definitions. This may occur due to programming error (possibly fatal), or if an arbitrary gradient resembles that from another block except with opposite sign or scaled by zero (which is probably ok). Expected parent block ID %d, found %d (block %d)', p_ij, p, n));
        end

        n = n + 1;
    end
end

%% Gradient heating related calculations

% Get block/row index corresponding to the beginning of 
% the segment instance with the largest combined (all axes) gradient energy.

% initialize max energy field
for i = 1:ceq.nSegments
    ceq.segments(i).Emax.val = 0;
    ceq.segments(i).Emax.n = 1;
end
   
% find segment instance with max energy
n = 1;
while n < ceq.nMax
    % Calculate total energy in segment instance
    i = ceq.loop(n, 1);  % segment index
    Etmp.gx = 0; Etmp.gy = 0; Etmp.gz = 0;
    nFirst = n;
    for j = 1:ceq.segments(i).nBlocksInSegment  
        Etmp.gx = Etmp.gx + ceq.loop(n, 11);
        Etmp.gy = Etmp.gy + ceq.loop(n, 12);
        Etmp.gz = Etmp.gz + ceq.loop(n, 13);
        n = n + 1;
    end
    Etmp.all = Etmp.gx + Etmp.gy + Etmp.gz;

    % update Emax field
    if Etmp.all > ceq.segments(i).Emax.val
        ceq.segments(i).Emax.n = nFirst;
        ceq.segments(i).Emax.val = Etmp.all;
    end
end
