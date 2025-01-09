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
end

% Determine max amplitude across blocks
for p = 1:length(parentBlockIndex)
    ceq.parentBlocks{p}.amp.rf = 0;
    ceq.parentBlocks{p}.amp.gx = 0;  
    ceq.parentBlocks{p}.amp.gy = 0;
    ceq.parentBlocks{p}.amp.gz = 0;

    for n = 1:ceq.nMax
        if parentBlockIDs(n) ~= p | parentBlockIDs(n) == -1
            continue; 
        end
        block = seq.getBlock(n);
        if ~isempty(block.rf)
            ceq.parentBlocks{p}.amp.rf = max(ceq.parentBlocks{p}.amp.rf, max(abs(block.rf.signal)));
        end
        for ax = {'gx','gy','gz'}
            g = block.(ax{1});
            if ~isempty(g)
                if strcmp(g.type, 'trap')
                    gamp = abs(g.amplitude);   % trapezoids are positive lobes by definition
                else
                    gamp = max(abs(g.waveform));
                end
                ceq.parentBlocks{p}.amp.(ax{1}) = max(ceq.parentBlocks{p}.amp.(ax{1}), gamp);
            end
        end
    end
end

% Set parent block waveform amplitudes to max
for p = 1:ceq.nParentBlocks
    b = ceq.parentBlocks{p};   % shorthand
    if ~isempty(b.rf)
        ceq.parentBlocks{p}.rf.signal = b.rf.signal/max(abs(b.rf.signal))*b.amp.rf;
    end
    for ax = {'gx','gy','gz'}
        g = b.(ax{1});
        if ~isempty(g)
            if strcmp(g.type, 'trap')
                ceq.parentBlocks{p}.(ax{1}).amplitude = b.amp.(ax{1});
            else
                if max(abs(b.(ax{1}).waveform)) > 0
                    ceq.parentBlocks{p}.(ax{1}).waveform = b.(ax{1}).waveform/max(abs(b.(ax{1}).waveform))*b.amp.(ax{1});
                end
            end
        end
    end
end


%% Get segment (block group) definitions
% Segments are defined by their first occurrence in the .seq file
previouslyDefinedSegmentIDs = [];
if ~arg.ignoreSegmentLabels
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
else
    % Each block becomes its own segment (as in TOPPE v5)
    for p = 1:ceq.nParentBlocks
        ceq.groups(p).groupID = p;
        ceq.groups(p).nBlocksInGroup = 1;
        ceq.groups(p).blockIDs = p;
    end

    % Add delay segment (dedicated segment that's always defined)
    ceq.groups(p+1).groupID = p+1;
    ceq.groups(p+1).nBlocksInGroup = 1;
    ceq.groups(p+1).blockIDs = 0;   % block ID zero is a flag indicating a delay block

    segmentIDs = parentBlockIDs;
    segmentIDs(parentBlockIDs==0) = p+1;

    segmentID2Ind = 1:(p+1);
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
            ceq.groups(iSeg).groupID = iSeg;
            ceq.groups(iSeg).nBlocksInGroup = length(Segments{segmentID});
            ceq.groups(iSeg).blockIDs = Segments{segmentID};
            iSeg = iSeg + 1;
        end
    end
end

ceq.nGroups = length(ceq.groups);

%% Get dynamic scan information
ceq.loop = zeros(ceq.nMax, 11);
physioTrigger = false;
for n = 1:ceq.nMax
    b = seq.getBlock(n);
    % set cardiac trigger
    T = getblocktype(b);
    physioTrigger = T(3);
    p = parentBlockIDs(n); 
    if p == 0  % delay block
        ceq.loop(n,:) = getdynamics(b, segmentID2Ind(segmentIDs(n)), p);
    elseif p > 0
        ceq.loop(n,:) = getdynamics(b, segmentID2Ind(segmentIDs(n)), p, ceq.parentBlocks{p}, physioTrigger);
    end
end


%% Remove zero-duration (label-only) blocks from ceq.loop
ceq.loop(parentBlockIDs == -1, :) = [];
ceq.nMax = size(ceq.loop,1);


%% Check that the execution of blocks throughout the sequence
%% is consistent with the segment definitions
n = 1;
while n < ceq.nMax
    i = ceq.loop(n, 1);  % segment id

    if (n + ceq.groups(i).nBlocksInGroup) > ceq.nMax
        break;
    end

    % loop through blocks in segment
    for j = 1:ceq.groups(i).nBlocksInGroup

        % compare parent block id in ceq.loop against block id in ceq.groups(i)
        p = ceq.loop(n, 2);  % parent block id
        p_ij = ceq.groups(i).blockIDs(j);
        if p ~= p_ij
            warning(sprintf('Sequence contains inconsistent segment definitions. This may occur due to programming error (possibly fatal), or if an arbitrary gradient resembles that from another block except with opposite sign or scaled by zero (which is probably ok). Expected parent block ID %d, found %d (block %d)', p_ij, p, n));
        end

        n = n + 1;
    end
end
