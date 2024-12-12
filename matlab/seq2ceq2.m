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

% Definitions:
% n, row        row index in .seq file
% i             segment array index, starting from 1
% j             block number within a segment, starting from 1

GAM = 4257.6;   % Hz/Gauss

minDelayBlockDuration = 20e-6 + 100*eps; 

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


%% Get all TRID labels and corresponding row indeces
%tridLabels = -1 * ones(1, ceq.nMax);
nTRIDlabels = 0;
for n = 1:ceq.nMax
    b = seq.getBlock(n);
    if isfield(b, 'label') 
        hasTRIDlabel = false;
        for ii = 1:length(b.label)
            if strcmp(b.label(ii).label, 'TRID')
                nTRIDlabels = nTRIDlabels + 1;
                tridLabels(nTRIDlabels) = b.label(ii).value;
                tridLabelsIndex(nTRIDlabels) = n;
                break;
            end
        end
    end
end
%tridLabels = tridLabels(1:nTRIDlabels);


%% Get list of (virtual) segments
%% These are distinct from segment 'instances'.

[uniqueTridLabels, I] = unique(tridLabels);
nBlocksPerTridLabel = diff(tridLabelsIndex); % this misses the last segment instance -- TODO
ceq.nSegments = length(uniqueTridLabels);
for i = 1:ceq.nSegments
    ceq.segments(i).nBlocksInSegment = nBlocksPerTridLabel(I(i));
    ceq.segments(i).segmentID = i;
    ceq.segments(i).rows = tridLabelsIndex(I(i)) + [0:ceq.segments(i).nBlocksInSegment-1];
end


%% Get parent blocks, by parsing first instance of each segment.
%% Also fill in the sequence of parent blocks for each segment.
%% Pure delay blocks are assigned parent block ID = 0

ceq.parentBlocks(1).row = -1;
ceq.nParentBlocks = 0;

for i = 1:ceq.nSegments

    for j = 1:ceq.segments(i).nBlocksInSegment

        n = ceq.segments(i).rows(j);

        b = seq.getBlock(n);
        T = getblocktype(b);

        % Skip blocks with "zero" duration, defined as 20us or less
        if b.blockDuration < minDelayBlockDuration
            continue;
        end

        % The first non-zero block is a parent block by definition.
        % This can be a pure delay block, or a block with rf/gradient/adc events.
        if ceq.parentBlocks(1).row == -1
            ceq.nParentBlocks = ceq.nParentBlocks + 1;
            ceq.parentBlocks(ceq.nParentBlocks).row = n;
            ceq.parentBlocks(ceq.nParentBlocks).block = b;
            ceq.parentBlocks(ceq.nParentBlocks).ID = ceq.nParentBlocks;
            ceq.segments(i).blockIDs(j) = ceq.nParentBlocks;
            continue;
        end

        % Check if block is similar to an existing parent block
        issame = false;
        for p = 1:ceq.nParentBlocks
            np = ceq.parentBlocks(p).row; 
            if compareblocks(seq, blockEvents(n,:), blockEvents(np,:), n, np)
                issame = true;
                ceq.segments(i).blockIDs(j) = ceq.parentBlocks(p).ID;
                pParent = p;
                break;
            end
        end

        % If not similar, add as a new parent block
        if ~issame
            if arg.verbose
                fprintf('\nFound new parent block on line %d\n', n);
            end
            ceq.nParentBlocks = ceq.nParentBlocks + 1;
            ceq.parentBlocks(ceq.nParentBlocks).row = n;
            ceq.parentBlocks(ceq.nParentBlocks).block = b;
            ceq.parentBlocks(ceq.nParentBlocks).ID = ceq.nParentBlocks;
            ceq.segments(i).blockIDs(j) = ceq.nParentBlocks;
        end
    end
end


%% Get dynamic scan information, including cardiac trigger
ceq.loop = zeros(ceq.nMax, 14);
physioTrigger = false;
activeSegmentID = [];
n = tridLabelsIndex(1);  % start of first segment instance
while n < ceq.nMax
    b = seq.getBlock(n);

    % set cardiac trigger
    T = getblocktype(b);
    physioTrigger = T(3);

    p = parentBlockIDs(n); 
    if p > -1 
        ceq.loop(n,:) = getdynamics(b, segmentID2Ind(segmentIDs(n)), p, physioTrigger);
        physioTrigger = false;
    end
end


%% Remove zero-duration (label-only) blocks from ceq.loop
ceq.loop(parentBlockIDs == -1, :) = [];
ceq.nMax = size(ceq.loop,1);


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
        msg = ['Sequence contains inconsistent segment definitions. ' ...
               'This may occur due to programming error (possibly fatal), ' ...
               'or if an arbitrary gradient resembles that from another block ' ...
               'except with opposite sign or scaled by zero (which is probably ok). ' ...
               'Often, a solution to this is to scale gradients to "eps" instead of ' ...
               'identically zero, when calling mr.scaleGrad().'];
        if p ~= p_ij
            %warning(sprintf('Sequence contains inconsistent segment definitions. This may occur due to programming error (possibly fatal), or if an arbitrary gradient resembles that from another block except with opposite sign or scaled by zero (which is probably ok). Often, a solution to this is to scale gradients to "eps" instead of identically zero, when calling mr.scaleGrad(). Expected parent block ID %d, found %d (block %d)', p_ij, p, n));
            warning(sprintf('%s\nExpected parent block ID %d, found %d (block %d)', msg, p_ij, p, n));
        end

        n = n + 1;
    end
end

%% Gradient heating related calculations

% Get block/row index corresponding to the beginning of the segment instance
% instance with the largest combined (all axes) gradient energy.

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
