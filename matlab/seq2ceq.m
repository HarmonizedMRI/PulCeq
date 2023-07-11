function ceq = seq2ceq(seqarg, varargin)
% function ceq = seq2ceq(seq)
%
% Convert a Pulseq file (http://pulseq.github.io/) to a PulCeq struct.
% See github/HarmonizedMRI/PulCeq/.
%
% Input
%   seqarg     a seq object, or name of a .seq file
%
% Output
%   ceq        struct, similar to github/HarmonizedMRI/PulCeq/src/pulCeq.h

%% parse inputs
% Defaults
arg.verbose = false;
arg.debug = false;
arg.nt      = [];

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
if isempty(arg.nt)
    ceq.nMax = size(blockEvents, 1);
else
    ceq.nMax = arg.nt;
end


%% Get parent blocks
% parent blocks = unique up to a scaling factor, or phase/frequency offsets.
% Contains waveforms with maximum amplitude across blocks.
% First find unique blocks, then determine and set max amplitudes.
% parentBlockIDs = vector of parent block IDs for all blocks

parentBlockIndex(1) = 1;  % first block is unique by definition

for n = 1:ceq.nMax
    if ~mod(n, 500) | n == ceq.nMax
        for inb = 1:20
            fprintf('\b');
        end
        fprintf('Block %d/%d', n, ceq.nMax);
    end

    issame = 0;

    % Pure delay blocks are handled separately
    b = seq.getBlock(n);
    if isdelayblock(b)
        parentBlockIDs(n) = 0;
        continue;
    end

    for p = 1:length(parentBlockIndex)
        n2 = parentBlockIndex(p);
        IsSame(p) = compareblocks(seq, blockEvents(n,:), blockEvents(n2,:), n, n2);
    end
    if sum(IsSame) == 0
        if arg.verbose
            fprintf('Found new block on line %d\n', n);
        end
        parentBlockIndex(p+1) = n;  % found a unique block, so add it to list
        parentBlockIDs(n) = p+1;
    else
        I = find(IsSame);
        parentBlockIDs(n) = I;
    end
end
fprintf('\n');

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
        if parentBlockIDs(n) ~= p
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
                ceq.parentBlocks{p}.(ax{1}).waveform = b.(ax{1}).waveform/max(abs(b.(ax{1}).waveform))*b.amp.(ax{1});
            end
        end
    end
end


%% Get block group definitions
currentCoreID = []; 
blockGroupIDs = zeros(1,ceq.nMax);  % keep track of which core each block belongs to
for n = 1:ceq.nMax
    if ~mod(n, 500) | n == ceq.nMax
        for inb = 1:20
            fprintf('\b');
        end
        fprintf('Block %d/%d', n, ceq.nMax);
    end

    b = seq.getBlock(n);

    if isfield(b, 'label')
        % wrap up the current core
        if ~isempty(currentCoreID)  
            Cores{currentCoreID} = [currentCoreID length(blockIDs) blockIDs];
        end

        % start new core
        currentCoreID = b.label.value;
        blockIDs = parentBlockIDs(n); 
    else
        % add block to this block group
        blockIDs = [blockIDs parentBlockIDs(n)];
    end

    blockGroupIDs(n) = currentCoreID;
end
fprintf('\n');

ceq.nGroups = length(Cores);
for i = 1:ceq.nGroups
    ceq.groups(i).groupID = Cores{i}(1);
    ceq.groups(i).nBlocksInGroup = Cores{i}(2);
    ceq.groups(i).blockIDs = Cores{i}(3:end);
end


%% Get dynamic scan information
ceq.loop = zeros(ceq.nMax, 10);
for n = 1:ceq.nMax
    b = seq.getBlock(n);
    p = parentBlockIDs(n); 
    if p == 0  % delay block
        ceq.loop(n,:) = getdynamics(b, blockGroupIDs(n), p);
    else
        ceq.loop(n,:) = getdynamics(b, blockGroupIDs(n), p, ceq.parentBlocks{p});
    end
end

