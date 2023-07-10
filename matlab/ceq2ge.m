function ceq2ge(ceq, sysGE, ofname, verbose)
%
% Write a Ceq struct to a set of files that can be executed
% on GE scanners using the TOPPE interpreter (v6)

if nargin < 4
    verbose = false;
end

gamma = 42576000;   % Hz/T

% Identify pure delay blocks (handled different from other blocks in the TOPPE interpreter)
k = 1;
for p = 1:ceq.nParentBlocks
    if ~isdelayblock(ceq.parentBlocks{p})
        activeParentBlocks.IDs(k) = p;
        k = k + 1;
    end
end

nActiveParentBlocks = length(activeParentBlocks.IDs);

% define .mod file names
for k = 1:nActiveParentBlocks
    p = activeParentBlocks.IDs(k);
    modFiles{k} = sprintf('module%d.mod', k);
end

% Write modules.txt
fid = fopen('modules.txt','w');
fprintf(fid,'Total number of unique cores\n');
fprintf(fid,'%d\n', nActiveParentBlocks);
fprintf(fid,'fname dur(us) hasRF hasADC trigpos\n');

for k = 1:nActiveParentBlocks
    p = activeParentBlocks.IDs(k);
    b = ceq.parentBlocks{p};
    hasRF(k) = ~isempty(b.rf);
    hasADC(k) = ~isempty(b.adc);
    if (hasRF(k) & hasADC(k))
        error('Block cannot contain both RF and ADC events');
    end
    dur = round(ceil(b.blockDuration/sysGE.raster)*sysGE.raster*1e6); % us
    fprintf(fid,'%s\t%d\t%d\t%d\t-1\n', ...
        modFiles{k}, dur, hasRF(k), hasADC(k));    
end
fclose(fid);

% Write .mod files
for k = 1:nActiveParentBlocks
    p = activeParentBlocks.IDs(k);

    % defaults
    rf = []; grad.x = []; grad.y = []; grad.z = [];
    isDelayBlock = true;

    b = ceq.parentBlocks{p};

    % npre = number of 4us samples to discard at beginning of RF/ADC window
    % rfres = number of 4us samples in RF/ADC window
    if hasRF(k)
        npre = ceil(b.rf.delay/sysGE.raster);
        rfres = ceil(b.rf.shape_dur/sysGE.raster);
        b1ScalingFile = modFiles{k};
    elseif hasADC(k)
        npre = ceil(b.adc.delay/sysGE.raster);
        rfres = ceil(b.adc.numSamples*b.adc.dwell/sysGE.raster);
        readoutFile = modFiles{k};
    else
        npre = 0;
    end
    nChop = [npre 0];

    % Interpolate waveforms and convert to Gauss and Gauss/cm
    if hasRF(k)
        tge = sysGE.raster/2 : sysGE.raster : b.rf.shape_dur;
        rf = interp1(b.rf.t, b.rf.signal, tge) / gamma * 1e4;  % Gauss
        rf = [zeros(npre,1); rf.'];
        isDelayBlock = false;
    end

    for ax = {'x','y','z'}
        g = b.(['g' ax{1}]);
        if ~isempty(g)
            isDelayBlock = false;
            if strcmp(g.type, 'grad')
                % Arbitrary gradient
                tge = sysGE.raster/2 : sysGE.raster : max(g.tt);
                grad.(ax{1}) = interp1(g.tt, g.waveform, tge) / gamma * 100;   % Gauss/cm
            else
                % Convert trapezoid to arbitrary gradient
                gtmp = [ linspace(0, 0, ceil(g.delay/sysGE.raster)) ...
                    linspace(0, g.amplitude, ceil(g.riseTime/sysGE.raster)+1)  ...
                    g.amplitude*ones(1, floor(g.flatTime/sysGE.raster)) ...
                    linspace(g.amplitude, 0, ceil(g.fallTime/sysGE.raster)+1) ]';
                grad.(ax{1}) = gtmp / gamma * 100;   % Gauss/cm
            end
        end
    end
    
    % write .mod file
    if ~isDelayBlock 
        toppe.writemod(sysGE, 'ofname', modFiles{k}, ...
            'rf', rf, 'gx', grad.x', 'gy', grad.y', 'gz', grad.z', ...
            'nChop', nChop);
    end
end

% Write scanloop.txt
toppe.write2loop('setup', sysGE, 'version', 6); 

for n = 1:ceq.nMax
   if ~mod(n, 2000) | n == ceq.nMax
        for inb = 1:30
            fprintf('\b');
        end
        fprintf('Block %d/%d', n, ceq.nMax);
    end

    i = ceq.loop(n, 1);   % block group ID
    p = ceq.loop(n, 2);   % parent block ID

    k = find(activeParentBlocks.IDs == p);
    if isempty(k)
        % delay block
        % TODO: how to deal with large number of different delay times in ceq struct
        toppe.write2loop('delay', sysGE, ...
            'textra', round(ceq.parentBlocks{p}.blockDuration/sysGE.raster)*sysGE.raster*1e3, ...  % msec
            'core', i);
        continue;
    end

    if ceq.parentBlocks{p}.amp.rf > 0
        RFamplitude = ceq.loop(n, 3) / ceq.parentBlocks{p}.amp.rf;  % scaling, [-1 1]
    else
        RFamplitude = 0;
    end
    RFphase = ceq.loop(n, 4);   % rad
    RFoffset = ceq.loop(n, 5);  % Hz

    slice = 1; echo = 1; view = 1;  % TODO

    ax = {'gx','gy','gz'};
    for idim = 1:length(ax)
        if ceq.parentBlocks{p}.amp.(ax{idim}) > 0
            amp.(ax{idim}) = ceq.loop(n, 5+idim) / ceq.parentBlocks{p}.amp.(ax{idim});
        else
            amp.(ax{idim}) = 0;
        end
    end

    %rfamp rfphs rffreq amp.gx amp.gy amp.gz recphs];

    DAQphase = ceq.loop(n, 9);

    toppe.write2loop(modFiles{k}, sysGE, ...
        'Gamplitude',  [amp.gx amp.gy amp.gz]', ...
        'RFamplitude', RFamplitude, ...
        'RFphase',     RFphase, ...
        'DAQphase',    DAQphase, ...
        'RFoffset',    RFoffset, ...
        'slice',       slice, ...
        'echo',        echo, ...
        'view',        view, ...
        'dabmode',     'on', ...
        'textra',      0, ...  % TODO
        'waveform',    1, ...
        'core', i);

end
toppe.write2loop('finish', sysGE);   % close file
fprintf('\n');

% Write .entry file
toppe.writeentryfile('toppeN.entry', ...
    'filePath', '/usr/g/research/pulseq/v6/seq2ge/', ...
    'b1ScalingFile', b1ScalingFile, ...
    'readoutFile', readoutFile);

% write block group file (cores.txt)
for i = 1:ceq.nGroups
    blockIDs = ceq.groups.blockIDs;

    for j = 1:length(blockIDs)
        k = find(blockIDs(j) == activeParentBlocks.IDs);
        if isempty(k)
            blockIDs(j) = 0;
        else
            blockIDs(j) = k;
        end
    end

    blockGroups{i} = blockIDs;
end

toppe.writecoresfile(blockGroups);

% Create 'sequence stamp' file for TOPPE.
% TODO: update plotseq to handle delay blocks (mod ID = 0)
% This file is listed in line 6 of toppeN.entry
%toppe.preflightcheck('toppeN.entry', 'seqstamp.txt', sysGE);

% Put TOPPE files in a .tar file (for convenience)
system(sprintf('tar cf %s toppeN.entry seqstamp.txt modules.txt scanloop.txt cores.txt', ofname));
for k = 1:nActiveParentBlocks
    system(sprintf('tar rf %s %s', ofname, modFiles{k}));
end

% clean up (unless in verbose mode)
if ~verbose
    system('rm toppeN.entry seqstamp.txt modules.txt scanloop.txt cores.txt');
    for k = 1:nActiveParentBlocks
        system(sprintf('rm %s', modFiles{k}));
    end
end


