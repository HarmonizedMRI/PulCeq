function ceq2ge(ceq, sysGE, ofname, verbose)
%
% Write a Ceq struct to a set of files that can be executed
% on GE scanners using the TOPPE interpreter (v6)

if nargin < 4
    verbose = false;
end

gamma = sysGE.gamma;   % Hz/T

raster = sysGE.raster*1e-6;   % sec

% define .mod file names
for p = 1:ceq.nParentBlocks
    modFiles{p} = sprintf('module%d.mod', p);
end

%% Write .mod files
for p = 1:ceq.nParentBlocks

    % defaults
    rf = []; grad.x = []; grad.y = []; grad.z = [];
    isDelayBlock = true;

    b = ceq.parentBlocks{p};

    hasRF(p) = ~isempty(b.rf);
    hasADC(p) = ~isempty(b.adc);

    % npre = number of 4us samples to discard at beginning of RF/ADC window
    % rfres = number of 4us samples in RF/ADC window
    if hasRF(p)
        npre = ceil(b.rf.delay/raster);
        rfres = ceil(b.rf.shape_dur/raster);
        b1ScalingFile = modFiles{p};
    elseif hasADC(p)
        npre = ceil(b.adc.delay/raster);
        rfres = ceil(b.adc.numSamples*b.adc.dwell/raster);
        readoutFile = modFiles{p};
    else
        npre = 0;
    end
    nChop = [npre 0];

    % Interpolate waveforms and convert to Gauss and Gauss/cm
    if hasRF(p)
        tge = raster/2 : raster : b.rf.shape_dur;
        rf = interp1(b.rf.t, b.rf.signal, tge, 'linear', 'extrap') / gamma * 1e4;  % Gauss
        rf = [zeros(npre,1); rf.'];
        if any(isnan(rf))
            keyboard
        end
        isDelayBlock = false;
    end

    for ax = {'x','y','z'}
        g = b.(['g' ax{1}]);
        if ~isempty(g)
            isDelayBlock = false;
            if strcmp(g.type, 'grad')
                % Arbitrary gradient
                tge = raster/2 : raster : max(g.tt);
                grad.(ax{1}) = interp1(g.tt, g.waveform, tge, 'linear', 'extrap') / gamma * 100;   % Gauss/cm
            else
                % Convert trapezoid to arbitrary gradient
                gtmp = [ linspace(0, 0, ceil(g.delay/raster)) ...
                    linspace(0, g.amplitude, ceil(g.riseTime/raster)+1)  ...
                    g.amplitude*ones(1, floor(g.flatTime/raster)) ...
                    linspace(g.amplitude, 0, ceil(g.fallTime/raster)+1) ]';
                grad.(ax{1}) = gtmp / gamma * 100;   % Gauss/cm
            end
        end
    end
    
    % write .mod file
    if ~isDelayBlock 
        toppe.writemod(sysGE, 'ofname', modFiles{p}, ...
            'rf', rf(:), 'gx', grad.x(:), 'gy', grad.y(:), 'gz', grad.z(:), ...
            'nChop', nChop);
    end
end

%% Write modules.txt
% Do this after creating .mod files, so .mod file duration can be set exactly.
fid = fopen('modules.txt','w');
fprintf(fid,'Total number of unique cores\n');
fprintf(fid,'%d\n', ceq.nParentBlocks);
fprintf(fid,'fname dur(us) hasRF hasADC trigpos\n');

for p = 1:ceq.nParentBlocks
    b = ceq.parentBlocks{p};
    if (hasRF(p) & hasADC(p))
        error('Block cannot contain both RF and ADC events');
    end
    rf = toppe.readmod(modFiles{p});
    dur = length(rf)*raster*1e6;  % us
    dur = max(dur, round(ceil(b.blockDuration/raster)*raster*1e6)); % us
    dur = dur + sysGE.psd_rf_wait*hasRF(p);  % conservative/lazy choice for now
    fprintf(fid,'%s\t%d\t%d\t%d\t-1\n', ...
        modFiles{p}, round(dur), hasRF(p), hasADC(p));    
end
fclose(fid);

%% Write scanloop.txt
% data frames (in P-file) are stored using indeces 'slice', 'echo', and 'view' 
sl = 1;
view = 1;
echo = 0; 
adcCount = 0;

toppe.write2loop('setup', sysGE, 'version', 6); 

for n = 1:ceq.nMax
    if ~mod(n, 2000) | n == ceq.nMax
        msg = sprintf('\tblock %d/%d', n, ceq.nMax);
        for ib = 1:strlength(msg)
            fprintf('\b');
        end
        fprintf(msg);
    end

    i = ceq.loop(n, 1);   % block group ID
    p = ceq.loop(n, 2);   % parent block ID

    if p == 0  % delay block
        toppe.write2loop('delay', sysGE, ...
            'textra', round(ceq.loop(n, 10)*1e6)/1e3, ... % msec
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

    % set slice/echo/view indeces (if block is an acquisition block)
    % view = 1, ..., system.maxView
    % sl   = 1, ..., system.maxSlice
    if hasADC(p)
        view = mod(adcCount, sysGE.maxView) + 1;
        sl   = floor(adcCount/sysGE.maxView) + 1;
        if sl > sysGE.maxSlice;
            error(sprintf('max number of slices ecxeeded (%d)', sysGE.maxSlice));
        end
        echo = floor(adcCount/(sysGE.maxView*sysGE.maxSlice));
        if echo > sysGE.maxEcho
            error(sprintf('max number of echoes ecxeeded (%d)', sysGE.maxEcho));
        end
        %fprintf('n: %d, view: %d, sl: %d, echo: %d\n', n, view, sl, echo);

        adcCount = adcCount+1;
    end

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

    toppe.write2loop(modFiles{p}, sysGE, ...
        'Gamplitude',  [amp.gx amp.gy amp.gz]', ...
        'RFamplitude', RFamplitude, ...
        'RFphase',     RFphase, ...
        'DAQphase',    DAQphase, ...
        'RFoffset',    RFoffset, ...
        'slice',       sl, ...
        'echo',        echo+1, ...  % write2loop starts indexing at 1
        'view',        view, ...
        'dabmode',     'on', ...
        'textra',      0, ...  
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
    blockGroups{i} = ceq.groups(i).blockIDs;
end
toppe.writecoresfile(blockGroups);

%% Create 'sequence stamp' file for TOPPE.
% TODO: update plotseq to handle delay blocks (mod ID = 0)
% This file is listed in line 6 of toppeN.entry
toppe.preflightcheck('toppeN.entry', 'seqstamp.txt', sysGE);

%% Put TOPPE files in a .tar file (for convenience)
system(sprintf('tar cf %s toppeN.entry seqstamp.txt modules.txt scanloop.txt cores.txt', ofname));
for p = 1:ceq.nParentBlocks
    system(sprintf('tar rf %s %s', ofname, modFiles{p}));
end

%% clean up (unless in verbose mode)
if ~verbose
    system('rm toppeN.entry seqstamp.txt modules.txt scanloop.txt cores.txt');
    for p = 1:ceq.nParentBlocks
        system(sprintf('rm %s', modFiles{p}));
    end
end

fprintf('\n');

