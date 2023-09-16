function ceq2ge(ceq, sysGE, ofname, varargin)
% function ceq2ge(ceq, sysGE, ofname, varargin)
%
% Write a Ceq struct to a set of files that can be executed
% on GE scanners using the TOPPE interpreter (v6)

% defaults
arg.verbose = false;
arg.ignoreTrigger = false;
arg.seqGradRasterTime = 10e-6;   % gradient raster time in .seq file

% Substitute specified system values as appropriate (from MIRT toolbox)
arg = vararg_pair(arg, varargin);

gamma = sysGE.gamma;   % Hz/T

raster = sysGE.raster*1e-6;   % sec
rasterUs = sysGE.raster;

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

    % number of 4us samples in block 
    %n = round(b.blockDuration/raster);

    % defaults
    % npre = number of 4us samples to discard at beginning of RF/ADC window
    % rfres = number of 4us samples in RF/ADC window
    npre = 0;
    nChopEnd = 0;   
    rfres = [];

    % Interpolate RF waveforms and convert to Gauss
    if hasRF(p)
        b1ScalingFile = modFiles{p};

        if b.rf.delay < sysGE.rfDeadTime*1e-6
            error(sprintf('Parent block %d: RF delay must be >= sysGE.rfDeadTime', p));
        end

        if b.rf.delay + b.rf.shape_dur + sysGE.rfRingdownTime*1e-6 > b.blockDuration
            error(sprintf('Parent block %d: RF ringdown extends past end of block', p));
        end

        tge = raster/2 : raster : b.rf.shape_dur;
        rf = interp1(b.rf.t, b.rf.signal, tge, 'linear', 'extrap') / gamma * 1e4;  % Gauss

        npre = round(b.rf.delay/raster);
        rfres = length(rf);

        rf = [zeros(npre,1); rf.'];
        
        isDelayBlock = false;
    end

    % Interpolate gradient waveforms and convert to Gauss/cm
    % We assume that there are only 3 types of gradients:
    % 1. trapezoid, starting and ending at 0
    % 2. extended trapezoid, specified on "corners" of shape. Assumed to start at zero if and only if delay > 0.
    % 3. arbitrary gradient specified on regular raster time. Assumed to start and end at zero ('first' and 'last' points). 
    for ax = {'x','y','z'}

        g = b.(['g' ax{1}]);

        if ~isempty(g)
            isDelayBlock = false;

            if strcmp(g.type, 'grad')
                % Arbitrary gradient
                % restore shape: if we had a
                % trapezoid converted to shape we have to find
                % the "corners" and we can eliminate internal
                % samples on the straight segments
                % but first we have to restore samples on the
                % edges of the gradient raster intervals
                % for that we need the first sample

                % check if we have an extended trapezoid, or an arbitrary gradient on a regular raster
                tt_rast=g.tt/arg.seqGradRasterTime+0.5;
                if all(abs(tt_rast-(1:length(tt_rast))')<1e-6)  % samples assumed to be on center of raster intervals
                    % arbitrary gradient on a regular raster
                    areaIn = sum(g.waveform)*arg.seqGradRasterTime;
                    g.first = 0;
                    g.last = 0;
                    wavtmp = [g.first g.waveform(:)' g.last];
                    tttmp = g.delay + [0 g.tt(:)' g.tt(end)+arg.seqGradRasterTime/2];
                else
                    % extended trapezoid: shape specified on "corners" of waveform
                    areaIn = sum( (g.waveform(1:(end-1)) + g.waveform(2:end))/2 .* diff(g.tt) );
                    wavtmp = g.waveform(:)';
                    tttmp = g.delay + g.tt(:)';
                end
            else
                % Convert trapezoid to arbitrary gradient
                areaIn = [g.riseTime/2 + g.flatTime + g.fallTime/2] * g.amplitude;
                if g.flatTime > 0
                    wavtmp = [0 1 1 0]*g.amplitude;       
                    tttmp = g.delay + [0 g.riseTime g.riseTime+g.flatTime g.riseTime+g.flatTime+g.fallTime ];
                else
                    wavtmp = [0 1 0]*g.amplitude;       
                    tttmp = g.delay + [0 g.riseTime g.riseTime+g.fallTime ];
                end
            end

            if g.delay > 0
                wavtmp = [0 wavtmp];
                tttmp = [0 tttmp];
            end

            tge = raster/2 : raster : tttmp(end);
            tmp = interp1(tttmp, wavtmp, tge);

            if any(isnan(tmp))
                msg = sprintf('NaN in gradient trapezoid waveform after interpolation (parent block %d)', p);
                error(msg);
            end

            % Scale to preserve area, and convert to Gauss/cm
            areaOut = sum(tmp) * raster;
            if areaIn > 1e-6 
                grad.(ax{1}) = tmp(:).' * areaIn/areaOut/ gamma * 100;   % Gauss/cm
            else
                % avoid division by 0
                grad.(ax{1}) = 0 * tmp(:).';
            end
        end
    end

    % ADC
    if hasADC(p)
        if b.adc.delay < sysGE.adcDeadTime*1e-6
            warning(sprintf('Parent block %d: ADC delay is < sysGE.adcDeadTime', p));
        end
        if b.adc.delay + b.adc.numSamples*b.adc.dwell > b.blockDuration
            error(sprintf('Parent block %d: ADC window extends past end of block', p));
        end
        npre = round(b.adc.delay/raster);
        rfres = round(b.adc.numSamples*b.adc.dwell/raster);
        readoutFile = modFiles{p};
    end

    % Set nChop, which is the number of samples to trim from beginning and end of RF/ADC window
    n = max([length(rf), length(grad.x), length(grad.y), length(grad.z)]);  % number of 4us samples in module
    if isempty(rfres)
        nChop = [0 0];
    else
        nChop = [npre n-npre-rfres];   % trim this many samples from beginning and end of RF/ADC window
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

    % trigger out
    if isfield(b, 'trig') & ~arg.ignoreTrigger
        if b.trig.delay < 100e-6
            warning('Requested trigger time too short. Setting to 100us');
            trigpos = 100;  % us
        else
            trigpos = round(b.trig.delay*1e6);   % us
        end
    else
        trigpos = -1;    % no trigger
    end

    rf = toppe.readmod(modFiles{p});
    dur = length(rf)*raster*1e6;  % us
    dur = max(dur, round(floor(b.blockDuration/raster)*raster*1e6)); % us
    dur = dur + sysGE.psd_rf_wait*hasRF(p);  % conservative/lazy choice for now
    fprintf(fid,'%s\t%d\t%d\t%d\t%d\n', ...
        modFiles{p}, round(dur), hasRF(p), hasADC(p), trigpos);    
end
fclose(fid);

%% write segment definition file (cores.txt) and determine TOPPE version
if ceq.nGroups > 0
    toppeVersion = 6;
    for i = 1:ceq.nGroups
        blockGroups{i} = ceq.groups(i).blockIDs;
    end
    toppe.writecoresfile(blockGroups);
else
    toppeVersion = 5;
end

%% Write scanloop.txt
% data frames (in P-file) are stored using indeces 'slice', 'echo', and 'view' 
sl = 1;
view = 1;
echo = 0; 
adcCount = 0;

toppe.write2loop('setup', sysGE, 'version', toppeVersion); 

fprintf('ceq2ge: Writing block %d/%d', 1, ceq.nMax); prev_n = 1; % Progress update trackers
for n = 1:ceq.nMax
    if ~mod(n, 2000) || n == ceq.nMax
        for ib = 1:strlength(sprintf('ceq2ge: Writing block %d/%d', prev_n, ceq.nMax))
            fprintf('\b');
        end
        prev_n = n;
        fprintf(sprintf('ceq2ge: Writing block %d/%d', n, ceq.nMax));
        if n == ceq.nMax, fprintf('\n'), end
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
    RFoffset = round(ceq.loop(n, 5));  % Hz

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

    % [rfamp rfphs rffreq amp.gx amp.gy amp.gz recphs]

    DAQphase = ceq.loop(n, 9);

    trigout = 1 * isfield(ceq.parentBlocks{p}, 'trig');

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
        'trigout',     trigout, ...
        'core', i);

end
toppe.write2loop('finish', sysGE);   % close file

% Write .entry file
toppe.writeentryfile('toppeN.entry', ...
    'filePath', '/usr/g/research/pulseq/v6/seq2ge/', ...
    'b1ScalingFile', b1ScalingFile, ...
    'readoutFile', readoutFile);

%% Create 'sequence stamp' file for TOPPE.
% TODO: update plotseq to handle delay blocks (mod ID = 0)
% This file is listed in line 6 of toppeN.entry
toppe.preflightcheck('toppeN.entry', 'seqstamp.txt', sysGE);

%% Put TOPPE files in a .tar file (for convenience)
if toppeVersion > 5
    system(sprintf('tar cf %s toppeN.entry seqstamp.txt modules.txt scanloop.txt cores.txt', ofname));
else
    system(sprintf('tar cf %s toppeN.entry seqstamp.txt modules.txt scanloop.txt', ofname));
end
for p = 1:ceq.nParentBlocks
    system(sprintf('tar rf %s %s', ofname, modFiles{p}));
end

%% clean up (unless in verbose mode)
if ~arg.verbose
    if toppeVersion > 5
        system('rm toppeN.entry seqstamp.txt modules.txt scanloop.txt cores.txt');
    else
        system('rm toppeN.entry seqstamp.txt modules.txt scanloop.txt');
    end
    for p = 1:ceq.nParentBlocks
        system(sprintf('rm %s', modFiles{p}));
    end
end

fprintf('Sequence file %s ready for execution on GE scanners\n', ofname);

