% write3DGRE.m
%
% 3D GRE demo sequence for Pulseq on GE v1.0 User Guide

% System/design parameters.
sys = mr.opts('maxGrad', 49, 'gradUnit','mT/m', ...
              'maxSlew', 199, 'slewUnit', 'T/m/s', ...
              'rfDeadTime', 100e-6, ...
              'rfRingdownTime', 60e-6, ...
              'adcDeadTime', 40e-6, ...
              'adcRasterTime', 2e-6, ...
              'rfRasterTime', 10e-6, ...   % GE: must be multiple of 2us
              'gradRasterTime', 4e-6, ...
              'blockDurationRaster', 4e-6, ...
              'B0', 3.0);

% Acquisition parameters
fov = [200e-3 200e-3 10e-3];   % FOV (m)
Nx = 200; Ny = 128; Nz = 4;    % Matrix size
slabThickness = 10e-3;         % slice thickness (m)
TR = 15e-3;                     % sec
dwell = 10e-6;                  % ADC sample time (s)
alpha = 90;                      % flip angle (degrees)
nCyclesSpoil = 2;               % number of spoiler cycles
Tpre = 2.0e-3;                  % prephasing trapezoid duration
rfSpoilingInc = 117;            % RF spoiling increment

% Create a new sequence object
seq = mr.Sequence(sys);           

% Create slab-selective pulse (waveform sampled on regular raster)
[rf] = mr.makeSincPulse(alpha*pi/180, 'Duration', 3e-3, ...
    'SliceThickness', slabThickness, 'apodization', 0.42, ...
    'timeBwProduct', 4, 'system', sys);

% Non-selective hard pulse (waveform sampled on corner points)
[rf2] = mr.makeBlockPulse(alpha/180*pi, sys, 'Duration', 0.5e-3);

% Define spin-warp gradients and ADC events
% Cut the redaout gradient into two parts for optimal spoiler timing
deltak = 1./fov;
Tread = Nx*dwell;

gyPre = mr.makeTrapezoid('y', sys, ...
    'Area', Ny*deltak(2)/2, ...   % PE1 gradient, max positive amplitude
    'Duration', Tpre);
gzPre = mr.makeTrapezoid('z', sys, ...
    'Area', Nz*deltak(3)/2, ...   % PE2 gradient, max positive amplitude
    'Duration', Tpre);

gxtmp = mr.makeTrapezoid('x', sys, ...  % readout trapezoid, temporary object
    'Amplitude', Nx*deltak(1)/Tread, ...
    'FlatTime', Tread);
gxPre = mr.makeTrapezoid('x', sys, ...
    'Area', -gxtmp.area/2, ...
    'Duration', Tpre);

adc = mr.makeAdc(Nx, sys, ...
    'Duration', Tread,...
    'Delay', gxtmp.riseTime);

% extend flat time so we can split at end of ADC dead time
gxtmp2 = mr.makeTrapezoid('x', sys, ...  % temporary object
    'Amplitude', Nx*deltak(1)/Tread, ...
    'FlatTime', Tread + adc.deadTime);   
[gx, ~] = mr.splitGradientAt(gxtmp2, gxtmp2.riseTime + gxtmp2.flatTime);

gzSpoil = mr.makeTrapezoid('z', sys, ...
    'Area', Nx*deltak(1)*nCyclesSpoil);
gxSpoil = mr.makeExtendedTrapezoidArea('x', gxtmp.amplitude, 0, gzSpoil.area, sys);

% y/z PE steps
pe1Steps = ((0:Ny-1)-Ny/2)/Ny*2;
pe2Steps = ((0:Nz-1)-Nz/2)/Nz*2;

% Calculate TR delay
TRmin = mr.calcDuration(rf) + mr.calcDuration(gxPre) ...
   + mr.calcDuration(gx) + mr.calcDuration(gxSpoil);
delayTR = TR - TRmin;

% make a spiral gradient
nleaf = 8; dt = 4e-10; 
[sp.wav] = getspiral(nleaf, sys.gradRasterTime, fov(1)*100, Nx);
sp.gx = mr.makeArbitraryGrad('x', real(sp.wav)*1e-4*sys.gamma*100, sys, ...
                        'delay', sys.adcDeadTime);
sp.gy = mr.makeArbitraryGrad('y', imag(sp.wav)*1e-4*sys.gamma*100, sys, ...
                        'delay', sys.adcDeadTime);

% make a cardiac trigger event
trig = mr.makeTrigger('physio1', 'duration', 20e-6);  % dummy duration -- ignored by GE interpreter


%%%%%%%%%%%% Start adding blocks to sequence %%%%%%%%%%%%%%%%%%

% add segment with different rf pulses
for ii = 1:0
    seq.addBlock(rf, mr.makeLabel('SET', 'TRID', 1), mr.makeDelay(4e-3));
    seq.addBlock(rf2, gxPre);
    seq.addBlock(gx);
    seq.addBlock(gxSpoil);
end

% Acquire 3D GRE sequence
% iZ < 0: Dummy shots to reach steady state
% iZ = 0: ADC is turned on and used for receive gain calibration on GE scanners
% iZ > 0: Image acquisition

nDummyZLoops = 0;

rf_phase = 0;
rf_inc = 0;

for iZ = -nDummyZLoops:Nz
    isDummyTR = iZ < 0;

    msg = sprintf('z encode %d of %d   ', iZ, Nz);
    for ibt = 1:(length(msg) + 2)
        fprintf('\b');
    end
    fprintf(msg);

    % add a cardiac trigger event for testing
    %seq.addBlock(trig);

    for iY = 1:Ny
        % Turn on y and z prephasing lobes, except during dummy scans and
        % receive gain calibration (auto prescan)
        yStep = (iZ > 0) * pe1Steps(iY) + eps;
        zStep = (iZ > 0) * pe2Steps(max(1,iZ)) + eps;

        % RF spoiling
        rf.phaseOffset = rf_phase/180*pi;
        adc.phaseOffset = rf_phase/180*pi;
        rf_inc = mod(rf_inc+rfSpoilingInc, 360.0);
        rf_phase = mod(rf_phase+rf_inc, 360.0);
        
        % Excitation
        % Mark start of segment (block group) by adding label.
        % Subsequent blocks in block group are NOT labelled.
        seq.addBlock(rf, mr.makeLabel('SET', 'TRID', 47-isDummyTR), trig);
        
        % Encoding
        seq.addBlock(gxPre, ...
            mr.scaleGrad(gyPre, yStep), ...
            mr.scaleGrad(gzPre, zStep));
        if isDummyTR
            seq.addBlock(gx);
        else
            seq.addBlock(gx, adc);
        end

        % rephasing/spoiling and TR delay
        seq.addBlock(gxSpoil, ...
            mr.scaleGrad(gyPre, -yStep), ...
            mr.scaleGrad(gzPre, -zStep));
        seq.addBlock(mr.makeDelay(delayTR));
    end
end

% add a spiral segment
for ii = 1:10
    seq.addBlock(sp.gx, sp.gy, mr.makeLabel('SET', 'TRID', 3), mr.makeDelay(50e-3));
    %seq.addBlock(mr.makeDelay(82.6e-3));
end


fprintf('Sequence ready\n');

% Check sequence timing
[ok, error_report]=seq.checkTiming;
if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

% Output for execution
seq.setDefinition('FOV', fov);
seq.setDefinition('Name', 'gre');
seq.write('test.seq');

%% Optional plots

% Plot sequence
Noffset = Ny*(nDummyZLoops+1);
%seq.plot('timerange',[Noffset Noffset+4]*TR, 'timedisp', 'ms');

return

% Plot k-space (2d)
[ktraj_adc,t_adc,ktraj,t_ktraj,t_excitation,t_refocusing] = seq.calculateKspacePP();
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D k-space plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
title('full k-space trajectory (k_x x k_y)');
