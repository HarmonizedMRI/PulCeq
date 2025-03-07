function writeceq(ceq, fn, varargin)
% function writeceq(ceq, fn)
%
% Options
%   pislquant     default: 1      number of ADC events at start of scan for setting Rx gain
% 
% Write Ceq struct to binary file

% parse optional inputs
arg.pislquant = 1;  
arg.NMAXBLOCKSFORGRADHEATCHECK = 64000; 
arg = vararg_pair(arg, varargin);

fid = fopen(fn, 'wb');  % big endian (network byte order)

fwrite(fid, 47, 'int16');  % for checking endianness when reading 

% write base blocks
fwrite(fid, ceq.nParentBlocks, 'int16');
for ii = 1:ceq.nParentBlocks
    sub_writeblock(fid, ceq.parentBlocks(ii).block);
end

% write segment definitions
fwrite(fid, ceq.nSegments, 'int16');
for ii = 1:ceq.nSegments
    sub_writesegment(fid, ceq.segments(ii), ceq.parentBlocks);  % write definition of one segment
end

% write loop
fwrite(fid, ceq.nMax, 'int32');
fwrite(fid, size(ceq.loop,2), 'int16');   % nColumnsInLoopArray
for ii = 1:size(ceq.loop,1)
    fwrite(fid, ceq.loop(ii,:), 'float32');  % write in row-major order
end

% get max B1 and gradient in sequence
maxB1 = 0;
maxGrad = 0;
for n = 1:ceq.nMax
    maxB1 = max(maxB1, abs(ceq.loop(n,3)));  % Hz
    maxGrad = max([maxGrad abs(ceq.loop(n, [6 8 10]))]);
end

% safety stuff. some are dummy values, TODO
fwrite(fid, 1, 'float32');  % maxRfPower, G^2 * sec
fwrite(fid, maxB1, 'float32');
fwrite(fid, maxGrad, 'float32');   % maxGrad
fwrite(fid, 0.0, 'float32');   % maxSlew
fwrite(fid, ceq.duration, 'float32');   % duration
fwrite(fid, ceq.nReadouts, 'int32');    % total number of ADC events in sequence
fwrite(fid, arg.pislquant, 'int32');    % number ADC events at start of scan for setting receive gain in Auto Prescan

% Determine number of blocks (rows) to use for
% sliding window gradient/RF safety check
n = 1;
while n < min(ceq.nMax, arg.NMAXBLOCKSFORGRADHEATCHECK)
    i = ceq.loop(n, 1);
    %nBlocksForGreadHeatCheck = nBlocksForGreadHeatCheck + ceq.segments(i).nBlocksInSegment;
    n = n + ceq.segments(i).nBlocksInSegment;
    if n > arg.NMAXBLOCKSFORGRADHEATCHECK
        n = n - ceq.segments(i).nBlocksInSegment;
        break;
    end
end
n = n - 1;
fwrite(fid, n, 'int32');   

%{
    fread(&(ceq->maxRfPower), sizeof(float), 1, fid);
    fread(&(ceq->maxB1), sizeof(float), 1, fid);
    fread(&(ceq->maxGrad), sizeof(float), 1, fid);
    fread(&(ceq->maxSlew), sizeof(float), 1, fid);
    fread(&(ceq->duration), sizeof(float), 1, fid);
    fread(&(ceq->nReadouts), sizeof(int), 1, fid);
    fread(&(ceq->nGain), sizeof(int), 1, fid);
%}

fclose(fid);

return


function sub_writeblock(fid, b)

fwrite(fid, b.ID, 'int32');
fwrite(fid, b.blockDuration, 'float32');
sub_writerf(fid, b.rf);
sub_writegrad(fid, b.gx);
sub_writegrad(fid, b.gy);
sub_writegrad(fid, b.gz);
sub_writeadc(fid, b.adc);
trig.type = 0;
sub_writetrig(fid, trig);

return

function sub_writerf(fid, rf)

% flag indicating rf event type: 0: none; 1: arbitrary; 2: extended trap rf

if isempty(rf)
    fwrite(fid, 0, 'int16');   
    return
end

nSamples = length(rf.signal);
magnitude = rf.signal/max(abs(rf.signal));    % normalized waveform

shape = sub_rf2shape(rf);

if nSamples > 3 & all(diff(diff(rf.t)) < 10*eps)
    % waveform sampled on center of raster times
    assert(rf.t(1) > 10*eps, 'raster gradient: first sample time cannot be at time 0');
    assert(rf.t(2) - rf.t(1) > 0, 'raster gradient: detected non-positive raster time');

    fwrite(fid, 1, 'int16');   
    fwrite(fid, 1, 'int16');   % complex flag
    sub_writearbitrary(fid, shape, true, true);

else
    % waveform specified on corner points
    assert(abs(rf.t(1)) < 10*eps, 'extended trapezoid rf: first corner point must be sampled at time 0');

    fwrite(fid, 2, 'int16');   
    fwrite(fid, 1, 'int16');   % complex flag
    sub_writearbitrary(fid, shape, true, false);
end

fwrite(fid, rf.shape_dur, 'float32');
fwrite(fid, rf.delay, 'float32');
E = sum(shape.magnitude(1:end-1).^2 .* diff(rf.t));  % normalized pulse energy
fwrite(fid, E, 'float32');

return


function sub_writegrad(fid, g)

% first value is flag indicating gradient type:
% 0: empty; 1: trap; 2: raster; 3: cornerpoints
if isempty(g)
    fwrite(fid, 0, 'int16');   % flag indicating no gradient event
    return;
end

if strcmp(g.type, 'trap')
    fwrite(fid, 1, 'int16');   
    fwrite(fid, g.delay,    'float32');           % sec
    fwrite(fid, g.riseTime, 'float32');           % sec
    fwrite(fid, g.flatTime, 'float32');           % sec
    fwrite(fid, g.fallTime, 'float32');           % sec
    return;
end

nSamples = length(g.waveform);
magnitude = g.waveform/max(abs(g.waveform));    % normalized waveform

if g.tt(1) > 0
    % waveform sampled on center of raster times
    assert(nSamples > 2, 'raster gradient must have more than 2 samples');
    assert(all(abs(diff(diff(g.tt))) < 10*eps), 'raster gradient sample times must be regularly spaced');
    assert(abs(g.tt(1)) > 10*eps, 'raster gradient: first sample time cannot be at time 0');

    raster = g.tt(2) - g.tt(1);
    assert(raster > 0, 'raster gradient: detected negative raster time');

    fwrite(fid, 2, 'int16');   
    fwrite(fid, g.delay,   'float32');           % sec
    fwrite(fid, nSamples,  'int32');
    fwrite(fid, raster,    'float32');           % sec
    fwrite(fid, magnitude, 'float32');           % normalized

    return;
end

% waveform specified on corner points
assert(abs(g.tt(1)) < 10*eps, 'first corner point must be sampled at time 0');

fwrite(fid, 3, 'int16');   
fwrite(fid, g.delay,   'float32');    % sec
fwrite(fid, nSamples,  'int32');
fwrite(fid, 0,         'float32');    % dummy raster time, ignored by interpreter
fwrite(fid, g.tt,      'float32');    % sec
fwrite(fid, magnitude, 'float32');

return

function sub_writeadc(fid, adc)

if isempty(adc)
    fwrite(fid, 0, 'int16');
else
    fwrite(fid, 1, 'int16');
    fwrite(fid, adc.numSamples, 'int32');
    fwrite(fid, adc.dwell, 'float32');
    fwrite(fid, adc.delay, 'float32');
end

return

function sub_writetrig(fid, trig)
    % assume no trigger for now. TODO
    fwrite(fid, trig.type, 'int16');          
return

function sub_writesegment(fid, s, parentBlocks)  % write definition of one segment
    fwrite(fid, s.ID, 'int16');
    fwrite(fid, s.nBlocksInSegment, 'int16');
    fwrite(fid, s.blockIDs, 'int16');
    %fwrite(fid, s.Emax.val, 'float32');
    fwrite(fid, s.Emax.n, 'int32'); % block/row index of first block in segment instance with max gradient energy

    % determine if segment contains any gradients (needed by gradient heating check in interpreter)
    nGradEvents = 0;
    for j = 1:s.nBlocksInSegment
        p = s.blockIDs(j);
        if p > 0
            for ax = {'gx','gy','gz'}
                g = parentBlocks(p).block.(ax{1});
                if ~isempty(g)
                    nGradEvents = nGradEvents + 1; 
                end
            end
        end
    end
    fwrite(fid, 1*(nGradEvents>0), 'int16');          

return

function shape = sub_rf2shape(rf)
% Convert rf event to PulseqShapeArbitrary struct

%{
typedef struct {
    int nSamples;           /* Number of waveform samples */
    float raster;           /* Sample duration (sec) */
    float* time;            /* Time coordinates for waveform samples (sec) */
    float* magnitude;       /* Magnitude waveform (normalized) */
    float* phase;           /* Phase waveform (rad), only for type COMPLEX */
} PulseqShapeArbitrary;
%}

shape.nSamples = length(rf.signal);
shape.raster = rf.t(2) - rf.t(1);
shape.time = rf.t;
shape.magnitude = abs(rf.signal)/max(abs(rf.signal(:)));
shape.phase = angle(rf.signal);

return

function sub_writearbitrary(fid, shape, complexflag, regular_raster)

    fwrite(fid, shape.nSamples, 'int32');
    fwrite(fid, shape.raster, 'float32');

    if ~regular_raster
        fwrite(fid, shape.time, 'float32');
    end

    fwrite(fid, shape.magnitude, 'float32');

    if complexflag
        fwrite(fid, shape.phase, 'float32');
    end

return
