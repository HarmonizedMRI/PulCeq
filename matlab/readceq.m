function ceq = readceq(fn)
% function ceq = readceq(fn)
%
% Read Ceq struct from file created with writeceq.m

fid = fopen(fn, 'rb');

ceq.nMax          = fread(fid, 1, 'int32');
ceq.nParentBlocks = fread(fid, 1, 'int16');
ceq.nSegments     = fread(fid, 1, 'int16');

% read parent blocks
for ii = 1:ceq.nParentBlocks
    ceq.parentBlocks{ii}.blockDuration = fread(fid, 1, 'float32');
    ceq.parentBlocks{ii}.rf = sub_readrf(fid);
    ceq.parentBlocks{ii}.gx = sub_readgrad(fid);
    ceq.parentBlocks{ii}.gy = sub_readgrad(fid);
    ceq.parentBlocks{ii}.gz = sub_readgrad(fid);
    ceq.parentBlocks{ii}.adc = sub_readadc(fid);
end

% read segment definitions
for ii = 1:ceq.nSegments
    ceq.segments(ii) = sub_readsegment(fid);
end

% read loop
ceq.loop = sub_readloop(fid, ceq.nMax);

fclose(fid);

return

function rf = sub_readrf(fid)

rf.type = fread(fid, 1, 'int16');
if rf.type
    n     = fread(fid, 1, 'int16');         % number of waveform samples
    rho   = fread(fid, n, 'float32');
    theta = fread(fid, n, 'float32');
    rf.signal = rho.*exp(1i*theta);
    rf.t = fread(fid, n, 'float32');
    rf.shape_dur   = fread(fid, 1, 'float32');
    rf.delay       = fread(fid, 1, 'float32');
    rf.freqOffset  = fread(fid, 1, 'float32');
    rf.phaseOffset = fread(fid, 1, 'float32');
else
    rf = [];
end

return

function g = sub_readgrad(fid)

type = fread(fid, 1, 'int16');   % 0: empty; 1: trap; 2: arbitrary gradient

switch type
    case 0
        g = [];
    case 1
        g.amplitude = fread(fid, 1, 'float32');
        g.riseTime = fread(fid, 1, 'float32');
        g.flatTime = fread(fid, 1, 'float32');
        g.fallTime = fread(fid, 1, 'float32');
        g.delay    = fread(fid, 1, 'float32');
end

return

function adc = sub_readadc(fid)

type = fread(fid, 1, 'int16');

switch type
    case 0
        adc = [];
    case 1
        adc.numSamples  = fread(fid, 1, 'int32');
        adc.dwell       = fread(fid, 1, 'float32');
        adc.delay       = fread(fid, 1, 'float32');
        adc.freqOffset  = fread(fid, 1, 'float32');
        adc.phaseOffset = fread(fid, 1, 'float32');
end

return

function s = sub_readsegment(fid, s)  
    s.segmentID = fread(fid, 1, 'int16');
    s.nBlocksInSegment = fread(fid, 1, 'int16');
    s.blockIDs = fread(fid, s.nBlocksInSegment, 'int16');
return

function l = sub_readloop(fid, nMax)

nColumsInLoopArray = fread(fid, 1, 'int16');
l = zeros(nMax, nColumsInLoopArray);
for ii = 1:nMax
    l(ii,:) = fread(fid, nColumsInLoopArray, 'float32');
end

return
