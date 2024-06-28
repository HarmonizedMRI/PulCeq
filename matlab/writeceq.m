function writeceq(ceq, fn)
% function writeceq(ceq, fn)
% 
% Write Ceq struct to binary file

fid = fopen(fn, 'wb');  % big endian (network byte order)

fwrite(fid, ceq.nMax, 'int32');
fwrite(fid, ceq.nParentBlocks, 'int16');
fwrite(fid, ceq.nSegments, 'int16');

% write parent blocks
for ii = 1:ceq.nParentBlocks

    b = ceq.parentBlocks{ii};

    fwrite(fid, b.blockDuration, 'float32');
    sub_writerf(fid, b.rf);
    sub_writegrad(fid, b.gx);
    sub_writegrad(fid, b.gy);
    sub_writegrad(fid, b.gz);
    sub_writeadc(fid, b.adc);
end

% write segment definitions
for ii = 1:ceq.nSegments
    sub_writesegment(fid, ceq.segments(ii));  % write definition of one segment
end

% write loop
sub_writeloop(fid, ceq.loop);

fclose(fid);

return

function sub_writerf(fid, rf)

if isempty(rf)
    fwrite(fid, 0, 'int16');   % flag indicating no rf event
else
    fwrite(fid, 1, 'int16');   % flag indicating that rf event is present
    fwrite(fid, 1, 'int16');                        % type. 1 = 'rf'
    fwrite(fid, length(rf.signal), 'int16');        % number of waveform samples
    fwrite(fid, abs(rf.signal),    'float32');      % Hz
    fwrite(fid, angle(rf.signal),  'float32');      % radians
    fwrite(fid, rf.t,              'float32');      % sec
    fwrite(fid, rf.shape_dur,      'float32');      % sec
    fwrite(fid, rf.delay,          'float32');      % sec
    fwrite(fid, rf.freqOffset,     'float32');      % Hz
    fwrite(fid, rf.phaseOffset,    'float32');      % radians
end

return

function sub_writegrad(fid, g)

% first value is flag indicating gradient type:
% 0: empty; 1: trap; 2: arbitrary gradient
if isempty(g)
    fwrite(fid, 0, 'int16');   % flag indicating no gradient event
else
    if strcmp(g.type, 'trap')
        fwrite(fid, 1, 'int16');   
        fwrite(fid, g.amplitude, 'float32');          % Hz/m
        fwrite(fid, g.riseTime, 'float32');           % sec
        fwrite(fid, g.flatTime, 'float32');           % sec
        fwrite(fid, g.fallTime, 'float32');           % sec
        fwrite(fid, g.delay,    'float32');           % sec
    else
        %fwrite(fid, 2, 'int16');   
        % TODO
    end
end

return

function sub_writeadc(fid, adc)

if isempty(adc)
    fwrite(fid, 0, 'int16');
else
    fwrite(fid, 1, 'int16');
    fwrite(fid, adc.numSamples, 'int32');
    fwrite(fid, adc.dwell, 'float32');
    fwrite(fid, adc.delay, 'float32');
    fwrite(fid, adc.freqOffset, 'float32');
    fwrite(fid, adc.phaseOffset, 'float32');
end

return

function sub_writesegment(fid, s)  % write definition of one segment
    fwrite(fid, s.segmentID, 'int16');
    fwrite(fid, s.nBlocksInSegment, 'int16');
    fwrite(fid, s.blockIDs, 'int16');
return

function sub_writeloop(fid, l)

for ii = 1:size(l,1)
    for jj = 1:size(l,2)  % write in row-major order
        fwrite(fid, l(ii,jj), 'float32');
    end
end

return
