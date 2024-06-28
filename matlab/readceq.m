function ceq = readceq(fn)

fid = fopen(fn, 'rb');

ceq.nMax          = fread(fid, 1, 'int32');
ceq.nParentBlocks = fread(fid, 1, 'int16');
ceq.nSegments     = fread(fid, 1, 'int16');

for ii = 1:ceq.nParentBlocks
    ceq.blockDuration = fread(fid, 1, 'float32');

    rf = sub_readrf(fid);
    ceq.parentBlocks{ii}.rf = rf;
end

fclose(fid);

function rf = sub_readrf(fid)

flg = fread(fid, 1, 'int16');
if flg
    rf.type = fread(fid, 1, 'int16');
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
