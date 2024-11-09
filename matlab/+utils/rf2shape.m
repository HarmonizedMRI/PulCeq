function shape = rf2shape(rf)
% Convert rf event to PulseqShapeArbitrary struct

%{
typedef struct {
    int nSamples;           /* Number of waveform samples */
    float raster;           /* Sample duration (sec) */
    float* time;            /* Time coordinates for waveform samples (sec) */
    float* magnitude;       /* Magnitude waveform (normalized) */
    float* phase;           /* Phase waveform (rad), only for type COMPLEX */
    float amplitude;        /* Hz */
} PulseqShapeArbitrary;
%}

shape.nSamples = length(rf.signal);
shape.raster = rf.t(2) - rf.t(1);
shape.time = rf.t;
tmp = abs(rf.signal);
shape.magnitude = tmp/max(tmp);
shape.phase = angle(rf.signal);
shape.amplitude = max(shape.magnitude);
