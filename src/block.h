#ifndef BLOCK_H
#define BLOCK_H

typedef struct {
    short nSamples;      /* number of waveform samples */
    float* magnitude;  
    short* phase;     
} PulseqShapeArbitrary;
    
typedef struct {
    short riseTime;    /* us */
    short flatTime;    /* us */
    short fallTime;    /* us */
} PulseqShapeTrap;

typedef struct {
    short type;                   /* NULL or ARBITRARY */
    short delay;                  /* us */
    float amplitude;              /* Gauss */
    PulseqShapeArbitrary wav;   /* arbitrary waveform in hardware units */
} PulseqRF;

typedef struct {
    short type;          /* NULL, TRAP, or ARBITRARY */
    short delay;         /* us */
    float amplitude;     /* Gauss/cm */
    union {
        PulseqShapeArbitrary wav;
        PulseqShapeTrap trap;
    } shape;
} PulseqGrad;

typedef struct {
    int   type;         /* NULL or ADC */
    short numSamples;   /* number of ADC samples */
    short dwell;        /* us. Must be multiple of 2us */
    short delay;        /* us */
} PulseqADC;

typedef struct {
    float      duration;     /* sec */
    PulseqRF   rf;
    PulseqGrad gx;
    PulseqGrad gy;
    PulseqGrad gz;
    PulseqADC  adc;
    PulseqTrig trig;
    int ID;                 /* block ID (optional) */
    int* v1;                /* to be defined as needed by user (optional) */
    float * v2;             /* to be defined as needed by user (optional) */
} PulseqBlock; 

int read_block(FILE* fid, PulseqBlock* blk);  /* cf. pulsegeq.readblock() */
int read_grad(FILE* fid, PulseqGrad* grad);   /* cf. pulsegeq.readgrad() */

#endif
