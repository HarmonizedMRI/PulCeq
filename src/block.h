#ifndef BLOCK_H
#define BLOCK_H

typedef struct {
    short nSamples;    /* number of 4us waveform samples */
    short* magnitude;  
    short* phase;     
} PULSEQ_SHAPE_ARBITRARY;
    
typedef struct {
    short riseTime;    /* us */
    short flatTime;    /* us */
    short fallTime;    /* us */
} PULSEQ_SHAPE_TRAP;

/*
typedef union {
    PULSEQ_SHAPE_ARBITRARY wav;
    PULSEQ_SHAPE_TRAP trap;
} PULSEQ_SHAPE;
*/

typedef struct {
    short type;                   /* NULL or ARBITRARY */
    short delay;                  /* us */
    float amplitude;              /* Gauss */
    PULSEQ_SHAPE_ARBITRARY wav;   /* arbitrary waveform in hardware units */
} PULSEQ_RF;

typedef struct {
    short type;          /* NULL, TRAP, or ARBITRARY */
    short delay;         /* us */
    float amplitude;     /* Gauss/cm */
    union {
        PULSEQ_SHAPE_ARBITRARY wav;
        PULSEQ_SHAPE_TRAP trap;
    } shape;
} PULSEQ_GRAD;

typedef struct {
    short type;         /* NULL or ADC */
    short numSamples;   /* number of ADC samples */
    short dwell;        /* us. Must be multiple of 2us */
    short delay;        /* us */
} PULSEQ_ADC;

typedef struct {
    short ID;
    short duration;
    PULSEQ_RF   rf;
    PULSEQ_GRAD gx;
    PULSEQ_GRAD gy;
    PULSEQ_GRAD gz;
    PULSEQ_ADC  adc;
} PULSEQ_BLOCK; 

int read_block(FILE* fid, PULSEQ_BLOCK* blk);  /* cf. pulsegeq.readblock() */
int read_grad(FILE* fid, PULSEQ_GRAD* grad);   /* cf. pulsegeq.readgrad() */

#endif
