#ifndef PULCEQ_H
#define PULCEQ_H

typedef struct {
    int nSamples;           /* number of waveform samples */
    float raster;           /* sample duration (sec) */
    float* magnitude;       /* magnitude waveform (normalized) */
    float* phase;           /* phase waveform (rad) */
} PulseqShapeArbitrary;
    
typedef struct {
    float riseTime;         /* sec */
    float flatTime;         /* sec */
    float fallTime;         /* sec */
} PulseqShapeTrap;

typedef struct {
    int type;                   /* NULL or ARBITRARY */
    float delay;                /* sec */
    float deadTime;             /* sec */
    float ringdownTime;         /* sec */
    PulseqShapeArbitrary wav;   /* arbitrary waveform, normalized amplitude */
} PulseqRF;

typedef struct {
    int type;                   /* NULL, TRAP, or ARBITRARY */
    float delay;                /* sec */
    union {
        PulseqShapeArbitrary wav;
        PulseqShapeTrap trap;
    } shape;
} PulseqGrad;

typedef struct {
    int   type;           /* NULL or ADC */
    int numSamples;       /* number of ADC samples */
    float dwell;          /* sec */
    float delay;          /* sec */
    float deadTime;       /* sec */
} PulseqADC;

typedef struct {
    int   type;           /* NULL or OUTPUT or ? */
    int   channel;        /* EXT1 or ? */
    float delay;          /* sec */
    float duration;       /* sec */
} PulseqTrig;

typedef struct {
    int ID;                 /* unique block ID */

    float      duration;    /* sec */
    PulseqRF   rf;
    PulseqGrad gx;
    PulseqGrad gy;
    PulseqGrad gz;
    PulseqADC  adc;
    PulseqTrig trig;

    /* vectors available for use as needed by the client program */
    int nVal1;                /* number of int values. Must be defined. */
    int* val1;                /* to be allocated dynamically by the client program */
    int nVal2;                /* number of float values. Must be defined. */
    float* nVal2;             /* to be allocated dynamically by the client program */
} PulseqBlock; 

/* Struct containing block IDs for all block groups */
typedef struct {
    int  nBlocksInGroup;    /* number of blocks in each group */
    int* blockIDs;          /* block id's for all groups */
} BlockGroup;

/* struct containing entire sequence definition */
typedef struct {
    int nParentBlocks; 
    PulseqBlock* parentBlocks;

    int nGroups;    
    BlockGroup* groups;      /* optional */

    float** loop    /* Dynamic scan settings: waveform amplitudes, phase offsets, etc. */
                    /* loop[n] = [rfamp gxamp gyamp gzamp rfFreqOffset rfPhaseOffset ...]
                       units:    [T     mT/m  mT/m  mT/m  Hz           rad           ...]
                    */

   int nMax;         /* number of blocks (rows in BLOCKS section) in .seq file */
} Seq;

#endif
