#include <string.h>
#include <math.h>

#include "sdkver.h" 

#ifndef HW_IO
#include <stdio.h>
#else
#include <stdioLib.h>
#endif

#if defined MR29_0_R01 || defined MR29_1_R01
#include <allocnode.h>
#endif
#if defined RX28_0_R05
#include "pulsegen/allocnode.h"
#endif

#include "jfn_globaldefs.h"   
#include "block.h"
#include "modules.h"  /* readshort() */
#include "constants.h"  /* readshort() */

/* Read a block file */
int read_block(FILE* fid, PULSEQ_BLOCK* block) {

    short type;
    short n;     
    short amp_int16;

    readshort(&(block->ID), 1, fid);
    readshort(&(block->duration), 1, fid);

    /* rf */
    readshort(&((block->rf).type), 1, fid);
    if ((block->rf).type == ARBITRARY) {
        readshort(&((block->rf).delay), 1, fid);     /* us */

        readshort(&(amp_int16), 1, fid);
        (block->rf).amplitude = 1.0 * amp_int16/RFSCALE;  /* Gauss */

        readshort(&n, 1, fid);
        (block->rf).wav.nSamples = n;       
        (block->rf).wav.magnitude = (short*)AllocNode(sizeof(short)*n);
        (block->rf).wav.phase = (short*)AllocNode(sizeof(short)*n);
        readshort((block->rf).wav.magnitude, n, fid);
        (block->rf).wav.magnitude[n-1]++;   /* EOS bit */
        readshort((block->rf).wav.phase, n, fid);
        /* TODO: turn on EOS bit for phase as well ? */

        fprintf(stderr,"\tread_block(): rf type, delay, amplitude, nSamples = %d, %d, %.5f, %d\n", 
            (block->rf).type, (block->rf).delay, (block->rf).amplitude, (block->rf).wav.nSamples);
    }

    /* gradients */
    read_grad(fid, &(block->gx));
    read_grad(fid, &(block->gy));
    read_grad(fid, &(block->gz));

    /* ADC */
    readshort(&((block->adc).type), 1, fid);
    if ((block->adc).type == ADC) {
        readshort(&((block->adc).numSamples), 1, fid);
        readshort(&((block->adc).dwell), 1, fid);  /* us */
        readshort(&((block->adc).delay), 1, fid);  /* us */
    }

    /* Print some values for debugging */
    fprintf(stderr, "\trf.type = %d\n", (block->rf).type);
    fprintf(stderr, "\tgx.type = %d\n", (block->gx).type);
    fprintf(stderr, "\tgy.type = %d\n", (block->gy).type);
    fprintf(stderr, "\tgz.type = %d\n", (block->gz).type);
    fprintf(stderr, "\tadc.type = %d\n", (block->adc).type);

    return JFN_SUCCESS;
}

/* read one gradient section from block file */
int read_grad(FILE *fid, PULSEQ_GRAD* g) {

    short amp_int16;
    short n;

    readshort(&(g->type), 1, fid);

    if ((g->type) != NULL) {
        readshort(&(g->delay), 1, fid);         /* us */

        readshort(&(amp_int16), 1, fid);
        g->amplitude = 1.0 * amp_int16/GSCALE;  /* Gauss/cm */
    }

    if ((g->type) == TRAP) {
        readshort(&((g->shape).trap.riseTime), 1, fid);  /* us */
        readshort(&((g->shape).trap.flatTime), 1, fid);  /* us */
        readshort(&((g->shape).trap.fallTime), 1, fid);  /* us */
    } else if ((g->type) == ARBITRARY) {
        readshort(&n, 1, fid);
        (g->shape).wav.nSamples = n;       
        (g->shape).wav.magnitude = (short*)AllocNode(sizeof(short)*n);
        readshort((g->shape).wav.magnitude, n, fid);
    }

    return JFN_SUCCESS;
}

