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
#include "seq.h"
#include "block.h"
#include "cores.h"
#include "constants.h"
#include "modules.h"   /* readshort */

/* Read sequence from file. Cf. pulsegeq.readGEseq() in Matlab toolbox */
int read_GE_seq(char* fname, SEQ* seq) { 

    FILE *fid;
    int p, i, j;
    short n1, n2;

    fprintf(stderr, "read_GE_seq(): reading %s \n", fname);

    if ((fid = fopen (fname, "r")) == NULL)
        return(JFN_FAILURE);

    /* Read parent blocks */
    readshort(&(seq->nParentBlocks), 1, fid);
    fprintf(stderr, "\tnParentBlocks = %d\n", seq->nParentBlocks);

    seq->parentBlocks = (PULSEQ_BLOCK*)AllocNode(sizeof(PULSEQ_BLOCK)*(seq->nParentBlocks));

    for (p = 0; p < seq->nParentBlocks; p++)  {
        read_block(fid, &(seq->parentBlocks[p]));
        fprintf(stderr, "\tread_GE_seq; parent block %d: ID, duration = %d, %d\n", p, 
            (seq->parentBlocks[p]).ID, (seq->parentBlocks[p]).duration); 
    }

    /* Read group definitions */
    readshort(&(seq->nGroups), 1, fid);

    seq->groups = (GROUP*)AllocNode(sizeof(GROUP)*(seq->nGroups));

    for (i = 0; i < seq->nGroups; i++)  { 
        read_group(fid, &(seq->groups[i]));
    }

    /* Read dynamic scan information */
    /* [coreID parentBlockID rfScale rfPhsScale rfFreq gxScale gyScale gzScale recPhsScale] */
    readshort(&n1, 1, fid);
    readshort(&n2, 1, fid);
    seq->nMax = (int)n1*MAXIAMP + (int)n2;
    fprintf(stderr, "\tseq->nMax = %d\n", seq->nMax);
    seq->dyn = (short**)AllocNode(sizeof(short*)*(seq->nMax));
    for (i = 0; i < seq->nMax; i++)  {
        seq->dyn[i] = (short*)AllocNode(sizeof(short)*DYNAMICSVECTORLENGTH);
        readshort(seq->dyn[i], DYNAMICSVECTORLENGTH, fid);
        fprintf(stderr, "\tdyn[%d] = ", i);
        for (j = 0; j < DYNAMICSVECTORLENGTH; j++) {
            fprintf(stderr, "%d ", seq->dyn[i][j]);
        }
        fprintf(stderr, "\n");
    }

    fclose(fid);

    return JFN_SUCCESS;
}

/* TODO: free memory: parentBlocks, groups, dyn */
int seq_freemem(SEQ* seq) {

    int p, j;

    /*
    for (p = 0; i < seq->nParentBlocks; p++)  {
        FreeNode(seq->parentBlocks[p]);
    }

    FreeNode(groups->groupIDs);
    */

    return JFN_SUCCESS;
}
