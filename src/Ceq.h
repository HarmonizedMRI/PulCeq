#ifndef CEQ_H
#define CEQ_H

#include "block.h"           /* PULSEQ_BLOCK struct */
#include "blockgroups.h"     /* GROUP struct */

/* struct containing entire sequence definition */
typedef struct {
   PULSEQ_BLOCK* parentBlocks;
   short nParentBlocks;    /* length of parentBlocks vector */

   BLOCKGROUPS groups;
   short nGroups;    

   short** dyn;      /* dynamic settings: waveform amplitudes, phase offsets, etc */
                     /* dimensions of dyn = nMax-by-9 */
   int nMax;         /* number of blocks (rows in BLOCKS section) in .seq file */
} CEQ;

int read_GE_seq(char* fname, CEQ* ceq);
int ceq_freemem(CEQ* ceq);

#endif
