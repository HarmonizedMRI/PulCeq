#ifndef SEQ_H
#define SEQ_H

#include "block.h"
#include "blockgroups.h"

/* struct containing entire sequence definition */
typedef struct {
   PulseqBlock* parentBlocks;
   short nParentBlocks; 

   BlockGroups groups;
   int nGroups;    

   float** dyn;      /* dynamic settings: waveform amplitudes, phase offsets, etc */
                     /* dimensions of dyn = nMax-by-9 */

   int nMax;         /* number of blocks (rows in BLOCKS section) in .seq file */
} Seq;

#endif
