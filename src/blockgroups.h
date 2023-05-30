#ifndef BLOCKGROUPS_H
#define BLOCKGROUPS_H

/* Struct containing block IDs for all block groups */
typedef struct {
    int   nGroups;           /* number of groups */
    int*  nBlocksInGroup;    /* number of blocks in each group, i.e.: */
                             /* size(blockIds[iGroup], 2) = nBlocksInGroup[iGroup] */
    int** blockIDs;          /* block id's for all groups. size(blockIds,1) = nGroups */
} BLOCKGROUPS;

int readBlockGroups(char *fname, BLOCKGROUPS* blockGroups);
int getDefaultBlockGroups(BLOCKGROUPS* blockGroups, int nBlocks);
void blockgroups_freemem(BLOCKGROUPS* blockGroups);


#endif
