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

#include "tp_defs.h"   
#include "modgroups.h"

/* Read cores.txt  */
int readModGroups(char* fname, TP_MODGROUPS* modGroups) {

    FILE *fid;
    int i, j;
    char line[200];

    fprintf(stderr,"readcorefile(): reading %s \n", fname);

    if ((fid = fopen (fname, "r")) == NULL)
        return(TPFAILURE);

    fgets(line, 200, fid);   /* skip line */
    fscanf(fid, "%d\n", &(modGroups->nGroups));
    fgets(line, 200, fid);   /* skip line */

    fprintf(stderr, "nGroups = %d\n", modGroups->nGroups);

    modGroups->nModulesInGroup = (int*)AllocNode(sizeof(int)*(modGroups->nGroups));
    modGroups->modIDs = (int**) AllocNode(sizeof(int*)*modGroups->nGroups);

    /* Step through rows in modgroups.txt and get module instance id's */
    for (i = 0; i < modGroups->nGroups; i++)  { 
        fscanf(fid, "%d", &(modGroups->nModulesInGroup[i]));
        /* fprintf(stderr, "nModulesInGroup[%d] = %d\n", i, modGroups->nModulesInGroup[i]); */
        modGroups->modIDs[i] = (int*)AllocNode(sizeof(int)*modGroups->nModulesInGroup[i]);
        for (j = 0; j < modGroups->nModulesInGroup[i]; j++)  {
            fscanf(fid, "%d", &(modGroups->modIDs[i][j]));
            fprintf(stderr, "modIDs[%d][%d] = %d\n", i, j, modGroups->modIDs[i][j]);
        }
    }

    fclose (fid);
    return TPSUCCESS;
}

/* Each module becomes its own group */
int getDefaultModGroups(TP_MODGROUPS* modGroups, int nModules) {

    int i, j;

    modGroups->nGroups = nModules;

    modGroups->nModulesInGroup = (int*)AllocNode(sizeof(int)*(modGroups->nGroups));
    modGroups->modIDs = (int**) AllocNode(sizeof(int*)*modGroups->nGroups);

    for (i = 0; i < modGroups->nGroups; i++)  { 
        modGroups->nModulesInGroup[i] = 1;
        modGroups->modIDs[i] = (int*)AllocNode(sizeof(int)*modGroups->nModulesInGroup[i]);
        for (j = 0; j < modGroups->nModulesInGroup[i]; j++)  {
            modGroups->modIDs[i][j] = i;
            fprintf(stderr, "modIDs[%d][%d] = %d\n", i, j, modGroups->modIDs[i][j]);
        }
    }

    return TPSUCCESS;
}

void modgroups_freemem(TP_MODGROUPS* modGroups) {

    int i;

    for (i = 0; i < modGroups->nGroups; i++)  {
        FreeNode(modGroups->modIDs[i]);
    }

    FreeNode(modGroups->modIDs);
    FreeNode(modGroups->nModulesInGroup);
}
