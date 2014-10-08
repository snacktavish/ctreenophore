/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *
 *  and
 *
 *  Programs of the PHYLIP package by Joe Felsenstein.
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *
 *  For any other enquiries send an Email to Alexandros Stamatakis
 *  Alexandros.Stamatakis@epfl.ch
 *
 *  When publishing work that is based on the results from RAxML-VI-HPC please cite:
 *
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models".
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

/** @file utils.c
 *  
 *  @brief Miscellaneous general utility and helper functions
 */
#ifdef WIN32
#include <direct.h>
#endif

#ifndef WIN32
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>
#endif

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>
#include <assert.h>
#include <errno.h>
#include "cycle.h"



#if ! (defined(__ppc) || defined(__powerpc__) || defined(PPC))
#if (defined(__AVX) || defined(__SSE3))
#include <xmmintrin.h>
#endif
/*
   special bug fix, enforces denormalized numbers to be flushed to zero,
   without this program is a tiny bit faster though.
#include <emmintrin.h> 
#define MM_DAZ_MASK    0x0040
#define MM_DAZ_ON    0x0040
#define MM_DAZ_OFF    0x0000
*/
#endif

#include "pll.h"

#define GLOBAL_VARIABLES_DEFINITION

#include "globalVariables.h"

static void pllTreeInitDefaults (pllInstance * tr, int tips);
static void initializePartitionsSequential(pllInstance *tr, partitionList *pr);


/** @ingroup instanceLinkingGroup
    @brief Initializes the PLL tree topology according to a parsed newick tree

    Set the tree topology based on a parsed and validated newick tree

    @param tree
      The PLL instance

    @param nt
      The \a pllNewickTree wrapper structure that contains the parsed newick tree

    @param useDefaultz
      If set to \b PLL_TRUE then the branch lengths will be reset to the default
      value.
*/
void
pllTreeInitTopologyNewick (pllInstance * tr, pllNewickTree * nt, int useDefaultz)
{
  pllStack * nodeStack = NULL;
  pllStack * head;
  struct item_t * item;
  int i, j, k;
  
/*
  for (i = 0; i < partitions->numberOfPartitions; ++ i)
   {
     partitions->partitionData[i] = (pInfo *) rax_malloc (sizeof (pInfo));
     partitions->partitionData[i]->partitionContribution = -1.0;
     partitions->partitionData[i]->partitionLH           =  0.0;
     partitions->partitionData[i]->fracchange            =  1.0;
   }
*/
  
  pllTreeInitDefaults (tr, nt->tips);

  i = nt->tips + 1;
  j = 1;
  nodeptr v;
  
  
  for (head = nt->tree; head; head = head->next)
  {
    item = (struct item_t *) head->item;
    if (!nodeStack)
     {
       pllStackPush (&nodeStack, tr->nodep[i]);
       pllStackPush (&nodeStack, tr->nodep[i]->next);
       pllStackPush (&nodeStack, tr->nodep[i]->next->next);
       ++i;
     }
    else
     {
       v = (nodeptr) pllStackPop (&nodeStack);
       if (item->rank)  /* internal node */
        {
          v->back           = tr->nodep[i];
          tr->nodep[i]->back = v; //t->nodep[v->number]
          pllStackPush (&nodeStack, tr->nodep[i]->next);
          pllStackPush (&nodeStack, tr->nodep[i]->next->next);
          double z = exp((-1 * atof(item->branch))/tr->fracchange);
          if(z < PLL_ZMIN) z = PLL_ZMIN;
          if(z > PLL_ZMAX) z = PLL_ZMAX;
          for (k = 0; k < PLL_NUM_BRANCHES; ++ k)
             v->z[k] = tr->nodep[i]->z[k] = z;

          ++ i;
        }
       else             /* leaf */
        {
          v->back           = tr->nodep[j];
          tr->nodep[j]->back = v; //t->nodep[v->number];

          double z = exp((-1 * atof(item->branch))/tr->fracchange);
          if(z < PLL_ZMIN) z = PLL_ZMIN;
          if(z > PLL_ZMAX) z = PLL_ZMAX;
          for (k = 0; k < PLL_NUM_BRANCHES; ++ k)
            v->z[k] = tr->nodep[j]->z[k] = z;
            
          //t->nameList[j] = strdup (item->name);
          tr->nameList[j] = (char *) rax_malloc ((strlen (item->name) + 1) * sizeof (char));
          strcpy (tr->nameList[j], item->name);
          
          pllHashAdd (tr->nameHash, tr->nameList[j], (void *) (tr->nodep[j]));
          ++ j;
        }
     }
  }
  
  tr->start = tr->nodep[1];
  
  pllStackClear (&nodeStack);

  if (useDefaultz == PLL_TRUE) 
    resetBranches (tr);
}

