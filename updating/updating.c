#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <pll.h>
 
/*usage for time being is with  $updating schirrm_query.phy schirrm_query.tre schirrm.model*/
/*PLAN:
Read in two or more trees. that differ only in palcement of query sequence.
Create virtual root at query tip
Calculate conditional likelihood vecotrs for nodes wich represent the same bipartitions
Compare changes in CLV across the two trees as distance from MRCA of query placements increases.
*/
void traverse(pllInstance * tr, nodeptr p);
char ** traverse_recursive(pllInstance * tr, 
                           nodeptr p, 
                           int * tip_count);

int main (int argc, char * argv[])
{
  pllInstance * inst;
  pllInstance * inst2; /* This is apparently dangerous because of global variables... BUT you can't have two trees in one instance either...*/
  pllInstanceAttr attr;
  pllAlignmentData * alignmentData;
  pllNewickTree * newick;
  partitionList * partitions;
  pllQueue * partitionInfo;
/*  double alpha = 0.75;*/
  double freqs[4] = { 0.25, 0.25, 0.25, 0.25 };

if (argc != 4)
   {
     fprintf (stderr, "usage: %s [phylip-file] [newick-file] [partition-file]\n", argv[0]);
     return (EXIT_FAILURE);
   }

  /* set PLL instance attributes */
  attr.rateHetModel = PLL_GAMMA;
  attr.fastScaling  = PLL_FALSE;
  attr.saveMemory   = PLL_FALSE;
  attr.useRecom     = PLL_FALSE;
  attr.randomNumberSeed = 0x12345;

  inst = pllCreateInstance (&attr);      /* Create the PLL instance */
  inst2 = pllCreateInstance (&attr);    
 

   /* Parse a PHYLIP file */
  alignmentData = pllParseAlignmentFile (PLL_FORMAT_PHYLIP, argv[1]);
 
   if (!alignmentData)
   {
     fprintf (stderr, "Error while parsing %s\n", argv[1]);
     return (EXIT_FAILURE);
   }
 

  /* Parse a NEWICK file */
  newick = pllNewickParseFile (argv[2]);
 
  if (!newick)
   {
     fprintf (stderr, "Error while parsing newick file %s\n", argv[2]);
     return (EXIT_FAILURE);
   }
 
  /* check whether the valid newick tree is also a tree that can be processed
     with our nodeptr structure, i.e if it is unrooted (ternary at root, binary
     at all other inner nodes */
  if (!pllValidateNewick (newick))     {
     fprintf (stderr, "Invalid phylogenetic tree\n");
     printf ("%d\n", errno);
     return (EXIT_FAILURE);
   }

  /* Parse the partitions file into a partition queue structure */
  partitionInfo = pllPartitionParse (argv[3]);
 
  /* Validate the partitions */
  if (!pllPartitionsValidate (partitionInfo, alignmentData))
   {
     fprintf (stderr, "Error: Partitions do not cover all sites\n");
     return (EXIT_FAILURE);
   }
 
  /* Commit the partitions and build a partitions structure */
  partitions = pllPartitionsCommit (partitionInfo, alignmentData);
 
  /* We don't need the the intermedia partition queue structure anymore */
  pllQueuePartitionsDestroy (&partitionInfo);

  /* Set the topology of the PLL tree from a parsed newick tree */
  pllTreeInitTopologyNewick (inst, newick, PLL_FALSE);
  pllTreeInitTopologyNewick (inst2, newick, PLL_FALSE);

  if (!pllLoadAlignment (inst, alignmentData, partitions))
   {
     fprintf (stderr, "Incompatible tree/alignment combination\n");
     return (EXIT_FAILURE);
   }
  pllLoadAlignment (inst2, alignmentData, partitions);
   /* Initialize the model. Note that this function will also perform a full
     tree traversal and evaluate the likelihood of the tree. Therefore, you
     have the guarantee that tr->likelihood is the valid likelihood */
  pllInitModel(inst, partitions);
  pllInitModel(inst2, partitions);
 /* we can actually free the parsed alignmentData, since they
   are deep-copied in the instance. We can also free the
   parsed tree as it is also copied to the instance */
  pllAlignmentDataDestroy (alignmentData);
  pllNewickParseDestroy (&newick);
  

  printf ("Equal base frequencies\n"); 
  printf ("Setting frequencies to:\n  A -> %f\n  C -> %f\n  G -> %f\n  T -> %f\n", freqs[0], freqs[1], freqs[2], freqs[3]);
  pllSetFixedBaseFrequencies(freqs, 4, 0, partitions, inst);

  printf ("+-------------------------------------------------------+\n");
  printf ("  log-likelihood..: %f\n", inst->likelihood);
  printf ("+-------------------------------------------------------+\n\n\n");

  pllEvaluateLikelihood (inst, partitions, inst->start, PLL_TRUE, PLL_FALSE);

  printf ("+-------------------------------------------------------+\n");
  printf ("  log-likelihood..: %f\n", inst->likelihood);
  printf ("+-------------------------------------------------------+\n");
  size_t states = (size_t) partitions->partitionData[0]->states;
  size_t width = (size_t) partitions->partitionData[0]->width; 
  size_t categories = 4;  /*I though t categories was gamma rate categories, but steting it to less than 4 causes a seq fault....*/
  printf("+-States is %lu-Width is %lu, categories is %lu,---numtips-is %i-------+\n\n\n",
  states, width, categories, inst->mxtips);
  int partition = 0;
  double * outProbs;
  
  outProbs = malloc(width * categories * states * sizeof(double));
  int numvals = 1;
  int np;
  for (np=6; np<7; np++)
  {
      printf ("+---------------------node by node : NODE %i-, %i-----------------------------+\n\n\n", np, inst->nodep[np]->number);
      pllEvaluateLikelihood (inst, partitions, inst->start, PLL_TRUE, PLL_FALSE);

      printf ("+-------------------------------------------------------+\n");
      printf ("  log-likelihood..: %f\n", inst->likelihood);
      pllGetCLV(inst, partitions, inst->nodep[np], partition, outProbs);
      int j;
      for (j=0;j < numvals;j++)
          {
            printf("%lf\n",outProbs[j]);
          }
      pllEvaluateLikelihood (inst, partitions, inst->nodep[2], PLL_TRUE, PLL_FALSE);
      printf ("+----------------Root at tip 2-------------------------------------+\n");
      printf ("  log-likelihood..: %f\n", inst->likelihood);
      pllGetCLV(inst, partitions, inst->nodep[np], partition, outProbs);
      for (j=0;j < numvals;j++)
          {
            printf("%lf\n",outProbs[j]);
          }

      printf ("+----------------Root at tip 3----------------------------------+\n");
      pllEvaluateLikelihood (inst, partitions, inst->nodep[3], PLL_TRUE, PLL_FALSE);
      printf ("  log-likelihood..: %f\n", inst->likelihood);
      pllGetCLV(inst, partitions, inst->nodep[np]->next, partition, outProbs);
      for (j=0;j < numvals;j++)
          {
            printf("%lf\n",outProbs[j]);
          }

      printf ("+-----------------------------------------------------+\n\n\n");
  }

/*Find bipart for node*/

    printf ("+----------------Find biparts for node %i----------------------------\n", np);
    traverse(inst, inst->nodep[3]);   


  /* free all allocated memory to eliminate memory leaks */
  pllPartitionsDestroy (inst, &partitions);
  pllDestroyInstance(inst);
  free(outProbs);

  return (EXIT_SUCCESS);
}

char ** traverse_recursive(pllInstance * tr, 
                           nodeptr p, 
                           int * tip_count)
{
  int i;

  /* lists (and counts) of tip names for left and right subtree */
  char ** ltips;
  char ** rtips;
  int ltip_count = 0;
  int rtip_count = 0;

  /* tip list for current inner node, which will be formed as the concatenation
     of left and right subtree tiplists */
  char ** tiplist;

  /* if it's a tip, create a list of size 1, add a pointer to its label, 
     set tip_count to 1 and return */
  if (isTip(p->number, tr->mxtips)) 
  {
    tiplist = (char **)malloc(sizeof(char *));
    *tiplist = tr->tipNames[p->number];
    *tip_count = 1;
    return tiplist;
  }

  /* if it's an inner node, recursively go over its two subtrees */
  ltips = traverse_recursive(tr, p->next->back, &ltip_count);
  rtips = traverse_recursive(tr, p->next->next->back, &rtip_count);

  /* unify the two tip lists into one  */
  tiplist = (char **)malloc((ltip_count+rtip_count) * sizeof(char *));
  for (i = 0; i < ltip_count; ++i)
    tiplist[i] = ltips[i];
  for (i = ltip_count; i < ltip_count+rtip_count; ++i)
    tiplist[i] = rtips[i-ltip_count];
  *tip_count = ltip_count + rtip_count;

  /* free the two lists as we will not need them any more */
  free(ltips);
  free(rtips);

  printf ("mRCA of %i: \n", p->number);
  for (i = 0; i < *tip_count; ++i)
    printf("%s ", tiplist[i]);
  printf("\n");

  /* partially evaluate the likelihood here, for node p, don't remember the
     function name/arguments and output the CLV */

  return tiplist;
}

/* start of recursive traversal, make sure p is an inner node */
void traverse(pllInstance * tr, nodeptr p)
{
  int tip_count = 0;
  char ** tiplist = NULL;

  /* traverse first direction */
  tiplist = traverse_recursive(tr, p->back, &tip_count);
  free(tiplist);

  /* reset number of tip labels in list and traverse second direction */
  tiplist = traverse_recursive(tr, p->next->back, &tip_count);
  free(tiplist);

  /* reset number of tip labels in list and traverse third direction */
  traverse_recursive(tr, p->next->next->back, &tip_count);
  free(tiplist);
}

/* function returning all tips on one side of a root at a node */
/*int tip(pllInstance *inst) 
{ 

    nodeptr node1=inst->nodep[6];
    nodeptr node2=node1->next;
    nodeptr node3=node2->next;
    if ( node1->x == 1 ){
          printf( "XFLAG" );
    } else {
            if (node1->back->number < inst->mxtips) 
              return node1->back->number < inst->mxtips;
    }
    if ( node2->x == 1 ){
          printf( "XFLAG 2nd try" );
    } else {
            if (node2->back->number < inst->mxtips) 
              return node2->back->number < inst->mxtips;
    }
    if ( node3->x == 1 ){
        printf( "XFLAG 3nd try" );
        } else {
            if (node3->back->number < inst->mxtips) 
              return node3->back->number < inst->mxtips;
    }
    

      printf("\n%i\n",node1->number);
      printf("\n%i\n",inst->mxtips);

      printf ("+-----------------------------------------------------+\n\n\n");
  return 0;
}
*/