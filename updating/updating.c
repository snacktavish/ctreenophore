#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <pll.h>
 
/*usage for time being is with  $updating schirrm_query.phy schirrm_query.tre schirrm.model*/
int main (int argc, char * argv[])
{
  pllInstance * inst;
  pllInstance * inst2; /* This is apparently dangerous because of global variables... BUT you can't have two trees in one instance either...*/
  pllInstanceAttr attr;
  pllAlignmentData * alignmentData;
  pllNewickTree * newick;
  partitionList * partitions;
  pllQueue * partitionInfo;
  nodeptr examplenode;
  double alpha = 0.75;
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

    printf ("Default settings with empirical frequencies\n"); 
    
    printf ("+-------------------------------------------------------+\n");
    printf ("  log-likelihood..: %f\n", inst->likelihood);
    printf ("+-------------------------------------------------------+\n\n\n");


  
  examplenode = inst->nodep[10];
  size_t states = (size_t)partitions->partitionData[0]->states;
  size_t width = (size_t) partitions->partitionData[0]->width;
  size_t categories = 4;
  printf("+-States is %lu-Width is %lu, categories is %lu,---numtips-is %i-------+\n\n\n",
    states, width, categories, inst->mxtips);
  double * outProbs;
  outProbs = malloc(width * categories * states * sizeof(double));
  int partition = 0;

  printf ("+-------------------------NEW------------------------------+\n\n\n");
  pllGetCLV(inst, partitions, examplenode, partition, outProbs);
  int i;
  for (i=0;i < 10;i++)
  {
    printf("%lf\n",outProbs[i]);
  }
  

  /* Perform some changes to the model */
  printf ("+-------------------------NEW------------------------------+\n\n\n");
  printf ("Setting Alpha to %f\n", alpha);
  pllSetFixedAlpha(alpha, 0, partitions, inst);
  printf ("Setting frequencies to:\n  A -> %f\n  C -> %f\n  G -> %f\n  T -> %f\n", freqs[0], freqs[1], freqs[2], freqs[3]);
  pllSetFixedBaseFrequencies(freqs, 4, 0, partitions, inst);
/* Evaluate again the likelihood */
 
  pllEvaluateLikelihood (inst, partitions, inst->start, PLL_TRUE, PLL_FALSE);
 
  printf ("+-------------------------------------------------------+\n");
  printf ("  log-likelihood..: %f\n", inst->likelihood);
  printf ("+-------------------------------------------------------+\n");
 
printf ("+-------------------------NEW------------------------------+\n\n\n");
  pllGetCLV(inst, partitions, examplenode, partition, outProbs);
  int j;
  for (j=0;j < 10;j++)
  {
    printf("%lf\n",outProbs[j]);
  }
  
printf ("+-------------------------NEW------------------------------+\n\n\n");

  /* free all allocated memory to eliminate memory leaks */
  pllPartitionsDestroy (inst, &partitions);
  pllDestroyInstance(inst);
  free(outProbs);
 
  return (EXIT_SUCCESS);
}