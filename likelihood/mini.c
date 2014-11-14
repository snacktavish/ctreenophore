#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <pll.h>

int main (int argc, char * argv[])
{
  pllNewickTree * newick;
  int val;
  pllInstance * tr;
  pllInstanceAttr attr;
  
  attr.rateHetModel     = 9;
  attr.fastScaling      = 9;
  attr.saveMemory       = 9;
  attr.useRecom         = 9;
  attr.randomNumberSeed = 0xDEADBEEF;
  attr.numberOfThreads  = 8; 

  tr = pllCreateInstance (&attr); 

  if (argc != 2)
   {
     fprintf (stderr, "usage: %s [newick-file]\n", argv[0]);
     return (EXIT_FAILURE);
   }
 
  /* Parse a NEWICK file */
  newick = pllNewickParseFile (argv[1]);

  if (!newick)
   {
     fprintf (stderr, "Error while parsing newick file %s\n", argv[1]);
     return (EXIT_FAILURE);
   }
  else
   {
    fprintf (stderr, "Succesfully parsed newick file %s\n", argv[1]);
   } 

  val=pllValidateNewick (newick);

  if (!pllValidateNewick (newick))
   {          
     fprintf (stderr, "Invalid phylogenetic tree\n");
     printf ("%d\n", errno);
     return (EXIT_FAILURE);
   }  
  fprintf (stderr, "Alloo %i\n", val);
  return (EXIT_SUCCESS);

  pllDestroyInstance (inst); 

}
