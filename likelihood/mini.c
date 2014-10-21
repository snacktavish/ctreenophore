#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <pll.h>

int main (int argc, char * argv[])
{
  pllNewickTree * newick;
  
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
  
  return (EXIT_SUCCESS);
}
