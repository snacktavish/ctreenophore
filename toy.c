#include <stdio.h>
#include <pll/pll.h>
 
int main (int argc, char * argv[])
{
  pllInstance * inst;
  pllInstanceAttr attr;
  pllAlignmentData * alignmentData;
 
  /* set PLL instance attributes */
  attr.rateHetModel = PLL_GAMMA;
  attr.fastScaling  = PLL_FALSE;
  attr.saveMemory   = PLL_FALSE;
  attr.useRecom     = PLL_FALSE;
  /*attr.randomNumber = 0x12345;  no attribute...*/
  inst = pllCreateInstance (&attr);      /* Create the PLL instance */
  pllNewickTree * pllNewickParseString (const char *newick); 
 
  pllDestroyInstance (inst);             /* Destroy the PLL instance */
}

