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
  attr.randomNumberSeed = 0x12345;
  inst = pllCreateInstance (&attr);      /* Create the PLL instance */
 
 
  pllDestroyInstance (inst);             /* Destroy the PLL instance */
}

