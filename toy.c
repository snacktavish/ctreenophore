#include <stdio.h>
#include <pll/pll.h>
 
int main (int argc, char * argv[])
{
  pllInstance * pllCreateInstance (pllInstanceAttr *attr);

typedef struct
{  
  int rateHetModel;
  int fastScaling;
  int saveMemory;
  int useRecom;
  long randomNumberSeed;
  int numberOfThreads;
} pllInstanceAttr

}

