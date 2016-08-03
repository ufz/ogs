// Test of the SWMM5 C Interfacing functions

// This is a command line executable that takes the name
// of a SWMM input file as its only command line argument
// and produces a time series listing to the console of the
// following system output results: total rainfall, total
// runoff and total outfall flow.

// This file must be compiled along with swmm5_iface.c
// and swmm5.h, and linked with swmm5.dll and swmm5.lib.

#include <stdio.h>
#include "swmm5_iface.h"

int main(int argc, char *argv[])
{
   int i, r;
   float x, y, z;
   char rptfile[] = "tmp.rpt";
   char outfile[] = "tmp.out";

// Check if a SWMM input file name is provided
   if (argc < 2)
   {
      printf("\nNo file name was provided.\n");
      return 0;
   }

// Run the SWMM analysis
   r = RunSwmmDll(argv[1], rptfile, outfile);
   if (r > 0)
   {
      printf("\nSWMM run was unsuccessful; error code = %d\n", r);
   }
   else
   {
   // Open outfile as a SWMM output file
      r = OpenSwmmOutFile(outfile);
      if (r == 1)
      {
         printf("\nInvalid results in SWMM output file.\n");
      }
      else if (r == 2)
      {
         printf("\nFile is not a SWMM output file.\n");
      }
      else
      {
         printf("\nTime       Total     Total     Total");
         printf("\nPeriod  Rainfall    Runoff   Outflow");
         printf("\n====================================");
         for (i=1; i<=SWMM_Nperiods; i++)
         {
             GetSwmmResult(3, 0, 1, i, &x);
             GetSwmmResult(3, 0, 4, i, &y);
             GetSwmmResult(3, 0, 11, i, &z);
             printf("\n%6d  %8.2f  %8.2f  %8.2f", i, x, y, z);
         }
         CloseSwmmOutFile();
      }
   }
   remove(rptfile);
   remove(outfile);
   return 0;
}