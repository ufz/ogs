

#include <string>
#include <iostream>
#include "ThirdParty/SWMMInterface/swmm5_iface.h"

int main(int argc, char *argv[])
{
   // Open outfile as a SWMM output file
    std::string outfile = argv[1];
    int r = OpenSwmmOutFile(const_cast<char*>(outfile.c_str()));
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
        for (int i=1; i<=SWMM_Nperiods; i++)
        {
            float x, y, z;
            GetSwmmResult(3, 0, 1, i, &x);
            GetSwmmResult(3, 0, 4, i, &y);
            GetSwmmResult(3, 0, 11, i, &z);
            printf("\n%6d  %8.2f  %8.2f  %8.2f", i, x, y, z);
        }
        CloseSwmmOutFile();
    }
   return 0;
}