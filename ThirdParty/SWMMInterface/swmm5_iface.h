// swmm5_iface.h
//
// Header file for SWMM 5 interfacing functions
//
// #include this file in any C module that references the functions
// contained in swmm5_iface.c.
//

#ifdef __cplusplus
extern "C" {
#endif

extern int    SWMM_Nperiods;           // number of reporting periods
extern int    SWMM_FlowUnits;          // flow units code
extern int    SWMM_Nsubcatch;          // number of subcatchments
extern int    SWMM_Nnodes;             // number of drainage system nodes
extern int    SWMM_Nlinks;             // number of drainage system links
extern int    SWMM_Npolluts;           // number of pollutants tracked
extern double SWMM_StartDate;          // start date of simulation
extern int    SWMM_ReportStep;         // reporting time step (seconds)

int    RunSwmmExe(char* cmdLine);
int    RunSwmmDll(char* inpFile, char* rptFile, char* outFile);
int    OpenSwmmOutFile(char* outFile);
int    GetSwmmResult(int iType, int iIndex, int vIndex, int period, float* value);
void   CloseSwmmOutFile(void);

#ifdef __cplusplus
}
#endif
