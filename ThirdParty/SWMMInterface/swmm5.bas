' Declarations of imported procedures from the EPASWMM DLL engine (SWMM5.DLL)

Declare Function swmm_run Lib "swmm5.dll" (ByVal F1 As String, ByVal F2 As String, ByVal F3 As String) As Long
Declare Function swmm_open Lib "swmm5.dll" (ByVal F1 As String, ByVal F2 As String, ByVal F3 As String) As Long
Declare Function swmm_start Lib "swmm5.dll" (ByVal saveFlag As Long) As Long
Declare Function swmm_step Lib "swmm5.dll" (elapsedTime As Double) As Long
Declare Function swmm_end Lib "swmm5.dll" () As Long
Declare Function swmm_report Lib "swmm5.dll" () As Long
Declare Function swmm_getMassBalErr Lib "swmm5.dll" (runoffErr As Single, flowErr As Single, qualErr As Single) As Long
Declare Function swmm_close Lib "swmm5.dll" () As Long
Declare Function swmm_getVersion Lib "swmm5.dll" () As Long
