unit swmm5;

{ Declarations of imported procedures from the EPASWMM DLL engine }
{ (SWMM5.DLL) }

interface

function   swmm_run(F1: PChar; F2: PChar; F3: PChar): Integer; stdcall;
function   swmm_open(F1: PChar; F2: PChar; F3: PChar): Integer; stdcall;
function   swmm_start(SaveFlag: Integer): Integer; stdcall;
function   swmm_step(var ElapsedTime: Double): Integer; stdcall;
function   swmm_end: Integer; stdcall;
function   swmm_getMassBalErr(var Erunoff: Single; var Eflow: Single;
               var Equal: Single): Integer; stdcall;
function   swmm_close: Integer; stdcall;

implementation

function   swmm_run;    external 'SWMM5.DLL';
function   swmm_open;   external 'SWMM5.DLL';
function   swmm_start;  external 'SWMM5.DLL';
function   swmm_step;   external 'SWMM5.DLL';
function   swmm_end;    external 'SWMM5.DLL';
function   swmm_getMassBalErr; external 'SWMM5.DLL';
function   swmm_close;  external 'SWMM5.DLL';
end.

