// swmm5_iface.pas
//
// Example code for interfacing SWMM 5 with Delphi Pascal programs.
//
// Remember to add this unit to the Uses clause of the calling program.

Unit SWMM5_IFACE

Interface

Uses SysUtils, Classes, Consts, WinTypes, WinProcs;

Var
  SWMM_Nperiods: Integer;              // number of reporting periods
  SWMM_FlowUnits: Integer;             // flow units code
  SWMM_Nsubcatch: Integer;             // number of subcatchments
  SWMM_Nnodes: Integer;                // number of drainage system nodes
  SWMM_Nlinks: Integer;                // number of drainage system links
  SWMM_Npolluts: Integer;              // number of pollutants tracked
  SWMM_StartDate: Double;              // start date of simulation
  SWMM_ReportStep: Integer;            // reporting time step (seconds)

  Function  RunSwmmExe(cmdLine : String): Integer;
  Function  RunSwmmDll(inpFile: String; rptFile: String; outFile: String): Integer;
  Function  OpenSwmmOutFile(outFile: String): Integer;
  Function  GetSwmmResult(iType, iIndex, vIndex, period: Integer; var value: Single):
            Integer;
  Procedure CloseSwmmOutFile;

Implementation

Uses swmm5;

Const
 SUBCATCH = 0;
 NODE     = 1;
 LINK     = 2;
 SYS      = 3;
 RECORDSIZE = 4;                       // number of bytes per file record

Var
  SubcatchVars: Integer;               // number of subcatch reporting variable
  NodeVars: Integer;                   // number of node reporting variables
  LinkVars: Integer;                   // number of link reporting variables
  SysVars: Integer;                    // number of system reporting variables
  Fout: File;                          // file handle
  StartPos: LongInt;                   // file position where results start
  BytesPerPeriod: LongInt;             // bytes used for results in each period


//-----------------------------------------------------------------------------
Function RunSwmmExe(cmdLine : String): Integer;
//-----------------------------------------------------------------------------
Var
  exitCode: Integer;
  pi: TProcessInformation;
  si: TStartupInfo;
Begin
  // --- Initialize data structures
  FillMemory(@si, sizeof(si), 0);
  FillMemory(@pi, sizeof(pi), 0);
  si.cb := sizeof(si);
  si.wShowWindow := SW_SHOWNORMAL;

  // --- launch swmm5.exe
  exitCode := CreateProcess(Nil, PChar(cmdLine), Nil,
			Nil, False, NORMAL_PRIORITY_CLASS,
			Nil, Nil, si, pi );

  // --- wait for program to end
  exitCode := WaitForSingleObject(pi.hProcess, INFINITE);

  // --- retrieve the error code produced by the program
  GetExitCodeProcess(pi.hProcess, exitCode);

  // --- release handles
  CloseHandle( pi.hProcess );
  CloseHandle( pi.hThread );
  Result := exitCode;
End;


//-----------------------------------------------------------------------------
Function RunSwmmDll(inpFile: String; rptFile: String; outFile: String): Integer;
//-----------------------------------------------------------------------------
Var
  err : Integer;
  elapsedTime : Double;
Begin

  // --- open a SWMM project
  err := swmm_open(inpFile, rptFile, outFile);
  If (err = 0) Then
  Begin

    // --- initialize all processing systems
    err := swmm_start(1);
    If (err = 0) Then
    Begin

      // --- step through the simulation
      Do
        ProcessMessages;
        err := swmm_step(elapsedTime);

        /////////////////////////////////////////////
        // --- call progress reporting function here,
        //     using elapsedTime as an argument
	/////////////////////////////////////////////

      While (elapsedTime > 0.0) And (err = 0);
    End;

    // --- close all processing systems
    swmm_end();
  End;

  // --- close the project
  swmm_close();
  Result := err;
End;


//-----------------------------------------------------------------------------
Function OpenSwmmOutFile(outFile: String): Integer;
//-----------------------------------------------------------------------------
Var
  magic1, magic2, errCode, version, position: Integer;
  offset, offset0: LongInt;
  err: Integer;
Begin
  // --- open the output file
  AssignFile(Fout, outFile);
  {$I-}
  Reset(Fout, 1);
  {$I+}
  if IOResult <> 0 Then
  Begin
    Result := 2;
    Exit;
  End;

  // --- check that file contains at least 14 records
  If FileSize(Fout) < 14*RECORDSIZE Then
  Begin
    Result := 1;
    CloseFile(Fout);
    Exit;
  End;

  // --- read parameters from end of file
  Seek(Fout, FileSize(Fout) -5*RECORDSIZE);
  BlockRead(Fout, position, RECORDSIZE);
  offset0 := position;
  BlockRead(Fout, position, RECORDSIZE);
  StartPos := position;
  BlockRead(Fout, SWMM_Nperiods, RECORDSIZE);
  BlockRead(Fout, errCode, RECORDSIZE);
  BlockRead(Fout, magic2, RECORDSIZE);

  // --- read magic number from beginning of file
  Seek(Fout, 0);
  BlockRead(Fout, magic1, RECORDSIZE);

  // --- perform error checks
  If magic1 <> magic2 Then err = 1
  Else If errCode <> 0 Then err = 1
  Else If SWMM_Nperiods = 0 Then err = 1
  Else err = 0;

  // --- quit if errors found
  if (err > 0 )
  {
    Result := err;
    CloseFile(Fout);
    Exit;
  }

  // --- otherwise read additional parameters from start of file
  BlockRead(Fout, version, RECORDSIZE); 
  BlockRead(Fout, SWMM_FlowUnits, RECORDSIZE);
  BlockRead(Fout, SWMM_Nsubcatch, RECORDSIZE);
  BlockRead(Fout, SWMM_Nnodes, RECORDSIZE);
  BlockRead(Fout, SWMM_Nlinks, RECORDSIZE);
  BlockRead(Fout, SWMM_Npolluts, RECORDSIZE);

  // Skip over saved subcatch/node/link input values
  offset := (SWMM_Nsubcatch+2) * RECORDSIZE     // Subcatchment area
             + (3*SWMM_Nnodes+4) * RECORDSIZE  // Node type, invert & max depth
             + (5*SWMM_Nlinks+6) * RECORDSIZE; // Link type, z1, z2, max depth & length
  offset := offset0 + offset;
  Seek(Fout, offset);

  // Read number & codes of computed variables
  BlockRead(Fout, SubcatchVars, RECORDSIZE); // # Subcatch variables
  Seek(Fout, FilePos(Fout)+(SubcatchVars*RECORDSIZE));
  BlockRead(Fout, NodeVars, RECORDSIZE);     // # Node variables
  Seek(Fout, FilePos(Fout)+(NodeVars*RECORDSIZE));
  BlockRead(Fout, LinkVars, RECORDSIZE);     // # Link variables
  Seek(Fout, FilePos(Fout)+(LinkVars*RECORDSIZE));
  BlockRread(Fout, SysVars, RECORDSIZE);     // # System variables

  // --- read data just before start of output results
  offset := StartPos - 3*RECORDSIZE;
  Seek(Fout, offset);
  BlockRead(Fout, SWMM_StartDate, 2*RECORDSIZE);
  BlockRead(Fout, SWMM_ReportStep, RECORDSIZE);

  // --- compute number of bytes of results values used per time period
  BytesPerPeriod := 2*RECORDSIZE +      // date value (a double)
                    (SWMM_Nsubcatch*SubcatchVars +
                     SWMM_Nnodes*NodeVars +
                     SWMM_Nlinks*LinkVars +
                     SysVars)*RECORDSIZE;

  // --- return with file left open
  Result := err;
End;


//-----------------------------------------------------------------------------
Function GetSwmmResult(iType, iIndex, vIndex, period: Integer; var value: Single):
         Integer;
//-----------------------------------------------------------------------------
Var
  offset: LongInt;
Begin
  // --- compute offset into output file
  Result := 0;
  value := 0.0;
  offset := StartPos + (period-1)*BytesPerPeriod + 2*RECORDSIZE;
  if ( iType = SUBCATCH ) Then
    offset := offset + RECORDSIZE*(iIndex*SubcatchVars + vIndex)
  else if (iType = NODE) Then
    offset := offset +  RECORDSIZE*(SWMM_Nsubcatch*SubcatchVars +
                                    iIndex*NodeVars + vIndex)
  else if (iType = LINK) Then
    offset := offset + RECORDSIZE*(SWMM_Nsubcatch*SubcatchVars +
                                   SWMM_Nnodes*NodeVars +
                                   iIndex*LinkVars + vIndex)
  else if (iType == SYS) Then
    offset := offset + RECORDSIZE*(SWMM_Nsubcatch*SubcatchVars +
                                   SWMM_Nnodes*NodeVars +
                                   SWMM_Nlinks*LinkVars + vIndex)
  else Exit;

  // --- re-position the file and read the result
  Seek(Fout, offset);
  BlockRead(Fout, value, RECORDSIZE);
  Result := 1;
End;


//-----------------------------------------------------------------------------
Procedure CloseOutFile;
//-----------------------------------------------------------------------------
Begin
  CloseFile(Fout);
End;
