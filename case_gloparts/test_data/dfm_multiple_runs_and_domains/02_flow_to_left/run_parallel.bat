@ echo off

    rem When using mpich2 for the first time on a machine:
    rem Execute "smpd -install" as administrator:
    rem     Preparation: Check that your Delft3D installation contains "...\x64\share\bin\smpd.exe". Optionally copy it to a local directory (it will run as a service).
    rem     "Start" -> "All programs" -> "Accessories", right-click "Command Prompt", "Run as Administrator"
    rem     In this command box:
    rem         cd ...\x64\share\bin
    rem (i.e.   cd "c:\Program Files (x86)\Deltares\Delft3D Flexible Mesh Suite HMWQ (2021.04)\plugins\DeltaShell.Dimr\kernels\x64\share\bin\"   )
    rem         smpd -install
    rem     When there is an smpd already running on the machine, it must be ended first, using the Microsoft Task Manager, 
    rem     or in the command  box: smpd -uninstall

call "p:\d-hydro\dimrset\latest\x64\dflowfm\scripts\run_dflowfm.bat" "--partition:ndomains=2:icgsolver=6" cb_2d.mdu

call "p:\d-hydro\dimrset\latest\x64\dimr\scripts\run_dimr_parallel.bat" 2 dimr_config.xml

    rem To prevent the DOS box from disappearing immediately: remove the rem on the following line
pause
