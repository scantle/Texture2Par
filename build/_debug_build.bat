set FCFLAGS=-O0 -g -cpp
set FCFLAGS=%FCFLAGS% -fbacktrace
set FCFLAGS=%FCFLAGS% -ffree-line-length-none
set FCFLAGS=%FCFLAGS% -ffixed-line-length-none
set FCFLAGS=%FCFLAGS% -fdec
set FCFLAGS=%FCFLAGS% -fallow-argument-mismatch
rem set FCFLAGS=%FCFLAGS% -fopenmp
set FCFLAGS=%FCFLAGS% -fcheck=all
set FCFLAGS=%FCFLAGS% -finit-real=snan
set FCFLAGS=%FCFLAGS% -finit-integer=-999999
set FLFLAGS=-fopenmp -g -Wl,--disable-dynamicbase

make PROGRAM=../bin/Texture2Par_dbg_nodyn FCFLAGS="%FCFLAGS%"
rem FCFLAGS="-O0 -g -cpp -fbacktrace -ffree-line-length-none -ffixed-line-length-none -fdec -fallow-argument-mismatch -fopenmp -fcheck=all -finit-real=snan -finit-integer=-999999" FLFLAGS="-fopenmp -g -Wl,--disable-dynamicbase
rem make -C build -B PROGRAM=../bin/Texture2Par_dbg_check FCFLAGS="-g -O0 -cpp -fbacktrace -fcheck=all -ffree-line-length-none -ffixed-line-length-none -fdec -fallow-argument-mismatch -fopenmp" FLFLAGS="-g -fopenmp