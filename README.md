# pc4sbml

Demonstrate the use of SBML metabolic models in PhysiCell (using libRoadrunner). We primarily target the app running on nanoHUB. See below for building on MinGW.

```
from the root dir:
$ cd src
$ make  # assumes you've downloaded libRoadrunner binary

To just test PhysiCell (without the Jupyter GUI):
$ myproj config.xml
```

To use the Jupyter notebook GUI:
```
$ cp myproj ../bin/

If you have modified `data/PhysiCell_settings.xml`, re-create the (two) Python modules used in the GUI:

#  From the root dir (e.g., ~/git/sbmlsim)
$ cd data
$ python xml2jupyter.py PhysiCell_settings.xml
$ cp user_params.py ../bin
$ cp microenv_params.py ../bin

# And you may need to update the initial.xml file if the (number,name, order of) substrates changes:

~/git/physicell_libroadrunner/src/output$ cp initial.xml ../../data/initial.xml 

# Install some additional modules and run the GUI:
pip install -U display_xml    # rf.  https://github.com/mpacer/display_xml
pip install -U hublib   # doesn't work on Windows it seems

~/git/physicell_libroadrunner$ jupyter notebook sbmlsim.ipynb
```

In this version, we just download the binary libRoadrunner libs/headers for
the relevant OS.

## Build on MinGW

My test consisted of downloading libRoadrunner binaries from https://sourceforge.net/projects/libroadrunner/files/libroadrunner-1.4.18/

specifically, I got:
https://sourceforge.net/projects/libroadrunner/files/libroadrunner-1.4.18/roadrunner-win64-vs14-cp35m.zip/download

I unzipped them into my `~/dev` directory. 

Then I attempted various incantations in `Make-mingw` to build. So far, I cannot get the linker to find the libRR functions. I've tried to run `make -f Make-mingw` in both a `Command Prompt` shell and a `MINGW64/Git Bash` shell. And I copied the `roadrunner_c_api.lib` to my local working `/src` dir, to try to eliminate problems. No luck so far. Here's an ugly verbose build/failure:
```
heiland@LAPTOP-P2P2QC36 MINGW64 ~/git/pc4sbml/src (master)
$ make -f Make-mingw
g++ -march=native  -fomit-frame-pointer -fopenmp -m64 -std=c++11 -I/c/Users/heiland/dev/roadrunner-win64-vs14-cp35m/include/rr/C  -v -lroadrunner_c_api  -o myproj BioFVM_vector.o BioFVM_mesh.o BioFVM_microenvironment.o BioFVM_solvers.o BioFVM_matlab.o BioFVM_utilities.o BioFVM_basic_agent.o BioFVM_MultiCellDS.o BioFVM_agent_container.o  pugixml.o PhysiCell_phenotype.o PhysiCell_cell_container.o PhysiCell_standard_models.o PhysiCell_cell.o PhysiCell_custom.o PhysiCell_utilities.o PhysiCell_SVG.o PhysiCell_pathology.o PhysiCell_MultiCellDS.o PhysiCell_various_outputs.o PhysiCell_pugixml.o PhysiCell_settings.o  heterogeneity.o main.cpp
Using built-in specs.
COLLECT_GCC=c:\Program Files\mingw-w64\x86_64-8.1.0-posix-seh-rt_v6-rev0\mingw64\bin\g++.exe
COLLECT_LTO_WRAPPER=c:/Program\ Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../libexec/gcc/x86_64-w64-mingw32/8.1.0/lto-wrapper.exe
Target: x86_64-w64-mingw32
Configured with: ../../../src/gcc-8.1.0/configure --host=x86_64-w64-mingw32 --build=x86_64-w64-mingw32 --target=x86_64-w64-mingw32 --prefix=/mingw64 --with-sysroot=/c/mingw810/x86_64-810-posix-seh-rt_v6-rev0/mingw64 --enable-shared --enable-static --disable-multilib --enable-languages=c,c++,fortran,lto --enable-libstdcxx-time=yes --enable-threads=posix --enable-libgomp --enable-libatomic --enable-lto --enable-graphite --enable-checking=release --enable-fully-dynamic-string --enable-version-specific-runtime-libs --disable-libstdcxx-pch --disable-libstdcxx-debug --enable-bootstrap --disable-rpath --disable-win32-registry --disable-nls --disable-werror --disable-symvers --with-gnu-as --with-gnu-ld --with-arch=nocona --with-tune=core2 --with-libiconv --with-system-zlib --with-gmp=/c/mingw810/prerequisites/x86_64-w64-mingw32-static --with-mpfr=/c/mingw810/prerequisites/x86_64-w64-mingw32-static --with-mpc=/c/mingw810/prerequisites/x86_64-w64-mingw32-static --with-isl=/c/mingw810/prerequisites/x86_64-w64-mingw32-static --with-pkgversion='x86_64-posix-seh-rev0, Built by MinGW-W64 project' --with-bugurl=https://sourceforge.net/projects/mingw-w64 CFLAGS='-O2 -pipe -fno-ident -I/c/mingw810/x86_64-810-posix-seh-rt_v6-rev0/mingw64/opt/include -I/c/mingw810/prerequisites/x86_64-zlib-static/include -I/c/mingw810/prerequisites/x86_64-w64-mingw32-static/include' CXXFLAGS='-O2 -pipe -fno-ident -I/c/mingw810/x86_64-810-posix-seh-rt_v6-rev0/mingw64/opt/include -I/c/mingw810/prerequisites/x86_64-zlib-static/include -I/c/mingw810/prerequisites/x86_64-w64-mingw32-static/include' CPPFLAGS=' -I/c/mingw810/x86_64-810-posix-seh-rt_v6-rev0/mingw64/opt/include -I/c/mingw810/prerequisites/x86_64-zlib-static/include -I/c/mingw810/prerequisites/x86_64-w64-mingw32-static/include' LDFLAGS='-pipe -fno-ident -L/c/mingw810/x86_64-810-posix-seh-rt_v6-rev0/mingw64/opt/lib -L/c/mingw810/prerequisites/x86_64-zlib-static/lib -L/c/mingw810/prerequisites/x86_64-w64-mingw32-static/lib '
Thread model: posix
gcc version 8.1.0 (x86_64-posix-seh-rev0, Built by MinGW-W64 project)
COLLECT_GCC_OPTIONS='-march=native' '-fomit-frame-pointer' '-fopenmp' '-m64' '-std=c++11' '-I' 'c:/Users/heiland/dev/roadrunner-win64-vs14-cp35m/include/rr/C' '-v' '-o' 'myproj.exe' '-shared-libgcc' '-mthreads' '-pthread'
 c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../libexec/gcc/x86_64-w64-mingw32/8.1.0/cc1plus.exe -quiet -v -I c:/Users/heiland/dev/roadrunner-win64-vs14-cp35m/include/rr/C -iprefix c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../lib/gcc/x86_64-w64-mingw32/8.1.0/ -D_MT -D_REENTRANT -U_REENTRANT main.cpp -march=skylake -mmmx -mno-3dnow -msse -msse2 -msse3 -mssse3 -mno-sse4a -mcx16 -msahf -mmovbe -maes -mno-sha -mpclmul -mpopcnt -mabm -mno-lwp -mfma -mno-fma4 -mno-xop -mbmi -msgx -mbmi2 -mno-pconfig -mno-wbnoinvd -mno-tbm -mavx -mavx2 -msse4.2 -msse4.1 -mlzcnt -mno-rtm -mno-hle -mrdrnd -mf16c -mfsgsbase -mrdseed -mprfchw -madx -mfxsr -mxsave -mxsaveopt -mno-avx512f -mno-avx512er -mno-avx512cd -mno-avx512pf -mno-prefetchwt1 -mclflushopt -mxsavec -mxsaves -mno-avx512dq -mno-avx512bw -mno-avx512vl -mno-avx512ifma -mno-avx512vbmi -mno-avx5124fmaps -mno-avx5124vnniw -mno-clwb -mno-mwaitx -mno-clzero -mno-pku -mno-rdpid -mno-gfni -mno-shstk -mno-avx512vbmi2 -mno-avx512vnni -mno-vaes -mno-vpclmulqdq -mno-avx512bitalg -mno-movdiri -mno-movdir64b --param l1-cache-size=32 --param l1-cache-line-size=64 --param l2-cache-size=8192 -mtune=skylake -quiet -dumpbase main.cpp -m64 -mthreads -auxbase main -std=c++11 -version -fomit-frame-pointer -fopenmp -o C:\Users\heiland\AppData\Local\Temp\cckCo2CT.s
GNU C++11 (x86_64-posix-seh-rev0, Built by MinGW-W64 project) version 8.1.0 (x86_64-w64-mingw32)
        compiled by GNU C version 8.1.0, GMP version 6.1.2, MPFR version 4.0.1, MPC version 1.1.0, isl version isl-0.18-GMP

GGC heuristics: --param ggc-min-expand=100 --param ggc-min-heapsize=131072
ignoring duplicate directory "c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/lib/gcc/../../lib/gcc/x86_64-w64-mingw32/8.1.0/include/c++"
ignoring duplicate directory "c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/lib/gcc/../../lib/gcc/x86_64-w64-mingw32/8.1.0/include/c++/x86_64-w64-mingw32"
ignoring duplicate directory "c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/lib/gcc/../../lib/gcc/x86_64-w64-mingw32/8.1.0/include/c++/backward"
ignoring duplicate directory "c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/lib/gcc/../../lib/gcc/x86_64-w64-mingw32/8.1.0/include"
ignoring nonexistent directory "C:/mingw810/x86_64-810-posix-seh-rt_v6-rev0/mingw64C:/msys64/mingw64/lib/gcc/x86_64-w64-mingw32/8.1.0/../../../../include"
ignoring duplicate directory "c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/lib/gcc/../../lib/gcc/x86_64-w64-mingw32/8.1.0/include-fixed"
ignoring duplicate directory "c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/lib/gcc/../../lib/gcc/x86_64-w64-mingw32/8.1.0/../../../../x86_64-w64-mingw32/include"
ignoring nonexistent directory "C:/mingw810/x86_64-810-posix-seh-rt_v6-rev0/mingw64/mingw/include"
#include "..." search starts here:
#include <...> search starts here:
 c:/Users/heiland/dev/roadrunner-win64-vs14-cp35m/include/rr/C
 c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../lib/gcc/x86_64-w64-mingw32/8.1.0/include/c++
 c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../lib/gcc/x86_64-w64-mingw32/8.1.0/include/c++/x86_64-w64-mingw32
 c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../lib/gcc/x86_64-w64-mingw32/8.1.0/include/c++/backward
 c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../lib/gcc/x86_64-w64-mingw32/8.1.0/include
 c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../lib/gcc/x86_64-w64-mingw32/8.1.0/include-fixed
 c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../lib/gcc/x86_64-w64-mingw32/8.1.0/../../../../x86_64-w64-mingw32/include
End of search list.
GNU C++11 (x86_64-posix-seh-rev0, Built by MinGW-W64 project) version 8.1.0 (x86_64-w64-mingw32)
        compiled by GNU C version 8.1.0, GMP version 6.1.2, MPFR version 4.0.1, MPC version 1.1.0, isl version isl-0.18-GMP

GGC heuristics: --param ggc-min-expand=100 --param ggc-min-heapsize=131072
Compiler executable checksum: 82f0c9785fd37a38ba7b7f8357369a82
COLLECT_GCC_OPTIONS='-march=native' '-fomit-frame-pointer' '-fopenmp' '-m64' '-std=c++11' '-I' 'c:/Users/heiland/dev/roadrunner-win64-vs14-cp35m/include/rr/C' '-v' '-o' 'myproj.exe' '-shared-libgcc' '-mthreads' '-pthread'
 c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../lib/gcc/x86_64-w64-mingw32/8.1.0/../../../../x86_64-w64-mingw32/bin/as.exe -v -I c:/Users/heiland/dev/roadrunner-win64-vs14-cp35m/include/rr/C --64 -o C:\Users\heiland\AppData\Local\Temp\ccD2LaVi.o C:\Users\heiland\AppData\Local\Temp\cckCo2CT.s
GNU assembler version 2.30 (x86_64-w64-mingw32) using BFD version (GNU Binutils) 2.30
COMPILER_PATH=c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../libexec/gcc/x86_64-w64-mingw32/8.1.0/;c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../libexec/gcc/;c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../lib/gcc/x86_64-w64-mingw32/8.1.0/../../../../x86_64-w64-mingw32/bin/
LIBRARY_PATH=c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../lib/gcc/x86_64-w64-mingw32/8.1.0/;c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../lib/gcc/;C:/Users/heiland/dev/roadrunner-win64-vs14-cp35m/lib/../lib/;c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../lib/gcc/x86_64-w64-mingw32/8.1.0/../../../../x86_64-w64-mingw32/lib/../lib/;c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../lib/gcc/x86_64-w64-mingw32/8.1.0/../../../../lib/;C:/Users/heiland/dev/roadrunner-win64-vs14-cp35m/lib/;c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../lib/gcc/x86_64-w64-mingw32/8.1.0/../../../../x86_64-w64-mingw32/lib/;c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../lib/gcc/x86_64-w64-mingw32/8.1.0/../../../
Reading specs from c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../lib/gcc/x86_64-w64-mingw32/8.1.0/libgomp.spec
COLLECT_GCC_OPTIONS='-march=native' '-fomit-frame-pointer' '-fopenmp' '-m64' '-std=c++11' '-I' 'c:/Users/heiland/dev/roadrunner-win64-vs14-cp35m/include/rr/C' '-v' '-o' 'myproj.exe' '-shared-libgcc' '-mthreads' '-pthread'
 c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../libexec/gcc/x86_64-w64-mingw32/8.1.0/collect2.exe -plugin c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../libexec/gcc/x86_64-w64-mingw32/8.1.0/liblto_plugin-0.dll -plugin-opt=c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../libexec/gcc/x86_64-w64-mingw32/8.1.0/lto-wrapper.exe -plugin-opt=-fresolution=C:\Users\heiland\AppData\Local\Temp\cc6VvYtJ.res -plugin-opt=-pass-through=-lmingwthrd -plugin-opt=-pass-through=-lmingw32 -plugin-opt=-pass-through=-lgcc_s -plugin-opt=-pass-through=-lgcc -plugin-opt=-pass-through=-lmoldname -plugin-opt=-pass-through=-lmingwex -plugin-opt=-pass-through=-lmsvcrt -plugin-opt=-pass-through=-lpthread -plugin-opt=-pass-through=-ladvapi32 -plugin-opt=-pass-through=-lshell32 -plugin-opt=-pass-through=-luser32 -plugin-opt=-pass-through=-lkernel32 -plugin-opt=-pass-through=-liconv -plugin-opt=-pass-through=-lmingwthrd -plugin-opt=-pass-through=-lmingw32 -plugin-opt=-pass-through=-lgcc_s -plugin-opt=-pass-through=-lgcc -plugin-opt=-pass-through=-lmoldname -plugin-opt=-pass-through=-lmingwex -plugin-opt=-pass-through=-lmsvcrt --sysroot=C:/mingw810/x86_64-810-posix-seh-rt_v6-rev0/mingw64 -m i386pep -Bdynamic -o myproj.exe c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../lib/gcc/x86_64-w64-mingw32/8.1.0/../../../../x86_64-w64-mingw32/lib/../lib/crt2.o c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../lib/gcc/x86_64-w64-mingw32/8.1.0/crtbegin.o -Lc:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../lib/gcc/x86_64-w64-mingw32/8.1.0 -Lc:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../lib/gcc -LC:/Users/heiland/dev/roadrunner-win64-vs14-cp35m/lib/../lib -Lc:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../lib/gcc/x86_64-w64-mingw32/8.1.0/../../../../x86_64-w64-mingw32/lib/../lib -Lc:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../lib/gcc/x86_64-w64-mingw32/8.1.0/../../../../lib -LC:/Users/heiland/dev/roadrunner-win64-vs14-cp35m/lib -Lc:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../lib/gcc/x86_64-w64-mingw32/8.1.0/../../../../x86_64-w64-mingw32/lib -Lc:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../lib/gcc/x86_64-w64-mingw32/8.1.0/../../.. -lroadrunner_c_api BioFVM_vector.o BioFVM_mesh.o BioFVM_microenvironment.o BioFVM_solvers.o BioFVM_matlab.o BioFVM_utilities.o BioFVM_basic_agent.o BioFVM_MultiCellDS.o BioFVM_agent_container.o pugixml.o PhysiCell_phenotype.o PhysiCell_cell_container.o PhysiCell_standard_models.o PhysiCell_cell.o PhysiCell_custom.o PhysiCell_utilities.o PhysiCell_SVG.o PhysiCell_pathology.o PhysiCell_MultiCellDS.o PhysiCell_various_outputs.o PhysiCell_pugixml.o PhysiCell_settings.o heterogeneity.o C:\Users\heiland\AppData\Local\Temp\ccD2LaVi.o -lstdc++ -lgomp -lmingwthrd -lmingw32 -lgcc_s -lgcc -lmoldname -lmingwex -lmsvcrt -lpthread -ladvapi32 -lshell32 -luser32 -lkernel32 -liconv -lmingwthrd -lmingw32 -lgcc_s -lgcc -lmoldname -lmingwex -lmsvcrt c:/Program Files/mingw-w64/x86_64-8.1.0-posix-seh-rt_v6-rev0/mingw64/bin/../lib/gcc/x86_64-w64-mingw32/8.1.0/crtend.o
heterogeneity.o:heterogeneity.cpp:(.text+0xc8e): undefined reference to `createRRInstance'
heterogeneity.o:heterogeneity.cpp:(.text+0xcaf): undefined reference to `__imp_loadSBML'
heterogeneity.o:heterogeneity.cpp:(.text+0xd0a): undefined reference to `__imp_getNumberOfReactions'
heterogeneity.o:heterogeneity.cpp:(.text+0xd25): undefined reference to `__imp_getNumberOfFloatingSpecies'
heterogeneity.o:heterogeneity.cpp:(.text+0xd40): undefined reference to `__imp_getNumberOfBoundarySpecies'
heterogeneity.o:heterogeneity.cpp:(.text+0xd5b): undefined reference to `__imp_getNumberOfGlobalParameters'
heterogeneity.o:heterogeneity.cpp:(.text+0xd76): undefined reference to `__imp_getNumberOfCompartments'
heterogeneity.o:heterogeneity.cpp:(.text+0xe83): undefined reference to `__imp_getFloatingSpeciesIds'
heterogeneity.o:heterogeneity.cpp:(.text+0xe8f): undefined reference to `__imp_stringArrayToString'
heterogeneity.o:heterogeneity.cpp:(.text+0xed0): undefined reference to `__imp_getFloatingSpeciesConcentrations'
heterogeneity.o:heterogeneity.cpp:(.text+0xff3): undefined reference to `createRRInstance'
heterogeneity.o:heterogeneity.cpp:(.text+0x1014): undefined reference to `__imp_loadSBML'
heterogeneity.o:heterogeneity.cpp:(.text+0x10b6): undefined reference to `createRRInstance'
heterogeneity.o:heterogeneity.cpp:(.text+0x10d7): undefined reference to `__imp_loadSBML'
heterogeneity.o:heterogeneity.cpp:(.text+0x1185): undefined reference to `createRRInstance'
heterogeneity.o:heterogeneity.cpp:(.text+0x11a0): undefined reference to `__imp_loadSBML'
heterogeneity.o:heterogeneity.cpp:(.text+0x160c): undefined reference to `__imp_getFloatingSpeciesConcentrations'
heterogeneity.o:heterogeneity.cpp:(.text+0x16f4): undefined reference to `__imp_setFloatingSpeciesConcentrations'
heterogeneity.o:heterogeneity.cpp:(.text+0x1725): undefined reference to `__imp_simulateEx'
collect2.exe: error: ld returned 1 exit status
make: *** [all] Error 1

heiland@LAPTOP-P2P2QC36 MINGW64 ~/git/pc4sbml/src (master)

```


