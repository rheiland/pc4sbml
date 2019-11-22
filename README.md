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

Then I attempted various incantations in `Make-mingw` to build. So far, I cannot get the linker to find the libRR functions. I've tried to run `make -f Make-mingw` in both a `Command Prompt` shell and a `MINGW64/Git Bash` shell.
```
C:\Users\heiland\git\pc4sbml\src>make -f Make-mingw
g++ -march=native  -fomit-frame-pointer -fopenmp -m64 -std=c++11 -I/c/Users/heiland/dev/roadrunner-win64-vs14-cp35m/include/rr/C  -L/c/Users/heiland/dev/roadrunner-win64-vs14-cp35m/lib  -lroadrunner_c_api -o myproj BioFVM_vector.o BioFVM_mesh.o BioFVM_microenvironment.o BioFVM_solvers.o BioFVM_matlab.o BioFVM_utilities.o BioFVM_basic_agent.o BioFVM_MultiCellDS.o BioFVM_agent_container.o  pugixml.o PhysiCell_phenotype.o PhysiCell_cell_container.o PhysiCell_standard_models.o PhysiCell_cell.o PhysiCell_custom.o PhysiCell_utilities.o PhysiCell_SVG.o PhysiCell_pathology.o PhysiCell_MultiCellDS.o PhysiCell_various_outputs.o PhysiCell_pugixml.o PhysiCell_settings.o  heterogeneity.o main.cpp
heterogeneity.o:heterogeneity.cpp:(.text+0xc8e): undefined reference to `createRRInstance'
heterogeneity.o:heterogeneity.cpp:(.text+0xcaf): undefined reference to `__imp_loadSBML'
heterogeneity.o:heterogeneity.cpp:(.text+0xd0a): undefined reference to `__imp_getNumberOfReactions'
heterogeneity.o:heterogeneity.cpp:(.text+0xd25): undefined reference to `__imp_getNumberOfFloatingSpecies'
...
```


I also tried, in a MinGW/Git Bash shell:
```
$ export LIBRARY_PATH
$ LIBRARY_PATH="C:/Users/heiland/dev/roadrunner-win64-vs14-cp35m/lib/"
```

