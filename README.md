# physicell_libroadrunner

Demonstrate the use of SBML metabolic models in PhysiCell (using libRoadrunner).

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
