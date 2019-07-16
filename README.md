# DESY beam test analysis package

The package for the DESY beam test data analysis.

## Working with macros
You could run the analysis with a simple macros (like at CERN analysis). You could fine the macros template under macros/TutorialMacro.cxx

### Compilation
The compilation is possible from the macros folder. You need to specify the macro you need as a source. You are free to create subdirectories in the macro folder with your analysis. Just call the compilation from the macro folder and specify the proper path.
Make sure that the proper compiler and ROOT version are used. ROOT >= v6.16 was tested.
To setup the environment at LXPLUS you can use setup.sh without changes.
```bash
cd macros
make SRC="YourMacro.cxx"
```

The output will be stored as bin/YourMacro.exe

### Running
The macros can read both single root file and list of root files. To run the macros call
```bash
./bin/YourMcro.exe -i input_file -o output_path
```
The ROOT help could be called with -h flag. The macro help will be called with flag '-m'.
The common view is the following
```
../bin/YourMacro.exe usage

   -i <input_file>      : input file name with a path
   -o <output_path>     : output files path
   -b                   : run in batch mode
   -v <verbose_lvel>    : verbosity level
   -d                   : test mode. run over first 30 events
   -h                   : print ROOT help
   -m                   : print ../bin/YourMacro.exe help
```
## HighLevel tool
For the DESY beam test analysis the HighLevel tool was created. The macros that we use before were good for quick checks. It's easy to create them from scratch and implement some quick analyses. During the long period of analysis we created a lot of macros and made them looking awful. The new flexible tool should allow us to create nice analysis algorithms.

The main idea is to separate the routine procedures (e.g. opening files, looping, writing the output etc.) from the analysis itself. So now the analyzer should only define the histos/canvases he want to store and the logic how to fill it.

### Tool structure
Rough scheme of the package
![](doc/html/_spatial_resol_ana_8cxx__incl.png)

[WIP]
1. Create your analysis inheriting from AnalysisBase
2. Add it in the Makefile for compilation. Cross-check the compilation.
3. Define histos you are interested in (inside Initialise()). Add them into _output_vector. The will be written to output file automaticaly.
4. Define the ProcessEvent() function. Fill your histogram for each event.
5. Enjoy the output!

## Geometry
The geometry information is stored in utils/Geom.hxx. The coordinate system and the pad size are the following:
```
Micromegas:

   31 ______________
     |              |
     |              |
 y   |              |
     |              |         Pad:
     |              |          __
   0 |______________|         |__| 1.1 cm
     0               35        1.0 cm
             x
```
