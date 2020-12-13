# Plotters

# Event display
Event display can be run with
```
root
.x EventDisplay.C("~/DATA/R2019_06_14-16_42_26-000.root")
```
The time binning of the data (510/511 bins) will be chosen automatically.

On the top left of the right frame you will see the the window with event number. It will be updated with next/previous button. Also you can out any number with keyboard and press Draw button to see the event you want.

P.S. if you know hot to move this window with number to bottom line with buttons, please, contribute:)

# Figures plotting with Jupyter

ROOT should be installed with the proper python version. Commands `python3` and `import ROOT` should work from terminal. That should work by default at LXPLUS

T2K style is optional and depends on the [package](https://gitlab.com/suvorov21/t2kstyle_py). If styling is not required -- simply comment the line.

Jupyter can be run either on LXPLUS or at your local machine with mounted path towards the analysis results. Change the `path` variable at the beginning of each notebook to set a path properly.

## Other tips

If Java script is not working, specify `%jsroot off` in the cell.

To make it working specify proper address in the ~/.jupyter/jupyter_notebook_config.py, e.g.
`c.NotebookApp.extra_static_paths=['/Users/suvorov/SOFT/CERN/ROOT/root-6.18.04-buildi_py3/js/']`
