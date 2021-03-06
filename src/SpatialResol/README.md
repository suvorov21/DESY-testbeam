# Spatial resolution analysis

Spatial resolution analysis is done with iterative process. 
The detailed description can be found at [arXiv:2106.12634](http://arxiv.org/abs/2106.12634).

# Usage

```
SpatialResol.exe
-i, --input  <f>    input file name or a list of files
-o, --output  <f>    output file name
--nthread   <int>   number of threads to use [1-5]

-t, --iter <int>    specify the iteration number.
--param, p <file>   files with parameters: selection., etc
--prev     <file>   specify the file from the previous iteration if different from default
--corr              whether to perform a fit w/ and w/o particular column and take geometrical mean

--start     <i>     start from event i
--end       <i>     end with event i

-b                  batch mode. Skip event plotting
-v          <i>     set verbosity level
-r                  rewrite the output file
-d                  debug mode. Run over 30 events.
```
# Analysis Flow

## Initialisation
- Define TTree vars

- Define histograms

## Process Event

```
#### step 1 ####################################################################
# clusterize and ommit 1st and last columns/rows
# E.g. form columns for horizontal clusters ot diagonals for tracks with angle close to 45 degree.
ClusteriseTrack()

get robust clusters (truncate, etc.)

#### step 2 ####################################################################

for cluster:clusters
    treat cross-talk
    get centre of charge position
    call fitter for cluster

#### step 3 ####################################################################

call track fitter for clusters_cleanm

extract track curvature

#### step 4 ####################################################################

if (correction)
    for cluster:clusters_clean
        fit track w/o each cluster

#### step 5 ####################################################################

for cluster:clusters_clean
    fill residuls

#### step 6 ####################################################################

for cluster:clusters
    for pad:cluster
        fill PRF

```

## Verbosity

```
v_progress = 1,
v_event_number,
v_analysis_steps,
v_fit_details,
v_residuals,
v_prf
```
