# Analysis Flow

## initialisation

## Process Event

```
#### step 1 ###############################################################################
# clusterize and ommit 1st and last columns/rows
if diagonal
    Diagonolize --> vector<TCluster> clusters
else
    Colonize --> vector<TCluster> clusters

# get robust clusters (truncate, etc.)

#### step 2 ####################################################################

for cluster:clusters
    get CoC position
    call fitter for cluster

if diagonal
    ommit 1st and last diags
    average 2 clusters in a row --> cluster_clean
else
    clusters --> clusters_clean

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
### Obsolete
```
for col:track
    init cluster params
    for pad:col
        make cluster
        fit cluster
    fit track

extract track curvature

if (correction)
    for col:track
        fit track w/o each column

for col:track
    fill residuals
    fill X scan

    for pad:col
        fill PRF
```