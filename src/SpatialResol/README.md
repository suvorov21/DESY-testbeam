# Analysis Flow

## initialisation

## Process Event
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

