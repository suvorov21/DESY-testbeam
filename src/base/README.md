# Analysis Base

# Analysis flow

1. Open the file
2. Guess the file format: 3D array or TRawEvent
3. Loop over entries in the input file
    1. If working with 3D array, subtract pedestal of 250
    2. Cast TRawEvent to TEvent and call reconstruction
    3. Call analysis Process() function to perform the analysis
4. Store output TObjects in the output file