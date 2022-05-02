# Analysis Base

# Analysis flow

1. Open the file
2. Guess the file format: 3D array or TRawEvent
3. Loop over entries in the input file
    1. Convert input format to TEvent:
       1. Associate row and column for each hit 
       2. Subtract pedestals
    3. Call reconstruction to identify used and unused hits
    4. Call analysis Process() function to perform the analysis
4. Store output TObjects in the output file