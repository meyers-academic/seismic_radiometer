# Newtonian noise code

This code computes newtonian noise given seismic maps from the radiometer. It is based on matlab code currently living on the CIT clusters, and UMN git.

## Files

Core files:

- `CoherentSeismicField.py` Class to read and store information about the seismic field -- radiometer map, eigenfunction information, speeds, etc.
- `TestMass.py` Class to compute properties of a test mass. This is where the work to compute acceleration from the seismic field is done.

Scripts:
- `calc_NN.py` Script to compute Newtonian noise, using the above two classes. **NOTE** Contains references to files on my local computer, but I think they are standard files a user should have...

Other files
- `NN_plotter` Code to make plots of seismic maps
- `notes.md` Some notes to self I used to organize writing the code. 