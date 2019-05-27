# Python Newtonian noise

## Workflow

- Read in the coherent map
- Read in other data: eigenfunctions,density, velocities
- Renormalize based on power spectrum (?), renormalize r wave to only use equator pixels
- Loop over polarizations, and sky directions. For each polarization/pixel, compute complex-valued acceleration response of test mass (need to sum different contributions). In GGN, the flags we need are: 7,8,9,11,12,13,14,15,16,17. (NOTE: This is different from matlab: we dont need to rescale by sqrt(2) to convert rms to amplitude, and we need to include a phase. DOUBLE CHECK that the flags I introduced (a) are correct and (b) include the phase)
- Sum accelerations coherently over pixel and polarization (NOTE: This is different from matlab: we  summed coherently over directions, and incoherently over polarizations)
- Keep track of accelerations from each polarization so we can make a budget later
- Compute newtonian noise using the accelerations
- Plot results