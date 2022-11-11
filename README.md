# Windowed Green function method for wave scattering by periodic arrays

This repository is solves the scattering problem by incident planewaves over a periodic domain in 2D and 3D for line arrays and surface arrays. 
The full-code is written in Julia and heavily relies on WaveProp packages (https://github.com/WaveProp) for geometry assembling, mesh construction and boundary integral equations discretizations.

WavePropBase and Nystrom dependencies should be first installed.

## Related work

The methodology has been presented and succesfully validated for the 2D problem with L-periodic line array in:

*Strauszer-Caussade T, Faria LM, Fernandez-Lado A, Pérez-Arancibia C. (2022) Windowed Green function method for wave scattering by periodic arrays of 2D obstacles. 
Stud Appl Math. ;1-39. https://doi.org/10.1111/sapm.12540*

## Warning!
Although the methodology converges for 2D it has not been fully tested and corrected for 3D. This corresponds to future work.
