# monte-carlo-rte-fluorescent-medium
Solution of radiative transfer equation in fluorescent medium including  luminescent cascade and non-fluorescent absorption

This code is used in the articles please cite them if you use the code in your study;

https://doi.org/10.1088/2053-1591/ab28b8 Improving photosynthetic efficiency using greenhouse coatings with scattering and fluorescent pigments

https://doi.org/10.1016/j.biosystemseng.2020.02.007 Improving crop production in solar illuminated vertical farms using fluorescence coatings

It solves RTE in fluorescent medium.

It can handle scattering in the medium due to phosphor particles.
It can handle luminescent cascade where absorption and emission spectra of fluorescence overlaps
It can handle non-flourescent absorbing medium by modifiying quantum yield accordingly and calculating refraction angle at the air-medium interface considering both n and k of the medium.
It can handle up-conversion and down-conversion

The code works on both MATLAB and Octave. 

main.m code is the main code to execute. Concentration, size and thickness of the phosphor can be altered inside this file.
