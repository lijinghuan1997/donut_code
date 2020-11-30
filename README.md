# donut_code

First, run the donut_equilibrium_model.m to get the data_final.mat file, which contains the electromagnetic profiles. 
Details about this model are shown in https://github.com/lijinghuan1997/cavitymodel.

Then, transform the data_final.mat into the text file. (This part is also contained in the matlab file.)
This text file is the electromagnetic field distritution, with the resolution 1 meter (the first column represents the distance from the cavity center).

Finally, we run the single-particle fortran code, which needs to load the information in the text file. 
The fortran code can output the trajectories and the corresponding PSDs of electrons, as well as the perturbed  into different text files (details in the outf section).

If you want to change the initial state of the single electron, you can change the energy/alpha(pitch angle)/beta(phase angle) around Line 265.
If you want to change the features of the evolving cavity, you can change the eq(R_0),eq2(Z_0),mm(corresponding to deltaE) around Line 266.
