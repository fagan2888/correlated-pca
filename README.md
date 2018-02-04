Folder contains the codes to reproduce the results of our paper illustrating the Guarantees for Singular Value Decomposition when the noise is anisotropic and possibly even data-dependent.


If you use the results please cite the following paper

[1] Namrata Vaswani and Praneeth Narayanamurthy, "Finte Sample Guarantees for PCA in non-isotropic and data-dependent noise", arxiv:1709.06255, 2017. 

Main contents of the folder

1. DemoPhaseTransition.m -- Main demo file to generate the phase transition plots for seeing the dependence of the number of samples required to achieve a desired error level versus either the rank or the signal dimension

2. DemoBoundValidation.m -- Main demo file to verify the validity, and tightness, of the bound predicted theoretically versus the actual average and maximum subspace errors obtained.

3. SimpleEVD.m -- main function to perform the SVD algorithm 
4. Calc_SubspaceError -- function to calculate the subspace errors between two subspaces
5. GenerateDatFile.m -- file to generate data files in accordance with TikZ format.

Please let me know if there are any clarifications/suggestions/mistakes at pkurpadn and iastate edu (make obvious changes)


