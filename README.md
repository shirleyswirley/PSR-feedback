# PSR (particle size remineralization) feedback

The code in this repository reproduces figures and tables in the following paper:
<br>Leung, S., Weber, T., Cram, J. A., & Deutsch, C. Variable phytoplankton size distributions reduce the sensitivity of global export flux to climate change. <i>Submitted to Biogeosciences.</i>

Please cite the above paper and the code itself (see here: http://doi.org/10.5281/zenodo.3783471) if you use any of it.

This code was written using MATLAB 9.6.0.1174912 (R2019a) Update 5.

How to run:
1. Download this repository.
2. Download model output and data from http://doi.org/10.5281/zenodo.3783473 and place it into this repository's folder.
3. Download colorbrewer package from https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab and place it into the utils folder.
4. Download m_map package from https://www.eoas.ubc.ca/~rich/map.html and place it into the utils folder.
5. To recreate the figures and tables from the above paper, simply start up MATLAB, make sure your current working directory is the repository's folder, and run the scripts in figure_code (order doesn't matter). Each of these scripts creates and saves/prints out a different figure in the paper.
