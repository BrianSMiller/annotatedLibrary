This zip file contains Matlab functions for subsampling from a long-term recording.

To use these functions, download and extract the compressed zip file onto your computer somewhere and make sure to add all of the folders within the zip file to your Matlab path. 

Before you can use the function you will need to create a parameter file. The easiest way to do this will be to adapt the existing parameter file,  parameters_casey2014.m. You'll want to change the code, input folder, and output folder to match the location of your acoustic data.

The script makes use of a small package under development called soundFolder, which I wrote for loading bits of acoustic data from a long-term recording. Right now soundFolder only works with folders full of wav files where the timestamp for the start time of the file is stored in the filename (e.g. files generated from PAMGuard, Ishmael, Raven, AAD moored recorders). However, with a little bit of effort it should be possible to adapt soundFolder to work on (H)ARP xWav files, and probably other types of data.

To run the function, simply load the parameter file, then run the subsampling function passing the parameters as the first argument.

params_casey2014; 							% Load params into the workspace
subsampleLowFrequencyFixedRate(params);		% Run the subsampling routine using above parameters

