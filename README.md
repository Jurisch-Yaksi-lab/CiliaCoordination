# CiliaCoordination

Motile cilia may display remarkable traveling waves of rapidly succeeding ciliary beats across the tissue, much like Mexican waves in a soccer stadium. Those waves are called [metachronal waves](https://en.wikipedia.org/wiki/Metachronal_rhythm).

Ciliacoordination is a collection of functions to quantify this collective beating of motile cilia. 

These codes support the analysis in the [preprint](https://doi.org/10.1101/2021.11.23.469646) on bioRxiv.

The raw data are available in [Mendeley](https://data.mendeley.com/datasets/th35c5833g/1), doi: 10.17632/th35c5833g.1 


### Summary
The **frequency analysis master** contains the main functions: 
- *fast fourier transform* to determine the ciliary beat frequency. 
- *pairwise coherence versus distance* to quantify synchronization. 
- *segmentation of frequency patches*
- *patchwise phase estimation* to visualize metachronal waves.
- *phase gradient calculations* to extract the wave directions and wave lengths. 

Note: although these codes can be run without aligning the data, it is highly recommended. 

**Additional codes** further accompany the preprint and are organized by figure. 
These codes are less polished than the frequency analysis master. 

All code was run on MATLAB version R2020B with toolboxes: image_toolbox, signal_toolbox, and statistics_toolbox.
Most codes can be run on a desktop computer, but the *pairwise coherence versus distance* code required the university's server farm.

### Running CiliaCoordination

The code can be run 1) per recording 2) in batch. 

The single recording code *Master_analysis* is a script organized in sections/steps.

  Step 0: In the first section, one selects a ciliary beat recording with a .mat extension. 
  Step 1: Define variables for analysing the recording. Mind to adjust the spatial resolution to the recorded ciliary beating. Unfortunately, the best values for lower frequency cutoff, minsize, and standard deviation depend on the tissue analyzed. 
  Step 2: Run the fast Fourier transform. There is an option to supply a mask. For instance, when there is an obvious piece of tissue debris present in the recording.
  Step 3: Run this to calculate coherence for a single pixel with all others in the recording. Note: it is necessary to supply the coordinates of the example pixel. The given settings work well usually. 
  Step 3.5: Compare the coherence of the last step with the powerspectral density in a scatterplot. 
  Step 4: Calculate the pairwise coherence of all signal pixels within the recording. Note: computationally heavy. It is better to run this on a cluster or server farm. Requires settings (window, noverlap, and nfft) given in step 3. 
  Step 5: Identify patches of ciliary beating at a similar frequency. Note: This step and onwards does not depend on step 3-4. 
  Step 6: Extract the phase for each frequency patch separately and extract the phase parameters wave direction / wavelength. It is possible to plot a panel of sanity checks for each patch by running check = true; 
  
The batch code *Master_batch* is a function running a part of the code for an entire folder. 

  Make sure that the setting are correct before running the batch code! 
  Running the code in batch will create a folder per recording. 
  The batch code will skip steps 0-1 and 3-4.


### Further references
* codes supporting the model are available in https://github.com/icemtel/reconstruct3d_opt, https://github.com/icemtel/stokes, and https://github.com/icemtel/carpet.
* there is complementary [toolbox](https://www.repository.cam.ac.uk/handle/1810/265273) for quantifying metachronal waves using different algorithms (Feriani, 2017)
