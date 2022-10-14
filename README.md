# CiliaCoordination

Motile cilia may display remarkable traveling waves of rapidly succeeding ciliary beats across the tissue, much like Mexican waves in a soccer stadium. Those waves called [metachronal waves](https://en.wikipedia.org/wiki/Metachronal_rhythm).

Ciliacoordination is a collection of functions to quantify this collective beating of motile cilia. 

These codes support the analysis in the [preprint](https://doi.org/10.1101/2021.11.23.469646) on bioRxiv.

The raw data are technically challenging to acquire, so we will make the accompanying raw data available from Mendeley. 


### Summary
The **frequency analysis master** contains the main functions: 
- *fast fourier transform* to determine the ciliary beat frequency. 
- *pairwise coherence versus distance* to quantify synchronization. 
- *segmentation of frequency patches*
- *patchwise phase estimation* to visualize metachronal waves.
- *phase gradient calculations* to extract the wave directions and wave lengths. 

Note: although these codes can be run without aligning the data, it is highly recommended. 


**Additional codes** further accompany the preprint and are organized by figure. 
These codes are less well maintained than the frequency analysis master. 

### Further references
* there is complementary [toolbox](https://www.repository.cam.ac.uk/handle/1810/265273) for quantifying metachronal waves (Feriani, 2017)
