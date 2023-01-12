We provide here additional information about the parameters for phase mapping for the identification of metachronal wave direction

# 1. How are frequency maps segmented in space into equal-frequency patches? 

We realized early on that the power spectrum has a higher frequency resolution than necessary (0.1Hz) and that small variations in frequency could represent noise rather than cilia beating at different frequencies. Hence, we recommend to bin the power spectrum (using a bin width of 0.54Hz). 
The precise choice of bin widths had only a small impact on the CBF segmentation (note that the spatial distribution of CBF remain similar). For bin widths grossly outside of this range, too few or too many CBF values result (see Figure below for a visual example). 
Based on our experiments with the multiciliated epithelium in the nose and brain of the zebrafish we would recommend to use a frequency binning of circa 0.5Hz


![image](https://user-images.githubusercontent.com/17269484/212029975-369f8fe5-eee4-44bb-a5b6-be38fd35953b.png)
_Legend: (A) Impact of the binning of the power spectrum on CBF values shown on CBF heatmaps (top) and histograms (bottom). We used a binning of 0.54Hz for all our analysis due to their minimal impact and good coverage of CBF values_

# 2. How are time signals treated? 
Time signals are not altered by our analysis. Binning is performed in the frequency domain only. 
The only threshold applied in our analysis was for the CBF analysis. CBF were identified between 15Hz and half of the frequency of acquisition (Fs = 100Hz; Nyquist = 50Hz), based on the Nyquist formula. We used 15Hz as a threshold for the CBF based on the analysis of >200 samples (see Reiten et al., 2017 for more information, which identified 15Hz as the lowest frequency in the nasal multiciliated epithelium. 
 
This thershold value needs to be adapted by future users to the tissue of interest based on their CBF. For instance, a cut-off of 10Hz was used for the ependymal cells of the zebrafish brain based on our earlier work (D'Gama et al., 2021).

# 3. How do the subtleties compare in a different system? 
We have used our analysis pipeline on different multiciliated tissue, i.e., the zebrafish nasal ciliated epithelium and the ependymal cells of the adult zebrafish brain, and different microscopes. 
Since ciliated cells in various tissues have different properties with regard to their CBF, apical size and numbers of cilia, we recommend experimenters to adapt these values for each ciliated system. Similarly, these parameters should be adapted to each acquisition system, and reflect the spatial resolution of the recordings. If you need additional help for implementing this pipeline to your data, please do not hesitate to contact us.

# 4. What is the optimal choice of patch-size?
For quantification of wave direction and wavelength, we recommend to opt for frequency patches that are as large as possible but display sufficiently homogeneous CBF. 
For our experiments, we choose the patch-size by 
1)	Binning the frequencies (as described above)
2)	Setting a minimum patch size of 400 pixels, which corresponds to circa 9µm2, roughly corresponding to the area swept over by one cilium. We have also used higher minimum patch size, eg 800pixel but as shown in teh Figure below, this was too stringent.

![image](https://user-images.githubusercontent.com/17269484/212031613-3010e383-fa86-4e3d-af89-f7b8d859d13b.png)
_Legend (B) Segmentations into frequency patches with or without binning and with a minimum patch size of 400 or 800 pixels. Note that a binning of 0.54Hz and minimum size of 400 pixels (9 µm2) provides the best segmentation of the ciliated epithelium with a reasonable number of patches._

# 5. What are the criteria for extracting phases? - What is the tolerance in frequency variation across regions?
Phases are extracted from the Fourier transform of the intensity time series of each pixel evaluated at the binned frequency of the segmented patch of that pixel. Using the same evaluation frequency for each pixel in a frequency patch is important to ensure that extracted phases are comparable across the patch. 
We tested this algorithm extensively on synthetic data of simulated coupled/uncoupled noisy phase oscillators and confirmed that we can robustly estimate their phase from time-series of finite length. 
We note a trade-off between precision and accuracy: using longer time series will increase the precision of phase estimation at the expense of reduced temporal resolution and accuracy. 

We additionally tested an alternative algorithm to extract cilia phase based on the Hilbert transform of intensity time series which gave results that were consistent with the Fourier-based algorithm. We have not pursued teh analysis with teh Hibler transform as results from this alternative algorithm were slightly noisier and less easy to interpret. 

# 6. Are 30s recordings ideal for this method?
For consistent results, the same duration of time series should be used for extracting phases and calculating coherence scores. 
The choice of optimal time duration of time series for the analyses represents a trade-off between increasing precision (longer duration of time series) versus accuracy (shorter time series). As a rule of thumb, time series should comprise a sufficient number of oscillation cycles (~500 in our case). For very long time series, additional post-processing steps can become necessary to reduce drift. 
As shown below in teh figure, using a time duration of 10-20s instead of 30s gave consistent results for the coherence score but background values were higher. Conversely, using a longer time duration of up to 240s did not improve the results for the coherence score. As longer recording time require a perfectly stable recording with zero drift which is very difficult to achieve, we recommend the users to select a time window of 30 sec for the coherence score to minimize background and impact of potential drift.

 ![image](https://user-images.githubusercontent.com/17269484/212032229-257c6d83-f4d8-44f6-8946-ba6373f08640.png)
 
_Legend: Coherence analysis for 3 reference pixels for different recording lengths (10-240s). Note that increasing recording length reduces the background values, but do not changes the overall coherence patterns. We recommend a duration of 30s to increase signal-to-noise ratio of the coherence score._


