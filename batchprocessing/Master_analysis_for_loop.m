%% Master Analysis file
% All analysis performed on light transmission recordings of ciliary beating in the zebrafish nose
% It is possible to run the analysis on aligned data. 
% Assumes that the input file is a .mat file, and only contains one variable with the recording

%% Step 0: Select the recording

% Never run this step in the loop

%% Step 1: Defining variables 
 
% Never run this step in the loop

%% Step 2: Perform the fast Fourier transform 

% Fourier Transform (~20s to run the fft & ~200s with plotting)
[PowerSpec,~,~,~,~,CBF.picSD,~,CBF.mask,CBF.nframe, PSD, PeakPos] = bmf_sw_cr_fft_analysis(data,CBF);

% It is possible to run the frequency analysis with a predefined mask. 
% [PowerSpec,~,~,~,~,CBF.picSD,~,CBF.mask,CBF.nframe, PSD, PeakPos] = bmf_sw_cr_fft_analysis4(data,CBF, mask);

% Save the results
save([CBF.targetP, CBF.name, '_CBF_parameters'], 'CBF');

%% Step 3: Coherence with example pixels 
 
% Never run this step in the loop

% % % % Decide on reference pixels: uncomment to choose
% % figure, imagesc(CBF.picSD); 
% % ref = ginput; 
% % ref = round(ref);
%  CBF.ref = ref; % Save the refence pixels
% % 
% % Set some variables specific to the coherence analysis

% 
% % Compute coherence all with one reference pixel (~120s per reference pixel)
% [Pxx, val] = cr_coherence_ref(data, CBF);
% 
% % Save results
% save([CBF.targetP, CBF.name, '_CBF_parameters'], 'CBF');

%% Step 3.1: Coherence versus spectral density scatterplot 

% Never run this step in the loop

%% Step 4: Pairwise Coherence versus distance
% Define the spatial resolution
CBF.f = 5; % subsample factor

% Bin into 100 bins.  
CBF.n = 100;  % Number of bins to segment the powerspectrum into
% 
% % May not work outside of Kavlifarm (~ 10min in kavlifarm using the parallel pool). Computationally heavy. 
% cr_coherence_all(data, Pxx, CBF)
% 
% save([CBF.targetP, CBF.name, '_CBF_parameters'], 'CBF');
% % 
%% Step 4.1: Pairwise Coherence versus distance - One patch

% Never run this step in the loop

%% Step 5: Frequency segmentation 

% Bin the powerspectrum (<10s)
[freqsBinned] = cr_bin_power_spectrum(CBF, PowerSpec);

% Segment the nose into frequency patches.
CBF.minsize = 400; 
I = freqsBinned.*CBF.mask;

% Run the frequency segmentation (~40s)
[CBF.lmatrix,complist] = identify_frequency_patches(I, CBF); 

% Save the results
save([CBF.targetP, CBF.name, '_CBF_parameters'], 'CBF');

%% Step 6: Calculate wavedirection and wavelength. 

% Decide whether you want to plot the patch check
check = false; % Do not run this during the loop

% Set the wavelength histogram variables
CBF.max_wavelength = 200;
CBF.bin_size = 30;

% Find the phases for all patches (~50s)
[phase_patch, pos] = calc_phase_patch(data,complist, CBF);

% Calculate the gradient, wave direction and wavelength (~40s)
[CBF] = cr_patchwise_analysis(phase_patch, pos, PowerSpec, check, CBF); 

% Save the results
save([CBF.targetP, CBF.name, '_CBF_parameters'], 'CBF');
