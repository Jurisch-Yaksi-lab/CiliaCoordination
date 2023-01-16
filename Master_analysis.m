%% Master Analysis file
% All analysis performed on light transmission recordings of ciliary beating in the zebrafish nose
% It is possible to run the analysis on aligned data. 
% Assumes that the input file is a .mat file, and only contains one variable with the recording

%% Step 0: Select the recording

% Select the file to analyze
[fileS,pathS] = uigetfile('*.mat', 'Select the file to analyze');
[~,CBF.name,~] = fileparts(fileS); % Give a name

% Select where to save the results
[pathP] = uigetdir(pwd,'Select where to save the results');
CBF.targetP = fullfile(pathP, CBF.name, filesep);
[~, ~] = mkdir(CBF.targetP);

var = who(matfile(fullfile(pathS,fileS)));
data = double(load(fullfile(pathS,fileS),var{1}).(var{1}));

save([CBF.targetP, CBF.name, '_CBF_parameters'], 'CBF');
%% Step 1: Defining variables 
 
% Define variables and save them in one common structure
CBF.x = size(data,1); 
CBF.y = size(data,2); 
CBF.w_min = 15; % Lower frequency cutoff 
CBF.caxis = [16 40]; % Upper and lower bound for any frequency plot. 
CBF.spatres = 0.15;  % [um/pixel] spatial resolution
CBF.minsize = 400; % minimum number of pixels to be considered signal
CBF.SD = 3; % maximum standard deviation for a block of 9 pixels to be considered signal

% Input the frequency of acquisition
CBF.Fs = 108; % Frequency of acquisition

% Save the results
save([CBF.targetP, CBF.name, '_CBF_parameters'], 'CBF');

%% Step 2: Perform the fast Fourier transform 

% Fourier Transform (~20s to run the fft & ~200s with plotting)
[PowerSpec,~,~,~,~,CBF.picSD,~,CBF.mask,CBF.nframe, PeakPos] = bmf_sw_cr_fft_analysis(data,CBF);

% It is possible to run the frequency analysis with a predefined mask. 
% [PowerSpec,~,~,~,~,CBF.picSD,~,CBF.mask,CBF.nframe, PeakPos] = bmf_sw_cr_fft_analysis(data,CBF, mask);

% Save the results
save([CBF.targetP, CBF.name, '_CBF_parameters'], 'CBF');

%% Step 3: Coherence with example pixels 
 
% % Decide on reference pixels: uncomment to choose
% figure, imagesc(CBF.picSD); 
% ref = ginput; 
% ref = round(ref);

% % Alternatively, load a matfile with the reference pixels presaved
[fileR,pathR] = uigetfile('*.mat', 'Select the file with reference pixels to analyze');
var = who(matfile(fullfile(pathR,fileR)));
ref = double(load(fullfile(pathR,fileR),var{1}).(var{1}));

% Save the refence pixels
 CBF.ref = ref; 
 
% Set some variables specific to the coherence analysis
CBF.window = hamming(100); % window 100, noverlap 80, nfft 100 work very well! 
CBF.noverlap = 80; 
CBF.nfft = 100; 

% Compute coherence all with one reference pixel (~120s per reference pixel)
[Pxx, val] = cr_coherence_ref(data, CBF);

% Save results
save([CBF.targetP, CBF.name, '_CBF_parameters'], 'CBF');

%% Step 3.5: Coherence versus spectral density scatterplot 

% Plot coherence versus power spectral contribution (~30s)
map = cr_coherence_vs_psd(val, PeakPos, PowerSpec, CBF);


%% Step 4: Pairwise Coherence versus distance
% Define the spatial resolution
CBF.f = 5; % subsample factor

% Bin into 100 bins.  
CBF.n = 100;  % Number of bins to segment the powerspectrum into

% May not work outside of Kavlifarm (~ 10min in kavlifarm using the parallel pool). Computationally heavy. 
cr_coherence_all(data, Pxx, CBF)

save([CBF.targetP, CBF.name, '_CBF_parameters'], 'CBF');
% % 
% % % %% Step 4.1: Pairwise Coherence versus distance - One patch
% % % 
% % % patch_no = 6; 
% % % % May not work outside of Kavlifarm. Computationally heavy. 
% % % cr_coherence_all_patch(data, Pxx, CBF, patch_no)
% % %  

%% Step 5: Frequency segmentation 

% Bin the powerspectrum (<10s)
[freqsBinned] = cr_bin_power_spectrum(CBF, PowerSpec);

% Segment the nose into frequency patches.
I = freqsBinned.*CBF.mask;

% Run the frequency segmentation (~40s)
[CBF.lmatrix,complist] = identify_frequency_patches(I, CBF); 

% Save the results
save([CBF.targetP, CBF.name, '_CBF_parameters'], 'CBF');

%% Step 6: Calculate wavedirection and wavelength. 

% Decide whether you want to plot the patch check
check = false; 

% Set the wavelength histogram variables
CBF.max_wavelength = 200;
CBF.bin_size = 30;

% Find the phases for all patches (~50s)
[phase_patch, pos] = calc_phase_patch(data,complist, CBF);

% Calculate the gradient, wave direction and wavelength (~40s)
[CBF] = cr_patchwise_analysis(phase_patch, pos, PowerSpec, check, CBF); 

% Save the results
save([CBF.targetP, CBF.name, '_CBF_parameters'], 'CBF');