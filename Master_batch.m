%% Master-file loop

clearvars; close all; clc

% Select the file to analyze
[pathS] = uigetdir(pwd, 'Select the folder with files to analyze');
CBF.sourceP = fullfile(pathS, filesep);

% Select where to save the results
[pathP] = uigetdir(pwd,'Select where to save the results');
CBF.folderP = fullfile(pathP, filesep);

%% Settings 

CBF.w_min = 15; % Lower frequency cutoff 
CBF.caxis = [16 40]; % Upper and lower bound for any frequency plot. 
CBF.spatres = 0.15;  % [um/pixel] spatial resolution
CBF.signal_size = 400; % minimum number of pixels to be considered signal
CBF.SD = 3; % maximum standard deviation for a block of 9 pixels to be considered signal
CBF.metadata_ID = '_aligned'; % sometimes the json file does not match the data file exactly

% CBF.window = hamming(100); % window 100, noverlap 80, nfft 100 work very well! 
% CBF.noverlap = 80; 
% CBF.nfft = 100; 

% % Coherence prepping
% % Define the spatial resolution
% CBF.f = 5; % subsample factor
% % Bin into 100 bins.  
% CBF.n = 100;  % Number of bins to segment the powerspectrum into

%% Run the codes

analysis_loop(CBF)