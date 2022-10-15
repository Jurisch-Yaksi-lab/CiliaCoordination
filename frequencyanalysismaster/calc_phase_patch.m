function [phase_patch,pos] = calc_phase_patch(data, complist, CBF)

% This function calculates the phase angles for a given patch at the peak
% frequency of that patch
% data = the raw data (usually aligned)
% complist = patch border information from the frequency segementation code
% CBF = the parameter structure

% output
% phase_patch = the phases calculated for one frequency per patch
% pos = position of the patch's peak frequency in the powerspectrum.  (i.e. don't convert it.)

%% Define some variables
NumPatches = length(complist.PixelIdxList);
phase_patch =  nan(CBF.x,CBF.y,NumPatches);
pos = nan(1,NumPatches);
nyquist = floor(CBF.nframe/2); % round down the frame number
iw_min=1+round(CBF.w_min/CBF.Fs*CBF.nframe); % minimum frequency

%% Calculate phase angles per patch

for n = 1:NumPatches
    
    % Allocate for speed
    data_patch = zeros(length(complist.PixelIdxList{n}),CBF.nframe);
    
    % Collect raw data in the nth patch
    for i = 1:CBF.nframe
        d = squeeze(data(:,:,i)); d = d(:); % Rearrange them
        data_patch(:,i)=d(CBF.lmatrix == n);
    end
    
    % Calculate the Fourier Transform
    data_fft=squeeze(fft(single(data_patch(:,1:CBF.nframe)),[],2));
    
    % Find the peak position
    PowerSpec = abs(data_fft(:,1:nyquist+1));
    [~,maxpos]=max(mean(PowerSpec(:,iw_min:end)));
    pos(n) = maxpos+iw_min-1;
    
    % Extract the phase only
    phi = angle(data_fft(:,pos(n)));
    
    % Give them an index so that I can easily plot them back
    L = nan(CBF.x,CBF.y);
    L(complist.PixelIdxList{n}) = phi; % Insert list of phases here;
    
    % Save the phase_angles for a given patch
    phase_patch(:,:,n) = L;
end

 phase_patch = phase_patch+pi; % shift phases to interval [0, 2*pi]

%% Save data -------------------------------------------------------------
save(fullfile(CBF.targetP,[CBF.name,'_results_calc_phase_patch.mat']), 'phase_patch', 'pos','-v7.3');


