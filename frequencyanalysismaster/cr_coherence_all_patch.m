function [] = cr_coherence_all_patch(data, Pxx, CBF, patch_no)
% A function to calculate the coherence of all pairwise pixel coherences in that video.  

%% Define variables 
nfreqs = CBF.nfft/2+1;
iw_min=1+round(CBF.w_min/CBF.Fs*CBF.n);


%% Calculate both the pairwise cross-spectrum and pairwise mscoherence
% Subsample 
% masked = CBF.mask(1:CBF.f:end,1:CBF.f:end); masked = ~isnan(masked);
m = (CBF.lmatrix == patch_no); % Bring this outside?
masked = double(m); masked(masked == 0) = NaN; 

% Define matrix size
nloops = length(masked(masked==1))-1;

% Extract the raw data for one patch only
[data_mask] = extract_raw(); 

% Extract the powerspectra for one patch only
[Pxx_mask] = extract_ps(); 

%% Calculate both pairwise mscoherence and pairwise distance

% Calculate pairwise coherence 
Cxy = pairwise_mscohere();

% Calculate the pairwise distance 
Dist_M = pairwise_distance();

%% Calculate the coherence matrix
[val,~] = max(Cxy(:,:,iw_min:end),[],3); % 
Vcoh = squeeze(val); 
Vcoh(isnan(Vcoh)) = 0; 
Vcoh = Vcoh + Vcoh'; % Copy values to the bottom half of the matrix
Vcoh(1:nloops+2:end) = ones(1,nloops+1); % Set all diagonal values to one. 

figure, imagesc(Vcoh), colorbar
saveas(gcf, [CBF.targetP,CBF.name, '_figure_coherence_matrix.png']);

%% Calculate the pairwise distance matrix

Dist = squeeze(Dist_M);
Dist(isnan(Dist)) = 0; 
Dist = Dist + Dist'; % Copy values to the bottom half of the matrix
Dist(1:nloops+2:end) = zeros(1,nloops+1); % Set all diagonal values to zero.  

figure, imagesc(Dist), colorbar
saveas(gcf, [CBF.targetP,CBF.name, '_figure_distance_matrix', sprintf('_%d', patch_no),'.png']);

%% Coherence versus distance 

C = Vcoh(~(triu(Vcoh,1) == 0)); 
D = Dist(~(triu(Dist,1) == 0)); 
% D = D*CBF.f*CBF.spatres; % Conversion to micrometer
D = D*CBF.spatres; % Conversion to micrometer


figure, scatter(D, C,30, 'filled', 'k', 'MarkerFaceAlpha', 0.02); 
xlabel('distance in µm'); ylabel('coherence'); 
title('coherence versus distance');

save([CBF.targetP, CBF.name, '_result_coherence_versus_distance', sprintf('_%d', patch_no),'.mat'], 'Vcoh', 'Dist');
saveas(gcf, [CBF.targetP,CBF.name, '_figure_coherence_vs_distance', sprintf('_%d', patch_no),'.png']);
% print('G:\Christa\CoherenceVSdistance', '-painters', '-depsc');

%% Nested functions

    function [data_mask] = extract_raw()
%         data_small = data(1:CBF.f:end,1:CBF.f:end,:);
        data_mask = zeros(nloops+1,CBF.nframe);
        
        for iframe = 1:CBF.nframe
            d = data(:,:,iframe);
            data_mask(:,iframe) = d(m);
        end
    end

    function [Pxx_mask] = extract_ps()
%         Pxx_small  = Pxx(1:CBF.f:end, 1:CBF.f:end,:);
        Pxx_mask = zeros(nloops+1,nfreqs);
        
        for ifreqs = 1:nfreqs
            p = Pxx(:,:,ifreqs);
            Pxx_mask(:,ifreqs) = p(m);
        end
    end

    function [Cxy] = pairwise_mscohere()
        
        % Calculate coherence, but do not calculate it twice per pixel
        Pxy_mask = nan(nloops+1,nloops+1,nfreqs);
        
        tic
        parfor jloop = 1:nloops
            
            % Every pixel is ref_pix once
            ref_pix = data_mask(jloop,:);
            
            test_pix = data_mask(jloop+1:end,:); %#ok<PFBNS> % calculate only between a pixel and all other pixels once.
            [Pxy1,~] = cpsd(ref_pix, test_pix',CBF.window,CBF.noverlap,CBF.nfft,CBF.Fs);
            Pxy_mask(jloop,:,:) = [nan(jloop,nfreqs); Pxy1'];
            
        end
        toc
        
        
        % Calculate the autospectra
        Pxx_M = nan(nloops+1,nloops+1,nfreqs);
        tic
        for jloop = 1:nloops
            
            ref_pix = Pxx_mask(jloop,:);
            test_pix = Pxx_mask(jloop+1:end,:);
            
            Pxx_M(jloop,:,:) = [nan(jloop,nfreqs); ref_pix.*test_pix];
            
        end
        toc
        
        % Coherence
        Cxy = abs(Pxy_mask).^2 ./ (Pxx_M);  % Should be from 0 to 1!
    end

    function [Dist_M] = pairwise_distance()

x = size(masked,1); 
y = size(masked,2); 
Dist_M = nan(nloops+1,nloops+1);

% Generate an index matrix
for i = 1:x
    for j = 1:y      
Idx{i,j} = [j,i];
    end    
end 

% Extract the distances for only those pixels that I subsampled
Idx_mask = Idx(m);

tic
for j = 1:nloops
    dist_k = [];
    ref_pix = Idx_mask{j,:};
    test_pix = Idx_mask(j+1:end,:);
    
    for k = 1:length(test_pix)
        
        test = test_pix{k};
        dist_k(k) = sqrt((ref_pix(1)-test(1))^2 + (ref_pix(2)-test(2))^2);
    end
    Dist_M(j,:) = [nan(j,1); dist_k'];
    
end
toc


    end 
end 