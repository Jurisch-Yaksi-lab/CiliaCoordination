function [Pxx, val] = cr_coherence_ref(data, CBF)
% A function to calculate the coherence of one reference pixel with all other pixels
% in that video.  

% The code is largely based on the mscohere() function in Matlab, but
% slightly optimized for speed. 
% https://se.mathworks.com/help/signal/ref/mscohere.html

% data = the raw data and will be converted to double. 
% Fs = the frequency of acquisition 
% ref = the locations of the reference pixels. 
%       It is a 2 collumn array that contains two integers. 
%       The first position is its collumn and the second position is its
%       row, which is the order of output when using a function like
%       ginput. 
% targetP = a path to the folder where the plots and data should be saved. 

%% Create a new folder for saving
mkdir([CBF.targetP, 'ref_pix_coherence']);

%% Set some variables 
nfreqs = CBF.nfft/2+1; % The number of frequency bins
freq_res = CBF.Fs/nfreqs; % The frequency resolution
pos_min = floor(CBF.w_min / freq_res); % The corresponding power spectrum position. 

%% Calulate the pwelch Fourier Transform 

% Allocate for speed
Pxx = zeros(CBF.x,CBF.y,nfreqs);

% Calculate the PSD for every pixel
for i = 1:CBF.x
    
    % Vectorize: each column is calculated independently
    pix = squeeze(data(i,:,:)); % [pixel intensity over time]
    
    % Run the Pwelch for every pixel and save in a matrix
    [Pxx1,~] = pwelch(pix',CBF.window, CBF.noverlap, CBF.nfft, CBF.Fs, 'psd');  % [power spectral density estimate]
    Pxx(i,:,:) = Pxx1';
end


%% Calculate the cross-spectrum between one reference pixel and all other pixels

val = zeros(CBF.x, CBF.y, length(CBF.ref));

for iref = 1:size(CBF.ref,1)
    r = CBF.ref(iref,:);
    [val(:,:,iref)] = cr_coherence_ref_cpsd();
end

% Save the data 
save([CBF.targetP, CBF.name, '_Pxx.mat'], 'Pxx', 'val');


    function [val]= cr_coherence_ref_cpsd()
        
        % Specify the reference pixel
        xx = r(1);
        yy = r(2);
        ref_pix = squeeze(data(yy,xx,:)); % The xx and yy need to be inverted here!
        
        % Allocate for speed
        Pxy = zeros(CBF.x,CBF.y,nfreqs);
        
        tic % Elapsed time is 62.040316 seconds.
        for irow= 1:CBF.x
            
            test_pix = squeeze(data(irow,:,:));
            
            % Cross-spectrum
            [Pxy1,~] = cpsd(ref_pix, test_pix',CBF.window,CBF.noverlap,CBF.nfft,CBF.Fs); % [Cross-spectral density]
            Pxy(irow,:,:) = Pxy1';
            
        end
        toc
        
        %% Calculate coherence
        % According to the magnitude squared coherence:  https://se.mathworks.com/help/signal/ref/mscohere.html
        
        % Reference pixel: match the matrix to the shape of Pxx
        Pxx0 = permute(repmat(squeeze(Pxx(yy,xx,:)), 1, CBF.x, CBF.y), [2,3,1]);
        
        % Coherence
        Cxy = abs(Pxy).^2 ./ (Pxx0.*Pxx); % [unitless coherence measure]
        
        %% Plot the coherence across the nose
        
        % Allocate for speed
        val = zeros(CBF.x,CBF.y);
        pos = zeros(CBF.x,CBF.y);
        
        % Collect the coherence value at its peak value
        for irow = 1:CBF.x
            for j = 1:CBF.y
                [val(irow,j),pos(irow,j)] = max(Cxy(irow,j,pos_min:end)); % Do not look at signals below 15 Hz.
            end
        end
        
        % Adjust the positions back to the original value
        pos = pos + (pos_min-1);
        
        figure, imagesc(val);
        c = colorbar; hold on, box off, axis off
        title(sprintf('coherence for [%d,%d]', xx,yy))
        c.Label.String = 'Coherence';
        caxis([0 1])
        set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
        set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
        xlabel('\mum'), ylabel('\mum'); axis image
        
        saveas(gcf, [CBF.targetP,'ref_pix_coherence\', sprintf('coherence_pix%d_%d.png', xx,yy)]);
        
        figure, imagesc(pos);
        c = colorbar; hold on, box off, axis off
        title(sprintf('frequency at max coherence for [%d,%d]', xx,yy))
        c.Label.String = 'frequency in Hz';
        caxis([16 35]); colormap jet;
        set(gca,'XTickLabel',get(gca,'XTick')*CBF.spatres)
        set(gca,'YTickLabel',get(gca,'YTick')*CBF.spatres)
        xlabel('\mum'), ylabel('\mum'); axis image
        
        
        saveas(gcf, [CBF.targetP,'ref_pix_coherence\', sprintf('coherence_freq_pix%d_%d.png', xx,yy)]);
    end

end
