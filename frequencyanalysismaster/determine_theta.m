function [theta, xfield, yfield] = determine_theta(phase, kernel_sigma)
% This function determines the direction angles THETA of gradients determined
% from a phase map. It takes care of the circular nature of phases and
% allows for a spatial smoothing before determining the gradients using a
% Gaussian 3x3 kernel. The strength of the smoothing is determined by
% KERNEL_SIGMA which can vary between 0 (no smoothing) and approximately 2 
% (much smoothing). Make sure that phases are defined between 0 and 2pi.

    %% Smooth phases locally by applying a 3x3 kernel
    if kernel_sigma > 0
        h = fspecial('gaussian', 3, kernel_sigma);
    else 
        h = fspecial('gaussian', 3, eps); % this is the identity kernel
        % The identity kernel will make sure that we can compare results across
        % kernels of different sigma values. NaN values of the noise region
        % will also erode the phase map as is the case for the convolutions for
        % kernel_sigma>0. This way we make sure that we really can compare
        % results.
    end
    % This convolves the phases and maps them to the [-pi,pi] interval
    phase = angle(imfilter(exp(1i*phase), h));
    % Make sure that we map phases back to the interval [0, 2pi] to calculate
    % the gradients correctly.
    phase = phase+pi;


    %% Determine gradients.
    % Deal with the circular nature of the phase. We don't want to unwrap it
    % here in space.
    % Strategy: Calculate the gradient for both, the phase and the phase+pi.
    %           Then, for each pixel, grab the smaller of the two gradients.
    % Important: I assume here the input phases 
    %            to be in the interval [0, 2pi]. Take care of that!
    [xfield1, yfield1] = imgradientxy(phase,'prewitt');
    [xfield2, yfield2] = imgradientxy(mod(phase+pi,2*pi),'prewitt');

    magnitude = @(x,y) sqrt(x.^2 + y.^2);
    [mag1, mag2] = deal(magnitude(xfield1,yfield1), magnitude(xfield2,yfield2));
    mask = mag1 <= mag2;

    % Construct gradient field
    [xfield, yfield] = deal(nan(size(phase)));

    xfield(mask) = (1/6).*xfield1(mask); % prewitt filter needs a normalization by (1/6)
    yfield(mask) = (1/6).*yfield1(mask);

    xfield(~mask) = (1/6).*xfield2(~mask);
    yfield(~mask) = (1/6).*yfield2(~mask);

    % Create 2D support points (i.e. a grid) where the velocity vectors are
    % located.
    [X,Y] = meshgrid(1:size(xfield,2),1:size(xfield,1)); % yes, arguments are swapped because of matlab: https://stackoverflow.com/questions/28418095

    % Determine directions
    direction = @(x,y) atan2(y,x);
    theta = direction(xfield, yfield);

end 