function [theta,THL] = est_theta_Fourier(b)
% Estimate main texture direction in an image.
%
% Input
% b:        square corrupted image, e.g. noisy or blurry (N x N)
%
% Output
% theta:    main texture direction
%
% By Rasmus Dalgas Kongskov, 18/12/2016, DTU Compute

N       = size(b); % size of image (assuming square image)
w       = fft2(b); % 2D Fourier transform

% modulus of coefficients and removing of constant region coefficients
aw      = abs(w(1:end,1:round(N/2)));
aw(1,1) = 0;

% find direction related to the ten maximum magnitude coefficients
THL = zeros(10,1);
for k = 1:10

    % find maximum magnitude coefficient
    [i,j]   = find(aw==max(aw(:)),1);
    
    % find frequency center
    if i>N/2
        c = [1,N+1];
    else
        c = [1,1];
    end
    
    % calculate main direction of image texture
    THL(k) = atan((j-c(1))/(i-c(2)))/pi*180;
    
    aw(i,j)=0;
    
end

% pick main direction from the median
theta = median(THL);

end