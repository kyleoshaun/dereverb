function [ mean_ssim, ssim_map ] = mssim_v2( image_1, image_2, weights, window, window_type )

% This script-file implements the structural similarity quality metric
% defined by Wang, 2004.
%
% Inputs:
%
%   image_1 - Reference image.
%   image_2 - Conditioned image.
%
%   weights - Matlab structure with fields corresponding to weight-term for
%       luminance, contrast and structure.  Given by alpha, beta and gamma, respectively.
%
%   window - Rectangular windows matrix defining the patch or area over which the
%       SSIM metric is calculated.
%
%   window_type - A logical flag determining the use of either the given rectangular
%       window (normalized) or a Gaussian window based on the largest dimension of the
%       given window argument.  A zero applies the rectagular-window.  The Gaussian
%       window is considered due to a comment made in Wang (2004) about "patchiness"
%       introduced through the use of the rectangular-window.
%
%
% Outputs:
%   ssim_map - Rectangular matrix contained results of local luminance, constrast and
%       structure calculations associated with the "patches."
%
%   mean_ssim - Average or mean of ssim_map values.
%   
%
% Notes:
%  a.)  Assumes unsigned 8-bit gray-scale for images.
%
%% See Also
%

% Author:  Michael R. Wirtzfeld
% Modification Date:  Thursday, April 26, 2012
% Creation Date:  Wednesday, April 25, 2012
% Modified by Ian Bruce to use mean(.,'all') instead of mean2(.), as the
% latter is only available with the image processing toolbox

%% Validate Arguments

IMAGE_SIZE_THRESHOLD = 5;
MINIMUM_WINDOW_SIZE = 3;

scaling_coefficients = [ 0.01 0.05 0.05 ];


% Validate number of arguments passed into function.
assert( ~(nargin ~= 5), '*** Insufficient Number of Arguments ***' );


% Test to ensure images are the same size.
if ( size(image_1) ~= size(image_2) )
    mean_ssim = NaN;
    ssim_map = NaN;
    return;
end


% Test size of given images.
[ rows, columns ] = size( image_1 );

if ( ( rows < IMAGE_SIZE_THRESHOLD ) || ( columns < IMAGE_SIZE_THRESHOLD ) )
    mean_ssim = NaN;
    ssim_map = NaN;
    return;
end


% Ensure minimum window size.
[ rows_window, columns_window ] = size( window );

if ( (rows_window * columns_window ) < MINIMUM_WINDOW_SIZE || ( rows_window > rows ) || ( columns_window > columns ) )
    mean_ssim = NaN;
    ssim_map = NaN;
    return;
end


% Ensure there are 3 coefficients for the limunance, contrast and structure terms.
assert( ~(length(scaling_coefficients) ~= 3)  , '*** Insufficient Number of Arguments ***' );


% Determine maximum bit range for scaling.
number_of_dimensions = length( size(image_1) );
    assert( ~(number_of_dimensions == 3), '*** Non-gray-scale Image ***' );
    
image_class_type = class( image_1 );
    assert( ~isa( image_class_type, 'unit8' ), '*** Invalid Class-type - Must be uint8 ***' );
    
maximum_value = 2^8 - 1;


% Type of Processing Window
if ( window_type == 0 )  % Square Window
    window = window ./ sum( window(:) );
else % Gaussian Window
    window = fspecial( 'gaussian', max(size(window)), 1.5 );
end



%% Processing

image_1 = double( image_1 );
image_2 = double( image_2 );

C_values = ( scaling_coefficients .* maximum_value ).^2;

% Mean Value Window Set
mu_1 = filter2( window, image_1, 'same' );
    mu_1_sq = mu_1.^2;
mu_2 = filter2( window, image_2, 'same' );
    mu_2_sq = mu_2.^2;
    
% Standard Deviation Value Window Set
sigma_1_sq = filter2( window, image_1.*image_1, 'same' ) - mu_1_sq;
    sigma_1 = abs( sigma_1_sq.^(0.5) );
sigma_2_sq = filter2( window, image_2.*image_2, 'same' ) - mu_2_sq;
    sigma_2 = abs( sigma_2_sq.^(0.5) );
%
sigma_12 = filter2( window, image_1.*image_2, 'same' ) - mu_1.*mu_2;


% Dimensions
luminance = ( 2 * mu_1 .* mu_2 + C_values(1) ) ./ ( mu_1_sq + mu_2_sq + C_values(1) );
    luminance = luminance .^ weights.alpha;

contrast = ( 2 * sigma_1 .* sigma_2 + C_values(2) ) ./ ( sigma_1_sq + sigma_2_sq + C_values(2) );
    contrast = contrast .^  weights.beta;

structure = ( sigma_12 + C_values(3) ) ./ ( sigma_1 .* sigma_2 + C_values(3) );
    structure = structure .^ weights.gamma;
    
    
% SSIM Map and MSSIM
ssim_map = luminance .* contrast .* structure;
mean_ssim = mean( ssim_map,'all');



%% References

% http://www.mathworks.com/help/toolbox/images/ref/blockproc.html
% http://www.mathworks.com/help/toolbox/images/f7-12726.html
% http://stackoverflow.com/questions/1637000/how-to-divide-an-image-into-blocks-in-matlab


