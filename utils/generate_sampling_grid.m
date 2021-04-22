%--- Description ---%
%
% Filename: generate_sampling_grid.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates a random sampling grid using either the uniform or Chebyshev (arcsine) measure
%
% Inputs: 
% samp_type - either 'uniform' (uniform measure) or 'Chebyshev' (Chebyshev measure)
% d - dimension
% m - number of sample points
% 
% Output:
% y_grid - the m x d array where the ith row is the ith sample point y_i

function y_grid = generate_sampling_grid(samp_type,d,m)

% uniform measure
if isequal(samp_type,'uniform')
    
y_grid = 2*rand(m,d)-ones(m,d); 

% Chebyshev measure
elseif isequal(samp_type,'chebyshev')

    y_grid = rand(m,d);
    y_grid = cos(pi .* y_grid); 

else

error('invalid samp_type')

end