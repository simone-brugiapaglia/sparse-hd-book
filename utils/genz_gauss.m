%--- Description ---%
%
% Filename: genz_gauss.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: computes the Gaussian function f_3 defined in Appendix A
%
% Inputs:
% y - m x d array of sample points
%
% Output:
% b - m x 1 array of function values at the sample points

function b = genz_gauss(y)

[m,d] = size(y);

b = 0;
for k = 1:d
b = b + (y(:,k) - (-1)^k/(k+1)).^2;
end
b = exp(-2*b/d);

end