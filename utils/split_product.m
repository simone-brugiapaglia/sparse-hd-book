%--- Description ---%
%
% Filename: split_product.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: computes the split product function f_2 defined in Appendix A
%
% Input:
% y - m x d array of sample points
%
% Output:
% b - an m x 1 array of function values at the sample points

function b = split_product(y)

[m,d] = size(y);

b = 1;

for k = 1:ceil(d/2)
    b = b./(1 - y(:,k)/4^k);
end

for k = (ceil(d/2)+1):d
    b = b.*cos(16*y(:,k)/2^k);
end

end