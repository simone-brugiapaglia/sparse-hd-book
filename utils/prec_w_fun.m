%--- Description ---%
%
% Filename: prec_w_fun.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: computes the function w(y) in the preconditioning scheme, 
% given by eqn (5.42)
%
% Input:
% y - an m x d grid of points, where each row is a sample point
%
% Output:
% b - an m x 1 array of function values at the sample points

function b = prec_w_fun(y)

[m,d] = size(y);

b = 1;
for k = 1:d
b = b*pi.*sqrt(1-y(:,k).^2);
end