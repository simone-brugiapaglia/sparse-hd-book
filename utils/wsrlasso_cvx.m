% Filename: wsrlasso_cvx.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: solves the (weighted) SR-LASSO problem using CVX
%
% Inputs:
% A - m x N measurement matrix
% b - m x 1 measurement vector
% w - N x 1 vector of weights (optional)
% lambda - SR-LASSO parameter
% cvx_opt - options for CVX
%
% Outputs:
% c - N x 1 array, the solution to the SR-LASSO problem
% stat - CVX solve status
%
% Note: if w is set as [], then the unweighted SR-LASSO problem is solved

function [c,stat] = wsrlasso_cvx(A,b,w,lambda,cvx_opt)

N = size(A,2);

if isequal(w,[])
    w = ones(N,1);
end

if ~exist('cvx_opt','var')
    cvx_opt = [];
end
cvx_opt = set_options(cvx_opt);

w = w(:); % the weights must form a column vector

if lambda >= 0
    cvx_begin
        cvx_quiet(~cvx_opt.verbose)
        cvx_solver(cvx_opt.solver)
        cvx_precision(cvx_opt.precision)
        
        variable z(N)
        minimize(lambda*norm(w.*z,1)+norm(A*z-b,2))
    cvx_end
    
    stat = cvx_status;
    c = z;
else
    error('LAMBDA must be nonnegative')
end
