% Filename: wqcbp_cvx.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: solves the (weighted) QCBP problem using CVX
%
% Inputs:
% A - m x N measurement matrix
% b - m x 1 measurement vector
% w - N x 1 vector of weights (optional)
% eta - QCBP parameter
% cvx_opt - options for CVX
%
% Outputs:
% c - N x 1 array, the solution to the QCBP problem
% stat - CVX solve status
%
% Note: if w is set as [], then the unweighted QCBP problem is solved

function [c,stat] = wqcbp_cvx(A,b,w,eta,cvx_opt)

N = size(A,2); 

if isequal(w,[])
    w = ones(N,1);
end

w = w(:); % the weights must form a column vector

if ~exist('cvx_opt','var')
    cvx_opt = [];
end
cvx_opt = set_options(cvx_opt);

if eta < 0
    error('eta must be nonnegative.')

elseif eta == 0
    cvx_begin  
        cvx_quiet(~cvx_opt.verbose)    
        cvx_solver(cvx_opt.solver)
        cvx_precision(cvx_opt.precision)
    
        variable z(N)
        minimize(norm(w.*z,1))
        subject to
            A*z == b       
    cvx_end
     
    stat = cvx_status;
    c = z;    

else
    cvx_begin 
        cvx_quiet(~cvx_opt.verbose)    
        cvx_solver(cvx_opt.solver)
        cvx_precision(cvx_opt.precision)
        
        variable z(N)
        minimize(norm(w.*z,1))
        subject to
            norm(A*z-b,2) <= eta      
    cvx_end
    
    stat = cvx_status;
    c = z;
end