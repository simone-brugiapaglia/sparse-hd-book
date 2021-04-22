%--- Description ---%
%
% Filename: generate_intrinsic_weights.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates the intrinsic weights for polynomial approximation
%
% Inputs:
% poly_type - either 'legendre' (Legendre polynomials), 'chebyshev'
% (Chebyshev polynomials) or 'preconditioned' (precoditioned Legendre
% polynomials)
% I - N x d array of multi-indices
%
% Output:
% u - N x 1 array of weights

function u = generate_intrinsic_weights(poly_type,I)

[d,N] = size(I); % get N (number of matrix columns) and d (dimension)
u = zeros(1,N);

if isequal(poly_type,'legendre')
    
    if d == 1
        u = sqrt(2.*I+ones(size(I)));      
    else
        u = prod(sqrt(2.*I+ones(size(I))));
    end
    
elseif isequal(poly_type,'chebyshev')
    
    for j = 1:N
        u(1,j) = sqrt(2)^(nnz(I(:,j)));
    end
       
elseif isequal(poly_type,'preconditioned')
    
    n = max(I(:)); % find max polynomial degree in I
    M = 100*n; y_grid = (linspace(-1,1,M))'; % generate a fine grid
    
    % compute the 1D preconditioned polynomials on the grid and find their
    % maxima
    A = sqrt(M)*generate_measurement_matrix('preconditioned',0:n,y_grid); 
    v = max(abs(A),[],1);
    
    % tensorize to get weights
    for j = 1:N
        u(1,j) = prod(v(I(:,j)+1));
    end
    
else
    error('Invalid poly_type');
end

u = u';
end


