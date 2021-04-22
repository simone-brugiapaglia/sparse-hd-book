%--- Description ---%
%
% Filename: generate_measurement_matrix.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates a random sampling grid using either the uniform or Chebyshev (arcsine) measure
%
% Inputs:
% poly_type - either 'legendre' (Legendre polynomials), 'chebyshev'
% (Chebyshev polynomials) or 'preconditioned' (preconditioned Legendre
% polynomials)
% I - s x d array of multi-indices
% y_grid - m x d array of sample points
%
% Output:
% A - normalized measurement matrix A

function A = generate_measurement_matrix(poly_type,I,y_grid)

[d,N] = size(I); % get N (number of matrix columns) and d (dimension)
m = size(y_grid,1); % get m (number of matrix rows)
A = zeros(m,N); % initialize A

n = max(I(:)); % find maximum polynomial degree

for i = 1:m
    y = y_grid(i,:); % select ith sample point
    
    % evaluate the 1D polynomials at the components of y
    if isequal(poly_type,'legendre') || isequal(poly_type,'preconditioned')
        L = legmat(y',n+1);
    elseif isequal(poly_type,'chebyshev')
        L = chebmat(y',n+1);
    else
        error('invalid poly_type')
    end
    
    % evaluate the dD polynomials via tensor products
    for j = 1:N
        Lij = zeros(d,1);
        for k = 1:d
            Lij(k,1) = L(k,I(k,j)+1);
        end
        A(i,j) = prod(Lij);
    end
end

% normalize A
A = A/sqrt(m);

% scale rows of A in preconditioned case
if isequal(poly_type,'preconditioned')
    A = sqrt(prec_w_fun(y_grid)).*A;
end

end

% compute the 1D Legendre matrix
function A = legmat(grid,k)
A = zeros(length(grid),k);
A(:,1) = 1;
A(:,2) = grid*sqrt(3);
for i = 2:k-1
    A(:,i+1)=(grid.* (2*i - 1).*A(:,i)./sqrt(i-1/2)- (i - 1).*A(:,i-1)./sqrt(i-3/2)).*sqrt(i+1/2)/i;
end
end

% compute the 1D Chebyshev matrix
function A=chebmat(grid,k)
A = zeros(length(grid),k);

for i = 1:k
    if i == 1
        A(:,i) = 1;
    else
        A(:,i) = sqrt(2).*cos((i-1).*acos(grid));
    end
end
end
