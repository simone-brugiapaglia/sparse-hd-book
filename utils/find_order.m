%--- Description ---%
%
% Filename: find_order.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: computes the order of the largest multi-index set of a given
% type within a specific maximum size N
%
% Inputs:
% index_type - either 'TP' (tensor product), 'TD' (total degree) or 'HC' (hyperbolic cross)
% d - dimension
% N - maximum desired size
%
% Output:
% n - the largest polynomial order so that the index set has size at most N

function n = find_order(index_type,d,N)

n = 1;
M = index_size(index_type,d,n);

% check the value n = 1 and return n = 0 if it exceeds N
if M > N
    n = 0;
    return;
end

% double n until index set size exceeds N
while M <= N
    
    M = index_size(index_type,d,n);
    n = 2*n;
end
n = n/2;

% bisection search for n until value found
nu = n;
nl = max(floor(n/2),1);
flg = 1;

while nu - nl > 1
    
    n = max(floor((nl+nu)/2),1);
    M = index_size(index_type,d,n);
    
    if M > N
        nu = n; flg = 1;
    else
        nl = n; flg = 0;
    end
    
end

if flg == 1
    n = n-1;
end

end

% determine the size of a given type of index set
function K = index_size(index_type,d,n)

if isequal(index_type,'TP')
    K = (n+1)^d;
    
elseif isequal(index_type,'TD')
    K = nchoosek(n+d,d);
    
elseif isequal(index_type,'HC')
    K = size(generate_index_set('HC',d,n),2);
    
else
    error('invalid index_type');
end

end

