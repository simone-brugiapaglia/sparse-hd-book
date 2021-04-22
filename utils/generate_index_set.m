%--- Description ---%
%
% Filename: generate_index_set.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: computes the (isotropic) tensor product, total degree or hyperbolic cross multi-index set
%
% Inputs:
% index_type - either 'TP' (tensor product), 'TD' (total degree) or 'HC' (hyperbolic cross)
% d - dimension
% n - polynomial order
%
% Output:
% I - the d x n array where the columns are the multi-indices from the
% desired multi-index set

function I = generate_index_set(index_type,d,n)

I = 0:n; % initialize 1D array

if d >= 2
    
    for k = 2:d
        J = [];
        
        for i = 0:n
            l = size(I,2);
            
            for j = 1:l
                z = I(:,j); % get the jth index in I
                
                % test whether to add [z ; i] to the new index set
                if index_type == 'TP'
                    crit = i <= n;
                elseif index_type == 'TD'
                    crit = i + sum(z) <= n;
                elseif index_type == 'HC'
                    crit = (i+1)*prod(z+1) <= n+1;
                else
                    error('invalid index_type');
                end
                
                % add [z ; i]
                if crit == 1
                    z = [z ; i ];
                    J = [J z];
                end
                
            end
        end
        
        I = J;
    end
    
end

end

