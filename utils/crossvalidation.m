%--- Description ---%
%
% Filename: crossvalidation.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: G-fold cross validation. The procedure corresponds to that
% described in Algorithm 1 in Appendix A.
%
% Inputs:
% A - Design matrix
% y - Vector of samples
% G - Number of groups
% Delta - Decoder; a function handle @(A,y,p)
% P - Finite set of parameters (array)
% T - Number of random tests
%
% Outputs:
% p_cv - Estimated parameter
% epsilon - Vector of validation errors

function [p_cv,epsilon] = crossvalidation(A,y,G,Delta,P,T)

[m,~] = size(A);
if length(y) ~= m
    error('A and y should have the same number of rows')
end

epsilon = zeros(T,G,length(P));

for t = 1:T
    
    %% Partition [m] in sets of cardinality floor(m/G) or floor(m/G)+1
    % Thanks to the equality
    %
    % (G-R)*(floor(m/G)) + R * (floor(m/G)+1),   with R = rem(m,G)
    %
    % [m] contains (G-R) groups of cardinality floor(m/G) and R of
    % cardinality floor(m/G)+1.
    R = rem(m,G);
    perm = randperm(m);
    group_start = 1;
    
    for g = 1 : G-R
        group_end = group_start + floor(m/G) - 1;
        I{g} = perm(group_start : group_end);
        group_start = group_end + 1;
    end
    
    for g = G-R+1 : G
        group_end = group_start + floor(m/G);
        I{g} = perm(group_start : group_end);
        group_start = group_end + 1;
    end
    
    for g = 1:G
        %% I{g} as validation set
        mv = length(I{g}); % # of validation samples
        Icomp = setdiff((1:m),I{g}); % complement of I{g}
        % Validation samples
        Av = sqrt(m/mv) * A(I{g},:);
        yv = sqrt(m/mv) * y(I{g});
        % Reconstruction samples
        Ar = sqrt(m/(m-mv)) * A(Icomp,:);
        yr = sqrt(m/(m-mv)) * y(Icomp);
        
        for i_p = 1:length(P)
            %% Reconstruction and validation with parameter p
            p = P(i_p);
            c = Delta(Ar,yr,p); % Reconstruction
            epsilon(t,g,i_p) = norm(Av * c - yv, 2)^2; % Validation
        end
    end
end

fprintf('\n')

% compute average validation error
mean_epsilon = sum(sum(epsilon,1),2) / (G*T);
[~, i_min] = min(mean_epsilon);

p_cv = P(i_min); % estimated parameter

end
