%--- Description ---%
%
% Filename: fig_52_data.m 
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates and saves the data for Figure 5.2
% 
% Input: 
% col_num - column number (either 1 or 2)

function fig_52_data(col_num)

%%% Define main parameters %%%

space = ' ';

poly_type = 'legendre'; % use Legendre polynomials

r = 16; % largest index set size is at most 2^r

d_values = [2 4 8 16]; % values of d
num_d = length(d_values);

fig_name = ['fig_52_',num2str(col_num)];

switch col_num
    case 1
        index_type = 'TD'; % use total degree index sets
    case 2
        index_type = 'HC'; % use hyperbolic cross index sets
    otherwise
        error('invalid number selected');
end

%%% Main loop %%

% arrays for storing the data
kappa_values_data = zeros(r,num_d);
s_values_data = zeros(r,num_d);

for j = 1:num_d
    d = d_values(j); % problem dimension
    
    for i = 1:r
        
        s_des = 2^i; % desired size s
        n = find_order(index_type,d,s_des); % find the largest allowed order n
        
        I = generate_index_set(index_type,d,n);
        s = size(I,2); % actual s
        s_values_data(i,j) = s;
        
        % compute corresponding value of kappa
        u = generate_intrinsic_weights(poly_type,I);
        kappa = sum(u.^2);
        kappa_values_data(i,j) = kappa;
        
        disp(['Figure 5.2_',num2str(col_num),space,'d = ',num2str(d),space,'desired s = ',num2str(s_des),space,'n = ',num2str(n),space,'actual s = ',num2str(s),space,'kappa = ',num2str(kappa),space,'s^2 = ',num2str(s^2)]);
        
    end
end

%%% Save data %%%

save(['../../data/chpt5/',fig_name,'_data.mat']);

end
