%--- Description ---%
%
% Filename: figs_73_74_data.m 
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates and saves the data for Figures 7.3 and 7.4
% 
% Inputs:
% fig_num - figure number (either 3 or 4)
% row_num - row number (either 1, 2 or 3)
% col_num - column number (either 1 or 2)

function figs_73_74_data(fig_num,row_num,col_num)

%%% Define main parameters %%%

space = ' ';

index_type = 'HC'; % use hyperbolic cross index set

err_grid_ratio = 10; % ratio between maximum index set size and error grid size
num_trials = 50; % number of random trials

theta_values = [0 0.5 1 1.5 2 3 4]; % exponents for the weights w_n = (u_n)^theta
num_theta = length(theta_values);

% CVX options
cvx_opt.precision = 'default';
cvx_opt.solver = 'mosek';
cvx_opt.verbose = false;

fig_name = ['fig_7',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num)];

% set the polynomial basis and sampling measure
if fig_num == 3
    poly_type = 'legendre';
    samp_type = 'uniform';
    
elseif fig_num == 4
    poly_type = 'chebyshev';
    samp_type = 'chebyshev';
    
else
    error('invalid figure selected');
end

% set the dimension, function, maximum desired number of basis functions and the values of m
if row_num == 1
    func = @iso_exp;
    
elseif row_num == 2
    func = @split_product;
    
elseif row_num == 3
    func = @genz_gauss;
    
else
    error('invalid row number');
    
end

if col_num == 1
    d = 8;
    N_max = 2000;
    m_values = 50:50:500;
    
elseif col_num == 2
    d = 16;
    N_max = 4500;
    m_values = 100:100:1000;
    
else
    error('invalid column number');
    
end

num_m = length(m_values);

%%% Construct index set and define N %%%

n = find_order(index_type,d,N_max); % find the maximum polynomial order for the given index_type
I = generate_index_set('HC',d,n); % compute index set
N = size(I,2); % number of basis functions

%%% Construct error grid, error matrix and error vector %%%

M = ceil(err_grid_ratio*N);
err_grid = generate_sampling_grid(samp_type,d,M);
A_err_grid = generate_measurement_matrix(poly_type,I,err_grid);
b_err_grid = func(err_grid)/sqrt(M);

%%% Compute a reference solution via least squares on the error grid %%%

c_ref = A_err_grid\b_err_grid;

% arrays for storing the errors
L2_error_data = zeros(num_m,num_trials,num_theta);
Linf_error_data = zeros(num_m,num_trials,num_theta);

%%% Main loop %%%

for i = 1:num_m
    m = m_values(i);
    
    % temporary arrays to allow parfor
    cond_num_single = zeros(num_trials,num_theta);
    L2_error_single = zeros(num_trials,num_theta);
    Linf_error_single = zeros(num_trials,num_theta);
    
    u = generate_intrinsic_weights(poly_type,I); % generate the intrinsic weights
    
    parfor t = 1:num_trials
        
        y_grid = generate_sampling_grid(samp_type,d,m); % generate sample points
        A = generate_measurement_matrix(poly_type,I,y_grid); % generate measurement matrix
        b = func(y_grid)/sqrt(m); % generate measurement vector
        
        eta = norm(A*c_ref-b); % use reference solution to find best eta
        
        % temporary array for parfor
        L2_error_temp = zeros(1,num_theta);
        Linf_error_temp = zeros(1,num_theta);
        
        for j = 1:num_theta
            
            theta = theta_values(j); w = u.^theta; % define weights
            
            [c,stat] = wqcbp_cvx(A,b,w,eta,cvx_opt); % solve the wQCBP problem
            
            L2_err = norm(A_err_grid*c - b_err_grid)/norm(b_err_grid); % compute L^2_rho-norm error
            Linf_err = norm(A_err_grid*c - b_err_grid,Inf)/norm(b_err_grid,Inf); % compute L^2_rho-norm error
            
            L2_error_temp(j) = L2_err;
            Linf_error_temp(j) = Linf_err;
            
            disp(['Figure 7.',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num),space,'m = ',num2str(m),space,'trial = ',num2str(t),space,'theta = ',num2str(theta),space,'L2 error = ',num2str(L2_err),space,'Linf error = ',num2str(Linf_err)]);
        end
        
        L2_error_single(t,:) = L2_error_temp;
        Linf_error_single(t,:) = Linf_error_temp;
    end
    
    L2_error_data(i,:,:) = L2_error_single;
    Linf_error_data(i,:,:) = Linf_error_single;
end

%%% Save data %%%
clear A_err_grid
save(['../../data/chpt7/',fig_name,'_data.mat'])

end