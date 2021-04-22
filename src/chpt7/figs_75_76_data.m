%--- Description ---%
%
% Filename: figs_75_76_data.m 
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates and saves the data for Figures 7.5 and 7.6
% 
% Inputs:
% fig_num - figure number (either 5 or 6)
% row_num - row number (either 1 or 2)
% col_num - column number (either 1 or 2)

function figs_75_76_data(fig_num,row_num,col_num)

%%% Define main parameters %%%

space = ' ';

index_type = 'HC'; % use hyperbolic cross index set

err_grid_ratio = 10; % ratio between maximum index set size and error grid size
num_trials = 50; % number of random trials

% CVX options
cvx_opt.precision = 'default';
cvx_opt.solver = 'mosek';
cvx_opt.verbose = false;

fig_name = ['fig_7',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num)];

% set the function
if fig_num == 5
    func = @iso_exp;
elseif fig_num == 6
    func = @split_product; 
else
    error('invalid figure selected');
end

% set the dimension, maximum desired number of basis functions and the values of m
if row_num == 1 && col_num == 1
    d = 1;
    N_max = 400;
    m_values = 10:10:100;
    
elseif row_num == 1 && col_num == 2
    d = 2;
    N_max = 800;
    m_values = 20:20:200;
    
elseif row_num == 2 && col_num == 1
    d = 8;
    N_max = 2000;
    m_values = 50:50:500;
    
elseif row_num == 2 && col_num == 2
    d = 16;
    N_max = 4500;
    m_values = 100:100:1000;
    
else
    error('invalid subfigure selected');
end

num_m = length(m_values);

%%% Construct index set and define N %%%

n = find_order(index_type,d,N_max); % find the maximum polynomial order for the given index_type
I = generate_index_set('HC',d,n); % compute index set
N = size(I,2); % number of basis functions

% arrays for storing the errors
num_cases = 3;
L2_error_data = zeros(num_m,num_trials,num_cases);
Linf_error_data = zeros(num_m,num_trials,num_cases);

%%% Construct error grid %%%
    
err_samp_type = 'uniform'; % sampling points for the error grid
M = ceil(err_grid_ratio*N);
err_grid = generate_sampling_grid(err_samp_type,d,M);

%%% Main loop %%%

for l = 1:num_cases
    
    % Legendre, uniform
    if l == 1
        samp_type = 'uniform';
        basis_type = 'legendre'; % functions used in the measurement matrix A
        err_basis_type = basis_type; % polynomials used in computing the error
        
        % Chebyshev, Chebyshev
    elseif l == 2
        samp_type = 'chebyshev';
        basis_type = 'chebyshev';
        err_basis_type = basis_type;
        
        % Preconditioned Legendre, Chebyshev
    else
        samp_type = 'chebyshev';
        basis_type = 'preconditioned';
        err_basis_type = 'legendre';
    end
    
    %%% Construct error matrix and error vector %%%
    
    A_err_grid = generate_measurement_matrix(err_basis_type,I,err_grid);
    b_err_grid = func(err_grid)/sqrt(M);
    
    %%% Compute a reference solution via least squares on the error grid %%%
    
    c_ref = A_err_grid\b_err_grid;
    
    for i = 1:num_m
        m = m_values(i);
        
        % temporary arrays to allow parfor
        L2_error_temp = zeros(1,num_trials);
        Linf_error_temp = zeros(1,num_trials);
        
        parfor t = 1:num_trials
            
            y_grid = generate_sampling_grid(samp_type,d,m); % generate sample points
            A = generate_measurement_matrix(basis_type,I,y_grid); % generate measurement matrix
            
            % generate measurement vector
            if isequal(basis_type,'preconditioned')
                b = sqrt(prec_w_fun(y_grid)).*func(y_grid)/sqrt(m);
            else
                b = func(y_grid)/sqrt(m);
            end
            
            eta = norm(A*c_ref-b); % use reference solution to find best eta
            
            [c,stat] = wqcbp_cvx(A,b,[],eta,cvx_opt); % solve the weighted QCBP problem
            
            L2_err = norm(A_err_grid*c - b_err_grid)/norm(b_err_grid); % compute L^2_rho-norm error
            Linf_err = norm(A_err_grid*c - b_err_grid,Inf)/norm(b_err_grid,Inf); % compute L^2_rho-norm error
            
            disp(['Figure 7.',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num),space,'polynomials = ',basis_type,space,'m = ',num2str(m),space,'trial = ',num2str(t),space,'L2 error = ',num2str(L2_err),space,'Linf error = ',num2str(Linf_err)]);
            
            L2_error_temp(1,t) = L2_err;
            Linf_error_temp(1,t) = Linf_err;
        end    
        
    L2_error_data(i,:,l) = L2_error_temp;
    Linf_error_data(i,:,l) = Linf_error_temp;
        
    end
     
end

%%% Save data %%%
clear A_err_grid
save(['../../data/chpt7/',fig_name,'_data.mat'])

end