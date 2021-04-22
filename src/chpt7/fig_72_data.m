%--- Description ---%
%
% Filename: fig_72_data.m 
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates and saves the data for Figure 7.2
% 
% Inputs:
% row_num - row number (either 1, 2 or 3)
% col_num - column number (either 1 or 2)

function fig_72_data(row_num,col_num)

%%% Define main parameters %%%

space = ' ';

index_type = 'HC'; % use hyperbolic cross index set

err_grid_ratio = 10; % ratio between maximum index set size and error grid size
num_trials = 50; % number of random trials

theta_values = [0 1]; % exponents for the weights w_n = (u_n)^theta
num_theta = length(theta_values);

fig_name = ['fig_72_',num2str(row_num),'_',num2str(col_num)];

% CVX options
cvx_opt.precision = 'default';
cvx_opt.solver = 'mosek';
cvx_opt.verbose = false;

% set the dimension, function, polynomial basis and sampling measure
switch row_num
    case 1
        d = 4;
        func = @iso_exp;
        N_max = 1000;
        m_values = 25:25:250;
    case 2
        d = 8;
        func = @product_cos;
        N_max = 2000;
        m_values = 50:50:500;
    case 3
        d = 16;
        func = @genz_gauss;
        N_max = 4500;
        m_values = 100:100:1000;
    otherwise
        error('invalid row number');
end

switch col_num
    case 1
        poly_type = 'legendre'; % use Legendre polynomials
        samp_type = 'uniform'; % use samples drawn randomly from the uniform distribution
    case 2
        poly_type = 'chebyshev'; % use Chebyshev polynomials
        samp_type = 'chebyshev'; % use samples drawn randomly from the Chebyshev distribution
    otherwise
        error('invalid column number');
end

num_m = length(m_values);

%%% Construct index set and define N %%%

n = find_order(index_type,d,N_max); % find the maximum polynomial order for the given index_type
I = generate_index_set(index_type,d,n); % compute index set
N = size(I,2); % number of basis functions

u = generate_intrinsic_weights(poly_type,I); % generate the intrinsic weights

%%% Construct error grid, error matrix and error vector %%%

M = ceil(err_grid_ratio*N);
err_grid = generate_sampling_grid(samp_type,d,M);
A_err_grid = generate_measurement_matrix(poly_type,I,err_grid);
b_err_grid = func(err_grid)/sqrt(M);

%%% Compute a reference solution via least squares on the error grid %%%

c_ref = A_err_grid\b_err_grid;

% arrays for storing the errors
L2_error_wqcbp_data = zeros(num_m,num_trials,num_theta);
Linf_error_wqcbp_data = zeros(num_m,num_trials,num_theta);
L2_error_wsrlasso_data = zeros(num_m,num_trials,num_theta);
Linf_error_wsrlasso_data = zeros(num_m,num_trials,num_theta);

%%% Main loop %%%

for i = 1:num_m
    m = m_values(i);
    lambda = 1/sqrt(25*m); % value of lambda for wsrlasso
    
    % temporary arrays to allow parfor
    L2_error_wqcbp_single = zeros(num_trials,num_theta);
    Linf_error_wqcbp_single = zeros(num_trials,num_theta);
    L2_error_wsrlasso_single = zeros(num_trials,num_theta);
    Linf_error_wsrlasso_single = zeros(num_trials,num_theta);
    
    parfor t = 1:num_trials
        
        y_grid = generate_sampling_grid(samp_type,d,m); % generate sample points
        A = generate_measurement_matrix(poly_type,I,y_grid); % generate measurement matrix
        b = func(y_grid)/sqrt(m); % generate measurement vector
        
        eta = norm(A*c_ref - b); % use reference solution to find best eta
        
        % temporary arrays to allow parfor
        L2_error_wqcbp_temp = zeros(1,num_theta);
        Linf_error_wqcbp_temp = zeros(1,num_theta);
        L2_error_wsrlasso_temp = zeros(1,num_theta);
        Linf_error_wsrlasso_temp = zeros(1,num_theta);
        
        for j = 1:num_theta
            
            theta = theta_values(j); w = u.^theta; % define weights
            
            % loop over the two solvers
            for soltype = 1:2
                
                if soltype == 1
                    solver_type = 'wqcbp';
                    [c,stat] = wqcbp_cvx(A,b,w,eta,cvx_opt);
                else
                    solver_type = 'wsrlasso';
                    [c,stat] = wsrlasso_cvx(A,b,w,lambda,cvx_opt);
                end
                
                L2_err = norm(A_err_grid*c - b_err_grid)/norm(b_err_grid); % compute L^2_rho-norm error
                Linf_err = norm(A_err_grid*c - b_err_grid,Inf)/norm(b_err_grid,Inf); % compute L^2_rho-norm error
                
                if soltype == 1
                    L2_error_wqcbp_temp(j) = L2_err;
                    Linf_error_wqcbp_temp(j) = Linf_err;
                else
                    L2_error_wsrlasso_temp(j) = L2_err;
                    Linf_error_wsrlasso_temp(j) = Linf_err;
                end
                
                disp(['Figure 7.1_',num2str(row_num),'_',num2str(col_num),space,'solver = ',solver_type,space,'m = ',num2str(m),space,'trial ',num2str(t),' of ',num2str(num_trials),space,'theta = ',num2str(theta),space,'L2 error = ',num2str(L2_err),space,'Linf error = ',num2str(Linf_err)]);
            end
        end
        
        L2_error_wqcbp_single(t,:) = L2_error_wqcbp_temp;
        Linf_error_wqcbp_single(t,:) = Linf_error_wqcbp_temp;
        L2_error_wsrlasso_single(t,:) = L2_error_wsrlasso_temp;
        Linf_error_wsrlasso_single(t,:) = Linf_error_wsrlasso_temp;
    end
    
    L2_error_wqcbp_data(i,:,:) = L2_error_wqcbp_single;
    Linf_error_wqcbp_data(i,:,:) = Linf_error_wqcbp_single;
    L2_error_wsrlasso_data(i,:,:) = L2_error_wsrlasso_single;
    Linf_error_wsrlasso_data(i,:,:) = Linf_error_wsrlasso_single;
end

%%% Save data %%%
clear A_err_grid
save(['../../data/chpt7/',fig_name,'_data.mat'])

end