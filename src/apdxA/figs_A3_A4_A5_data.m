%--- Description ---%
%
% Filename: figs_A3_A4_A5_data.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates the data for Figures A.3, A.4 and A.5
% 
% Inputs:
% fig_num - figure number (either 3, 4 or 5)
% row_num - row number (either 1 or 2)
% col_num - column number (either 1 or 2)

function figs_A3_A4_A5_data(fig_num,row_num,col_num)

%%% Define main parameters %%%

space = ' ';

index_type = 'HC'; % use hyperbolic cross index set

err_grid_ratio = 10; % ratio between maximum index set size and error grid size
num_trials = 50; % number of random trials

num_solvers = 5; % wBP, wQCBP oracle, wQCBP CV, wsqrt-LASSO, wsqrt-LASSO CV

% CVX options
cvx_opt.precision = 'default';
cvx_opt.solver = 'mosek';
cvx_opt.verbose = false;

fig_name = ['fig_A',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num)];

% set the function
if fig_num == 3
    func = @iso_exp;
    
elseif fig_num == 4
    func = @product_cos;
    
elseif fig_num == 5
    func = @genz_gauss;
    
else
    error('invalid row number');
    
end

% set the polynomial basis and sampling measure
switch row_num
    case 1
        poly_type = 'legendre';
        samp_type = 'uniform';
    case 2
        poly_type = 'chebyshev';
        samp_type = 'chebyshev';
    otherwise
        error('invalid figure selected');
end

% set noise level (l2-norm of additive noise)
switch col_num
    case 1
        noise_level = 0;
    case 2
        noise_level = 1e-2;
    otherwise
        error('invalid figure selected');
end

% set the dimension, maximum desired number of basis functions and the values of m
d = 8;
N_max = 2000;
m_values = 50:50:500;
num_m = length(m_values);

% Set cross validation parameters (see Algorithm 1 in the book)
K = 5; % Number of groups in the K-fold cross-validation
T = 1; % Number of random tests
P = 10.^(linspace(-6,0,25)); % Defines the set of parameters to test in cv.

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
L2_error_data = zeros(num_m,num_trials,num_solvers);
Linf_error_data = zeros(num_m,num_trials,num_solvers);

%%% Main loop %%%

for i = 1:num_m
    m = m_values(i);
    
    % temporary arrays to allow parfor
    L2_error_single = zeros(num_trials,num_solvers);
    Linf_error_single = zeros(num_trials,num_solvers);
    
    u = generate_intrinsic_weights(poly_type,I); % generate the intrinsic weights
    
    parfor t = 1:num_trials
        
        y_grid = generate_sampling_grid(samp_type,d,m); % generate sample points
        A = generate_measurement_matrix(poly_type,I,y_grid); % generate measurement matrix
        b = func(y_grid)/sqrt(m); % generate measurement vector
        
        % add noise (if any)
        e = randn(size(b));  % gaussian noise
        e = noise_level * e / norm(e);  % rescale noise
        b = b + e;  % add noise to measurements
        
        
        % temporary array for parfor
        L2_error_temp = zeros(1,num_solvers);
        Linf_error_temp = zeros(1,num_solvers);
        
        % oracle eta (wQCBP) and theoretical lambda (wsqrt-LASSO)
        eta_oracle = norm(A*c_ref-b);
        lambda_theory = 1/sqrt(25*m);
        
        for j = 1:num_solvers
            switch j
                case 1 % wBP
                    
                    [c, ~] = wqcbp_cvx(A,b,u,0,cvx_opt);
                    
                case 2 % wQCBP oracle
                    
                    [c, ~] = wqcbp_cvx(A,b,u,eta_oracle,cvx_opt);
                    
                case 3 % wQCBP CV
                    
                    % cross validation
                    Delta = @(A,b,p) wqcbp_cvx(A,b,u,p,cvx_opt); % decoder
                    eta_cv = crossvalidation(A,b,K,Delta,P,T);
                    
                    % recovery via wQCBP
                    [c, ~] = Delta(A,b,eta_cv);
                    
                case 4 % wsqrt-LASSO theory
                    
                    [c, ~] = wsrlasso_cvx(A,b,u,lambda_theory,cvx_opt);
                    
                case 5 % wsqrt-LASSO CV
                    
                    % cross validation
                    Delta = @(A,b,p) wsrlasso_cvx(A,b,u,p,cvx_opt); % decoder
                    lambda_cv = crossvalidation(A,b,K,Delta,P,T);
                    
                    % recovery via wsqrt-LASSO
                    [c, ~] = Delta(A,b,lambda_cv);
                    
            end
            
            L2_err = norm(A_err_grid*c - b_err_grid)/norm(b_err_grid); % compute L^2_rho-norm error
            Linf_err = norm(A_err_grid*c - b_err_grid,Inf)/norm(b_err_grid,Inf); % compute L^2_rho-norm error
            
            L2_error_temp(j) = L2_err;
            Linf_error_temp(j) = Linf_err;
            
            disp(['Figure A.',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num),space,'m = ',num2str(m),space,'trial = ',num2str(t),space,'solver = ',num2str(j),space,'L2 error = ',num2str(L2_err),space,'Linf error = ',num2str(Linf_err)]);
        end
        
        L2_error_single(t,:) = L2_error_temp;
        Linf_error_single(t,:) = Linf_error_temp;
    end
    
    L2_error_data(i,:,:) = L2_error_single;
    Linf_error_data(i,:,:) = Linf_error_single;
end

%%% Save data %%%
clear A_err_grid
save(['../../data/apdxA/',fig_name,'_data.mat'])

end