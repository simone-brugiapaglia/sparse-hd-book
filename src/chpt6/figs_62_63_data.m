%--- Description ---%
%
% Filename: figs_62_63_data.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse polynomial approximation of high-dimensional functions", SIAM
%
% Description: generates the data for Figures 6.2 and 6.3
% 
% Inputs: 
% fig_num - figure number (either 1 or 2)
% row_num - row number (either 1, 2 or 3)
% col_num - row number (either 1 or 2)

function figs_62_63_data(fig_num,row_num,col_num)

%%% Define main parameters %%%

space = ' ';

index_type = 'HC'; % use hyperbolic cross index set

d = 8; % dimension
N_max = 2000; % max number of basis functions
m = 500;

func = @iso_exp; % set the function to approximate

err_grid_ratio = 10; % ratio between maximum index set size and error grid size
num_trials = 50; % number of random trials

num_solvers = 3;

% noise levels to use
sigma_values = [10^(-1) 10^(-2) 10^(-3) 0];
num_noise = length(sigma_values);

% parameter values to use
num_params = 19;

% parameter values
switch row_num
    case 1
        param_values = 10.^(-7:0.5:2);
    case 2
        param_values = 10.^(-8:0.5:1);
    case 3
        param_values = 10.^(-3:0.25:1.5);
    otherwise
        error('invalid row number');
end

%%% Set solver options %%%

% CVX options
cvx_opt.precision = 'default';
cvx_opt.solver = 'mosek';
cvx_opt.verbose = false;

fig_name = ['fig_6',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num)];

% set the polynomial basis and sampling measure
if fig_num == 1
    poly_type = 'legendre';
    samp_type = 'uniform';
    
elseif fig_num == 2
    poly_type = 'chebyshev';
    samp_type = 'chebyshev';
    
else
    error('invalid figure selected');
end

%%% Main loop %%%

%%% Construct index set and define N %%%

n = find_order(index_type,d,N_max); % find the maximum polynomial order for the given index_type
I = generate_index_set('HC',d,n); % compute index set
N = size(I,2); % number of basis functions

%%% Construct error grid, error matrix and error vector %%%

M = ceil(err_grid_ratio*N);
err_grid = generate_sampling_grid(samp_type,d,M);
A_err_grid = generate_measurement_matrix(poly_type,I,err_grid);
b_err_grid = func(err_grid)/sqrt(M);

% arrays for storing the errors
error_data = zeros(num_params,num_trials,num_noise);
    
parfor t = 1:num_trials
    
    error_temp = zeros(num_params,num_noise); % temporary array for parfor
    
    y_grid = generate_sampling_grid(samp_type,d,m); % generate sample points
    A = generate_measurement_matrix(poly_type,I,y_grid); % generate measurement matrix
    b = func(y_grid)/sqrt(m); % generate measurement vector
    
    e = 2*rand(m,1)-1; % random noise vector
    
    for j = 1:num_noise
        
        % add random noise to the measurement vector
        sigma = sigma_values(j);
        b_noisy = b + sigma*e/norm(e);
        
        for i = 1:num_params
                
            param = param_values(i);
            
            switch col_num
                case 1
                    w = [];
                case 2
                    u = generate_intrinsic_weights(poly_type,I); % generate the intrinsic weights
                    w = u;
                otherwise
                    error('invalid column number');
            end
            
            if row_num == 1
                [c,stat] = wqcbp_cvx(A,b_noisy,w,param,cvx_opt); % solve the wQCBP problem
            elseif row_num == 2
                [c,stat] = wlasso_cvx(A,b_noisy,w,param,cvx_opt); % solve the wLASSO problem
            elseif row_num == 3
                [c,stat] = wsrlasso_cvx(A,b_noisy,w,param,cvx_opt); % solve the wLASSO problem
            else
                error('invalid column number');
            end
            
            err = norm(A_err_grid*c - b_err_grid)/norm(b_err_grid); % compute L^2_rho-norm error
            error_temp(i,j) = err;
            
            disp(['Figure 6.',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num),space,'parameter case = ',num2str(i),' of ',num2str(num_params),space,'noise case = ',num2str(j),' of ',num2str(num_noise),space,'trial = ',num2str(t),' of ',num2str(num_trials),space,'L2 error = ',num2str(err)]);
        end

    end
    
    error_data(:,t,:) = error_temp;
      
end

%%% Save data %%%
clear A A_err_grid b_err_grid
save(['../../data/chpt6/',fig_name,'_data.mat'])

end