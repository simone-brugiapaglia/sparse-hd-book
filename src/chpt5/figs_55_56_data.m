%--- Description ---%
%
% Filename: figs_55_56_data.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse polynomial approximation of high-dimensional functions", SIAM
%
% Description: generates the data for Figures 5.4 and 5.5
% 
% Inputs: 
% fig_num - figure number (either 5 or 6)
% col_num - column number (either 1 or 2)

function figs_55_56_data(fig_num,row_num)

%%% Define main parameters %%%

space = ' ';

poly_type = 'legendre'; % use Legendre polynomials
samp_type = 'uniform'; % use samples drawn randomly from the uniform distribution

func = @iso_exp; % function to approximate

err_grid_ratio = 10; % ratio between maximum index set size and error grid size
num_trials = 50; % number of random trials

max_rows = 10000; % maximum number of rows in the measurement matrix
num_n_des = 20; % desired number of n values to use

fig_name = ['fig_5',num2str(fig_num),'_',num2str(row_num)];

switch fig_num
    case 5
        index_type = 'TD';
    case 6
        index_type = 'HC';
    otherwise
        error('invalid figure number selected');
end

switch row_num
    
    case 1
        d = 2; % problem dimension
        scale_fun = @(t) [t ceil(1.5*t) 2*t ceil(t*log(t)) ceil(0.5*t^2)]; % scalings of m with s
        
        % values for n
        if fig_num == 5
            n_values = 1:15;
        else
            n_values = 2:2:34;
        end
        
    case 2
        d = 4; % problem dimension
        scale_fun = @(t) [t ceil(1.25*t) ceil(1.5*t) 2*t ceil(t*log(t))]; % scalings of m with s
        
        % values for n
        if fig_num == 5
            n_values = 1:11;
        else
            n_values = 3:3:51;
        end
        
    case 3
        d = 8; % problem dimension
        scale_fun = @(t) [t ceil(1.25*t) ceil(1.5*t) 2*t ceil(t*log(t))]; % scalings of m with s
        
        % values for n
        if fig_num == 5
            n_values = 1:5;
        else
            n_values = 1:16;
        end
        
    otherwise
        error('invalid row number selected');
end

%%% Construct error grid %%%

n = max(n_values); I = generate_index_set(index_type,d,n);
s = size(I,2); M = ceil(err_grid_ratio*s);
err_grid = generate_sampling_grid(samp_type,d,M);
b_err_grid = func(err_grid)/sqrt(M);

%%% Main loop %%%

num_cases = length(scale_fun(1)); % number of cases
num_n = length(n_values); % number of n values

% arrays for storing the data
s_values_data = zeros(1,num_n);
cond_num_data = zeros(num_n,num_trials,num_cases);
L2_error_data = zeros(num_n,num_trials,num_cases);
Linf_error_data = zeros(num_n,num_trials,num_cases);

for i = 1:num_n
    
    n = n_values(i);
    I = generate_index_set(index_type,d,n); % generate index set
    A_err_grid = generate_measurement_matrix(poly_type,I,err_grid); % generate error matrix
    
    s = size(I,2); % Define s
    s_values_data(i) = s;
    
    % temporary arrays for parfor
    cond_num_single = zeros(num_trials,num_cases);
    L2_error_single = zeros(num_trials,num_cases);
    Linf_error_single = zeros(num_trials,num_cases);
    
    for l = 1:num_cases
        
        z = scale_fun(s); m = z(l); % select oversampling amount
        
        % temporary array for parfor
        cond_num_temp = zeros(num_trials,1);
        L2_error_temp = zeros(num_trials,1);
        Linf_error_temp = zeros(num_trials,1);
        
        parfor t = 1:num_trials
            
            y_grid = generate_sampling_grid(samp_type,d,m); % generate sample points
            A = generate_measurement_matrix(poly_type,I,y_grid); % generate measurement matrix
            b = func(y_grid)/sqrt(m); % generate measurement vector
            
            kappa = cond(A); % compute condition number
            
            c = A\b; % compute least-squares fit
            L2_err = norm(A_err_grid*c - b_err_grid)/norm(b_err_grid); % compute L^2_rho-norm error
            Linf_err = norm(A_err_grid*c - b_err_grid,Inf)/norm(b_err_grid,Inf); % compute L^2_rho-norm error
            
            cond_num_temp(t) = kappa;
            L2_error_temp(t) = L2_err;
            Linf_error_temp(t) = Linf_err;
            
        end
        
        cond_num_single(:,l) = cond_num_temp;
        L2_error_single(:,l) = L2_error_temp;
        Linf_error_single(:,l) = Linf_error_temp;
        
        disp(['Figure 5.',num2str(fig_num),'_',num2str(row_num),space,'n = ',num2str(n),space,'case = ',num2str(l),space,'s = ',num2str(s),space,'m = ',num2str(m),space,'cond(A) = ',num2str(mean(cond_num_temp)),space,'L2 error = ',num2str(mean(L2_error_temp))]);
        
    end
    
    cond_num_data(i,:,:) = cond_num_single;
    L2_error_data(i,:,:) = L2_error_single;
    Linf_error_data(i,:,:) = Linf_error_single;
    
end

%%% Save data %%%
clear A A_err_grid b_err_grid
save(['../../data/chpt5/',fig_name,'_data.mat'])

end
