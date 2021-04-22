%--- Description ---%
%
% Filename: fig_54_data.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse polynomial approximation of high-dimensional functions", SIAM
%
% Description: generates the data for Figure 5.4
% 
% Input: 
% col_num - column number (either 1 or 2)

function fig_54_data(col_num)

%%% Define main parameters %%%

space = ' ';

d = 1; % problem dimension
poly_type = 'legendre'; % use Legendre polynomials
samp_type = 'uniform'; % use samples drawn randomly from the uniform distribution

s_values = 2:2:50; % range of problem sizes
err_grid_ratio = 10; % ratio between maximum index set size and error grid size
num_trials = 50; % number of random trials

scale_fun = @(t) [t 2*t ceil(t*log(t)) ceil(0.5*t^2) t^2 ceil(t^2*log(t))];

%%% Select which subfigure to generate %%%

switch col_num
    case 1
        func = @(y) exp(3*y);
    case 2
        func = @(y) 1./(1+y.^2);
    otherwise
        error('invalid figure number selected');
end

fig_name = ['fig_54_',num2str(col_num)];

%%% Construct error grid %%%

s = max(s_values); I = 0:(s-1); % find maximum index set size
M = ceil(err_grid_ratio*s);
err_grid = generate_sampling_grid(samp_type,d,M);
A_err_grid = generate_measurement_matrix(poly_type,I,err_grid);
b_err_grid = func(err_grid)/sqrt(M);

%%% Main loop %%%

num_cases = length(scale_fun(1)); % number of cases
num_s = length(s_values); % number of s values

% arrays for storing the data
L2_error_data = zeros(num_s,num_trials,num_cases);
Linf_error_data = zeros(num_s,num_trials,num_cases);

for i = 1:num_s
   
    s = s_values(i); I = 0:(s-1); % generate index set
    
    % temporary arrays for parfor
    L2_error_single = zeros(num_trials,num_cases);
    Linf_error_single = zeros(num_trials,num_cases);
    
    for l = 1:num_cases
        
       z = scale_fun(s); m = z(l); % select oversampling amount
        
        % temporary array for parfor
        L2_error_temp = zeros(num_trials,1); 
        Linf_error_temp = zeros(num_trials,1); 
        
        parfor t = 1:num_trials
            
            y_grid = generate_sampling_grid(samp_type,d,m); % generate sample points
            A = generate_measurement_matrix(poly_type,I,y_grid); % generate measurement matrix
            b = func(y_grid)/sqrt(m); % generate measurement vector
            
            c = A\b; % compute least-squares fit
            L2_err = norm(A_err_grid(:,1:s)*c - b_err_grid)/norm(b_err_grid); % compute L^2_rho-norm error
            Linf_err = norm(A_err_grid(:,1:s)*c - b_err_grid,Inf)/norm(b_err_grid,Inf); % compute L^2_rho-norm error
            
            L2_error_temp(t) = L2_err;
            Linf_error_temp(t) = Linf_err;
            
        end
        
        L2_error_single(:,l) = L2_error_temp;
        Linf_error_single(:,l) = Linf_error_temp;
       
        disp(['Figure 5.4_',num2str(col_num),space,'s = ',num2str(s),space,'case = ',num2str(l),space,'s = ',num2str(s),space,'m = ',num2str(m),space,'L2 error = ',num2str(mean(L2_error_temp))]);
        
    end
    
    L2_error_data(i,:,:) = L2_error_single;
    Linf_error_data(i,:,:) = Linf_error_single;
    
end

%%% Save data %%%
clear A A_err_grid b_err_grid
save(['../../data/chpt5/',fig_name,'_data.mat'])

end