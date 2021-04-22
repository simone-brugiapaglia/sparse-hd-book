%--- Description ---%
%
% Filename: fig_53_data.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse polynomial approximation of high-dimensional functions", SIAM
%
% Description: generates the data for Figure 5.3
% 
% Input: 
% col_num - column number (either 1 or 2)

function fig_53_data(col_num)

%%% Define main parameters %%%

space = ' ';

d = 1; % problem dimension
poly_type = 'legendre'; % use Legendre polynomials
samp_type = 'uniform'; % use samples drawn randomly from the uniform distribution

s_values = 2:2:50; % range of problem sizes
err_grid_ratio = 10; % ratio between maximum index set size and error grid size
num_trials = 50; % number of random trials

fig_name = ['fig_53_',num2str(col_num)]; % filenames

%% Main loop %%

switch col_num
    case 1
        scale_fun = @(t) [t ceil(1.5*t) 2*t 4*t ceil(t*log(t)) ceil(2*t*log(t))]; % scalings of m with s
    case 2
        scale_fun = @(t) [ceil(0.5*t^2) t^2 ceil(1.5*t^2) 2*t^2 ceil(t^2*log(t)) ceil(2*t^2*log(t))]; % scalings of m with s
    otherwise
        error('invalid column number');
end

num_cases = length(scale_fun(1));
num_s = length(s_values);

cond_num_data = zeros(num_s,num_trials,num_cases);

for i = 1:num_s
    
    s = s_values(i); I = 0:(s-1); % generate index set
    
    cond_num_single = zeros(num_trials,num_cases); % temporary array for parfor
    
    for l = 1:num_cases
        
        z = scale_fun(s); m = z(l); % select oversampling amount
        
        cond_num_temp = zeros(num_trials,1); % temporary array for parfor
        
        parfor t = 1:num_trials
            
            % generate sampling grid and matrix
            y_grid = generate_sampling_grid(samp_type,d,m);
            A = generate_measurement_matrix(poly_type,I,y_grid);
            
            % compute condition number and store value
            kappa = cond(A);
            cond_num_temp(t) = kappa;
            
        end
        
        cond_num_single(:,l) = cond_num_temp;
        
        disp(['Figure 5.3_',num2str(col_num),space,'s = ',num2str(s),space,'case = ',num2str(l),space,'s = ',num2str(s),space,'m = ',num2str(m),space,'cond(A) = ',num2str(mean(cond_num_temp))]);
    end
    
    cond_num_data(i,:,:) = cond_num_single;
    
end

%%% Save data %%%
clear A
save(['../../data/chpt5/',fig_name,'_data.mat'])

end
