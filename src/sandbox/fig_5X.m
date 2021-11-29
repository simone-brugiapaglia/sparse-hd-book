%--- Description ---%
%
% Filename: fig_XXXX.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse polynomial approximation of high-dimensional functions", SIAM
%
% Description: test file for the optimal sampling procedure!!


%%%%%% TEMP %%%%%%

clear all; close all; clc;
addpath(genpath('../../utils'))

fig_num = 11;
row_num = 1;
col_num = 2;

%%%%%%%%%%%%%%%%%%

%%% Define main parameters %%%

space = ' ';

poly_type = 'legendre';

err_grid_ratio = 10; % ratio between maximum index set size and error grid size
num_trials = 2; % number of random trials

index_type = 'HC';

scale_fun = @(t) [2*t ceil(t*log(t))]; % scalings of m with s
num_cases = length(scale_fun(1)); % number of cases

fig_name = ['fig_5',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num)];

if fig_num == 11
    func = @iso_exp; % function to approximate
elseif fig_num == 12
    func = @split_product; % function to approximate
else
    error('invalid figure number');
end

if row_num == 1 && col_num == 1
    d = 1;
    n_values = 3:3:75;
    
elseif row_num == 1 && col_num == 2
    d = 2;
    n_values = 13:13:234;
    
elseif row_num == 2 && col_num == 1
    d = 4;
    n_values = 3:3:51;
    
elseif row_num == 2 && col_num == 2
    d = 8;
    n_values = 1:16;
    
else
    error('invalid row or column number selected');
end

%%% Construct error grid %%%

err_samp_type = 'uniform'; % use a uniform grid
n = max(n_values); I = generate_index_set(index_type,d,n);
s = size(I,2); M = ceil(err_grid_ratio*s);
err_grid = generate_sampling_grid(err_samp_type,d,M);
b_err_grid = func(err_grid)/sqrt(M);

%% Main Loop

num_n = length(n_values);
num_polys = 2;

% arrays for storing the data
s_values_data = zeros(1,num_n);
cond_num_data = zeros(num_n,num_trials,num_cases*num_polys);
L2_error_data = zeros(num_n,num_trials,num_cases*num_polys);
Linf_error_data = zeros(num_n,num_trials,num_cases*num_polys);

for p = 1:num_polys
    
    basis_type = 'legendre'; % functions used in the measurement matrix A
    err_basis_type = basis_type; % polynomials used in computing the error
    
    for i = 1:num_n
        n = n_values(i);
        I = generate_index_set(index_type,d,n); % generate index set
        A_err_grid = generate_measurement_matrix(err_basis_type,I,err_grid); % generate error matrix
        [Q,~] = qr(A_err_grid,0); % compute QR decomposition of the error matrix
        
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
            
            %parfor t = 1:num_trials
            
            for t = 1:num_trials
                
                %%% generate sample points, measurement matrix and measurement vector
                
                % uniform
                if p == 1
                    samp_type = 'uniform';
                    y_grid = generate_sampling_grid(samp_type,d,m);
                    A = generate_measurement_matrix(basis_type,I,y_grid);
                    b = func(y_grid)/sqrt(m);
                    
                    
                    % optimal
                elseif p == 2
                    mu = sum(abs(Q.^2)/s,2); % vector of probabilities
                    J = datasample(1:M,m,'Replace',true,'Weights',mu); % sample m indices according to mu
                    w = 1./sqrt(mu);
                    A = w(J).*A_err_grid(J,:)*sqrt(M/m);
                    b = w(J).*b_err_grid(J)*sqrt(M/m);
                    
                end
                
                
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
            
            disp(['Figure 5.',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num),space,'polynomials = ',basis_type,space,'n = ',num2str(n),space,'case = ',num2str(l),space,'s = ',num2str(s),space,'m = ',num2str(m),space,'cond(A) = ',num2str(mean(cond_num_temp)),space,'L2 error = ',num2str(mean(L2_error_temp))]);
            
        end
        
        rng = p:num_polys:num_cases*num_polys;
        
        cond_num_data(i,:,rng) = cond_num_single;
        L2_error_data(i,:,rng) = L2_error_single;
        Linf_error_data(i,:,rng) = Linf_error_single;
        
    end
    
end

%%

perm = [];
for p = 1:num_polys
    perm = [perm, p, num_polys + p];
end

hPlot = plot_book_style(s_values_data, Linf_error_data(:,:,perm), 'shaded_pairs','mean_std_log10');
set(gca, 'yscale', 'log')
axis tight

legend(hPlot,'LU ($m = 2s$)','LU ($m = s \log(s)$)','Opt ($m = 2s$)','Opt ($m = s \log(s)$)','LC ($m = 2s$)','LC ($m = s \log(s)$)','location','northeast');
xlabel('$s$', 'interpreter', 'latex')
ylabel('Relative $L^{\infty}(\mathcal{U})$-error', 'interpreter', 'latex')

set_axis_param
set_fonts


%%% Save data %%%
%clear A A_err_grid b_err_grid
%save(['../../data/chpt5/',fig_name,'_data.mat'])

