%--- Description ---%
%
% Filename: fig_11_plot.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates Figure 1.1

close all; clear all; clc;
addpath(genpath('../../utils'))

for row_num = 1:2
    
    fig_name = ['fig_11_',num2str(row_num)];
    
    switch row_num
        case 1
            func = @iso_exp;
            d = 4;
        case 2
            func = @split_product;
            d = 8;
    end
    
    poly_type = 'legendre'; % use Legendre polynomials
    index_type = 'TD'; % use hyperbolic cross index set
    samp_type = 'uniform'; % use samples drawn randomly from the uniform distribution
    
    percent = 10; % percentage largest coefficients to plot
    
    N_max = 200; % desired number of basis functions
    grid_ratio = 500; % ratio between maximum index set size and error grid size
    
    n = find_order(index_type,d,N_max); % find the maximum polynomial order for the given index_type
    I = generate_index_set(index_type,d,n); % compute index set
    N = size(I,2); % actual number of basis functions
    
    %%% Construct grid, measurement matrix and measurement vector %%%
    
    M = ceil(grid_ratio*N);
    y_grid = generate_sampling_grid(samp_type,d,M);
    A = generate_measurement_matrix(poly_type,I,y_grid);
    b = func(y_grid)/sqrt(M);
    
    %%% Compute a reference solution via least squares on the error grid %%%
    c = A\b;
    
    y_range = [min(abs(c))/10,max(abs(c))*10];
    
    %%% Plot results and save %%%
    [ms, lw, fs, colors, markers, AlphaLevel] = get_fig_param();
    
    fig = figure(2*row_num-1);
    
    j = 1;
    semilogy(abs(c),'o','markersize',ms,'MarkerEdgeColor',colors{j},'LineWidth',lw,'Color',colors{j});
    hold on
    
    [c_sort,J] = sort(abs(c),'descend');
    s = ceil(percent*N/100);
    Js = J(1:s);
    
    j = 2;
    semilogy(Js,abs(c(Js)),'o','markersize',ms,'MarkerFaceColor',colors{j},'MarkerEdgeColor',colors{j},'LineWidth',lw,'Color',colors{j});
    hold off
    
    axis tight
    set_axis_param
    set_fonts
    
    saveas(fig,['../../figs/chpt1/',fig_name,'_1'],'epsc');
    
    fig = figure(2*row_num);
    
    j = 1;
    semilogy(c_sort,'o','markersize',ms,'MarkerEdgeColor',colors{j},'LineWidth',lw,'Color',colors{j});
    hold on
    
    j = 2;
    semilogy(1:s,c_sort(1:s),'o','markersize',ms,'MarkerFaceColor',colors{j},'MarkerEdgeColor',colors{j},'LineWidth',lw,'Color',colors{j});
    hold off
    
    axis tight
    set_axis_param
    set_fonts
    
    saveas(fig,['../../figs/chpt1/',fig_name,'_2'],'epsc');
    
end
