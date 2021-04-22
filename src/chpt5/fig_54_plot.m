%--- Description ---%
%
% Filename: fig_54_plot.m 
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates the plots for Figure 5.4

clear all; close all; clc;
addpath(genpath('../../utils'))

%%% Loop over all subfigures %%%

for col_num = 1:2
    
    fig_name = ['fig_54_',num2str(col_num)];
    load(['../../data/chpt5/',fig_name,'_data'])
    
    fig = figure(col_num);
    
    hPlot = plot_book_style(s_values, L2_error_data, 'shaded', 'mean_std_log10');
    set(gca, 'yscale', 'log')
    set(gca, 'ylim', [eps, 10])
    
    xlabel('$s$','interpreter', 'latex')
    ylabel('Relative $L^2_\varrho(\mathcal{U})$-error','interpreter', 'latex')
    
    legend(hPlot,'$m = s$','$m = 2 s$','$m = s \log(s)$','$m = 0.5 s^2$','$m = s^2$','$m = s^2 \log(s)$','location','southwest');
    
    set_axis_param
    set_fonts

    saveas(fig,['../../figs/chpt5/',fig_name],'epsc');
    
end
