%--- Description ---%
%
% Filename: fig_53_plot.m 
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates the plots for Figure 5.3

clear all; close all; clc;
addpath(genpath('../../utils'))

%%% Loop over all subfigures %%%

for col_num = 1:2
    
    fig_name = ['fig_53_',num2str(col_num)];
    load(['../../data/chpt5/',fig_name,'_data'])
    
    fig = figure(col_num);
    
    hPlot = plot_book_style(s_values, cond_num_data, 'shaded', 'mean_std_log10');
    set(gca, 'yscale', 'log')  
    axis tight
    
    switch col_num
        case 1
            legend(hPlot,'$m = s$','$m =  1.5 s $','$m = 2 s$','$m = 4 s$','$m =  s \log(s) $','$m =  2 s \log(s)  $','location','northwest');
        case 2
            legend(hPlot,'$ m =  0.5 s^2 $','$m = s^2$','$m =  1.5 s^2 $','$m = 2 s^2$','$m =  s^2 \log(s) $','$m =  2 s^2 \log(s) $','location','northwest');
    end
    
    xlabel('$s$','interpreter', 'latex')
    ylabel('$\mathrm{cond}(${\boldmath$A$}$)$','interpreter', 'latex')
    
    if col_num == 2
        ylim([1,1e2])
    end
    
    set_axis_param
    set_fonts
    
    saveas(fig,['../../figs/chpt5/',fig_name],'epsc');
    
end
