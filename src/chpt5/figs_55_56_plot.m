%--- Description ---%
%
% Filename: figs_55_56_plot.m 
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates the plots for Figures 5.5 and 5.6

clear all; close all; clc;
addpath(genpath('../../utils'))

%%% Loop over both figures and all subfigures %%%

for fig_num = 5:6
    for row_num = 1:3
        
        fig_name = ['fig_5',num2str(fig_num),'_',num2str(row_num)];
        load(['../../data/chpt5/',fig_name,'_data'])
        
        %%% plot cond(A) versus s %%%
        
        fig = figure(2*row_num-1);
        
        hPlot = plot_book_style(s_values_data, cond_num_data, 'shaded', 'mean_std_log10');
        set(gca, 'yscale', 'log')
        axis tight
        
        if row_num == 1
            legend(hPlot, '$m = s$','$m = 1.5 s$','$m = 2 s$','$m = s \log(s)$','$m = 0.5 s^2$','location','northwest');
        elseif row_num == 2
            legend(hPlot, '$m = s$','$m = 1.25 s$','$m = 1.5 s$','$m = 2 s$','$m = s \log(s)$','location','northwest');
        else
            legend(hPlot, '$m = s$','$m = 1.25 s$','$m = 1.5 s$','$m = 2 s$','$m = s \log(s)$','location','northwest');
        end
        xlabel('$s$','interpreter', 'latex')
        ylabel('$\mathrm{cond}(${\boldmath$A$}$)$','interpreter', 'latex')
        
        set_axis_param
        set_fonts
        
        saveas(fig,['../../figs/chpt5/',fig_name,'_1'],'epsc');
        
        %%% plot relative error versus s %%%
        
        fig = figure(2*row_num);
        
        hPlot = plot_book_style(s_values_data, L2_error_data, 'shaded', 'mean_std_log10');
        set(gca, 'yscale', 'log')
        axis tight
        
        if row_num == 1
            h = legend(hPlot, '$m = s$','$m = 1.5 s$','$m = 2 s$','$m = s \log(s)$','$m = 0.5 s^2$','location','northeast');
        elseif row_num == 2
            h = legend(hPlot, '$m = s$','$m = 1.25 s$','$m = 1.5 s$','$m = 2 s$','$m = s \log(s)$','location','northeast');
        else
            h = legend(hPlot, '$m = s$','$m = 1.25 s$','$m = 1.5 s$','$m = 2 s$','$m = s \log(s)$','location','northeast');
        end
        xlabel('$s$','interpreter', 'latex')
        ylabel('Relative $L^2_\varrho(\mathcal{U})$-error','interpreter', 'latex')
        
        set_axis_param
        set_fonts
        
        saveas(fig,['../../figs/chpt5/',fig_name,'_2'],'epsc');
        
    end
end
