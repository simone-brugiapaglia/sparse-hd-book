%--- Description ---%
%
% Filename: figs_61_62_plot.m 
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates the plots for Figures 6.2 and 6.3

clear all; close all; clc;
addpath(genpath('../../utils'))

%%% Loop over both figures and all subfigures %%%

for fig_num = 2:3
    for row_num = 1:3
        for col_num = 1:2
            
            fig_name = ['fig_6',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num)];
            load(['../../data/chpt6/',fig_name,'_data'])
            
            fig = figure;
            
            hPlot = plot_book_style(param_values, error_data, 'shaded', 'mean_std_log10');
            set(gca, 'yscale', 'log', 'xscale', 'log')
            axis tight
            
            xlabel('Tuning parameter','interpreter', 'latex')
            ylabel('Relative $L^2_\varrho(\mathcal{U})$-error','interpreter', 'latex')
            
            legend(hPlot, '$\sigma = 10^{-1}$','$\sigma = 10^{-2}$','$\sigma = 10^{-3}$','$\sigma = 0$','location','southeast');
            
            set_axis_param
            set_fonts

            saveas(fig,['../../figs/chpt6/',fig_name],'epsc');
            
        end
    end
end
