%--- Description ---%
%
% Filename: figs_73_74_plot.m 
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates the plots for Figures 7.3 and 7.4

clear all; close all; clc;
addpath(genpath('../../utils'))

%%% Loop over all figures and subfigures %%%

for fig_num = 3:4
    for row_num = 1:3
        for col_num = 1:2
            
            fig_name = ['fig_7',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num)];
            disp(fig_name)
            load(['../../data/chpt7/',fig_name,'_data'])
            disp(fig_name)
            
            %%% L2 error plot %%%
            
            fig = figure;
            
            hPlot = plot_book_style(m_values, L2_error_data, 'shaded', 'mean_std_log10');
            set(gca, 'yscale', 'log')
            axis tight

            
            legend(hPlot, '$\theta = 0.0$','$\theta = 0.5$','$\theta = 1.0$','$\theta = 1.5$','$\theta = 2.0$','$\theta = 3.0$','$\theta = 4.0$','location','northeast');
            xlabel('$m$','interpreter', 'latex')
            ylabel('Relative $L^2_{\varrho}(\mathcal{U})$-error','interpreter', 'latex')
            
            set_axis_param
            set_fonts

            saveas(fig,['../../figs/chpt7/',fig_name,'_L2'],'epsc');
            
            %%% Linf error plot %%%
            
            fig = figure;
            
            hPlot = plot_book_style(m_values, Linf_error_data, 'shaded', 'mean_std_log10');
            set(gca, 'yscale', 'log')
            axis tight
            
            legend(hPlot, '$\theta = 0.0$','$\theta = 0.5$','$\theta = 1.0$','$\theta = 1.5$','$\theta = 2.0$','$\theta = 3.0$','$\theta = 4.0$','location','northeast');
            xlabel('$m$','interpreter', 'latex')
            ylabel('Relative $L^\infty(\mathcal{U})$-error','interpreter', 'latex')
            
            set_axis_param
            set_fonts
            
            saveas(fig,['../../figs/chpt7/',fig_name,'_Linf'],'epsc');
            
        end
    end
end