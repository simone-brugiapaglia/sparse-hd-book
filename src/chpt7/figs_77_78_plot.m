%--- Description ---%
%
% Filename: figs_77_78_plot.m 
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates the plots for Figures 7.7 and 7.8

clear all; close all; clc;
addpath(genpath('../../utils'))

%%% Loop over all figures and subfigures %%%

for fig_num = 7:8
    for row_num = 1:3
        for col_num = 1:2
            
            fig_name = ['fig_7',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num)];
            load(['../../data/chpt7/',fig_name,'_data'])
            
            fig = figure;
            
            hPlot = plot_book_style(m_values, Linf_error_data, 'shaded', 'mean_std_log10');
            set(gca, 'yscale', 'log')
            axis tight
            
            legend(hPlot, 'LU','CC','LC','location','northeast');
            xlabel('$m$','interpreter', 'latex')
            ylabel('Relative $L^\infty_{\varrho}(\mathcal{U})$-error','interpreter', 'latex')
            
            set_axis_param
            set_fonts
            
            saveas(fig,['../../figs/chpt7/',fig_name],'epsc');
            
        end
    end
end
