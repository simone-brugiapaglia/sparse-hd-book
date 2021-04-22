%--- Description ---%
%
% Filename: figs_A3_A4_A5_plot.m 
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates the plots for Figures A.3, A.4 and A.5

clear all; close all; clc;
addpath(genpath('../../utils'))

%%% Loop over both figures and all subfigures %%%

for fig_num = 3:5
    for row_num = 1:2
        for col_num = 1:2
            
            fig_name = ['fig_A',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num)];
            load(['../../data/apdxA/',fig_name,'_data'])
            
            %%% L2 error plot %%%
            
            fig = figure;
            
            hPlot = plot_book_style(m_values, L2_error_data, 'shaded', 'mean_std_log10');
            set(gca, 'yscale', 'log')
            axis tight
            
            legend(hPlot, 'wBP','wQCBP (oracle)','wQCBP (CV)','wSR-LASSO (theory)','wSR-LASSO (CV)','location','northeast');
            xlabel('$m$','interpreter', 'latex')
            ylabel('Relative $L^2_{\varrho}(\mathcal{U})$-error','interpreter', 'latex')
            
            set_axis_param
            set_fonts
            
            fig_name = ['fig_A',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num)];
            saveas(fig,['../../figs/apdxA/',fig_name,'_L2'],'epsc');
            
            %%% Linf error plot %%%
            
            fig = figure;
            
            hPlot = plot_book_style(m_values, Linf_error_data, 'shaded', 'mean_std_log10');
            set(gca, 'yscale', 'log')
            axis tight
            
            legend(hPlot, 'wBP','wQCBP (oracle)','wQCBP (CV)','wSR-LASSO (theory)','wSR-LASSO (CV)','location','northeast');
            xlabel('$m$','interpreter', 'latex')
            ylabel('Relative $L^\infty(\mathcal{U})$-error','interpreter', 'latex')
            
            set_axis_param
            set_fonts
            
            saveas(fig,['../../figs/apdxA/',fig_name,'_Linf'],'epsc');
            
        end
    end
end