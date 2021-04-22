%--- Description ---%
%
% Filename: fig_72_plot.m 
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates the plots for Figure 7.2

clear all; close all; clc;
addpath(genpath('../../utils'))

%%% Loop over all subfigures %%%

for row_num = 1:3
    for col_num = 1:2
        
        fig_name = ['fig_72_',num2str(row_num),'_',num2str(col_num)];
        load(['../../data/chpt7/',fig_name,'_data'])
                
        %%% L2 error plot %%%
       
        fig = figure;
        
        % Create data array
        for j = 1:num_theta
            error_data(:,:,2*j-1) = L2_error_wqcbp_data(:,:,j);
            error_data(:,:,2*j)   = L2_error_wsrlasso_data(:,:,j);
        end
        
        hPlot = plot_book_style(m_values, error_data, 'shaded_pairs', 'mean_std_log10');
        set(gca, 'yscale', 'log')
        axis tight
        
        legend(hPlot,'QCBP','SR-LASSO','wQCBP','wSR-LASSO','location','northeast');
        xlabel('$m$','interpreter', 'latex')
        ylabel('Relative $L^2_\varrho(\mathcal{U})$-error','interpreter', 'latex')
        
        set_axis_param
        set_fonts
        
        saveas(fig,['../../figs/chpt7/',fig_name,'_L2'],'epsc');
        
        %%% Linf error plot %%%
        
        fig = figure;

        % Create data array
        for j = 1:num_theta
            error_data(:,:,2*j-1) = Linf_error_wqcbp_data(:,:,j);
            error_data(:,:,2*j)   = Linf_error_wsrlasso_data(:,:,j);
        end
        
        hPlot = plot_book_style(m_values, error_data, 'shaded_pairs', 'mean_std_log10');
        set(gca, 'yscale', 'log')
        axis tight
        
        legend(hPlot,'QCBP','SR-LASSO','wQCBP','wSR-LASSO','location','northeast');
        xlabel('$m$','interpreter', 'latex')
        ylabel('Relative $L^\infty(\mathcal{U})$-error','interpreter', 'latex')
        
        set_axis_param
        set_fonts
        
        saveas(fig,['../../figs/chpt7/',fig_name,'_Linf'],'epsc');
        
    end
end