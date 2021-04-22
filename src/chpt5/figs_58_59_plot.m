%--- Description ---%
%
% Filename: fig_58_59_plot.m 
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates the plots for Figures 5.8 and 5.9

clear all; close all; clc;
addpath(genpath('../../utils'))

%%% Loop over both figures and all subfigures %%%

for fig_num = 8:9
    for row_num = 1:2
        for col_num = 1:2
            
            fig_name = ['fig_5',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num)];
            load(['../../data/chpt5/',fig_name,'_data'])
            
            fig = figure(2*(row_num-1)+col_num);
            
            % Create permutation to plot data in the correct order
            perm = [];
            for p = 1:num_polys
                perm = [perm, p, num_polys + p];
            end
            
            % Plot data
            hPlot = plot_book_style(s_values_data, Linf_error_data(:,:,perm), 'shaded_pairs','mean_std_log10');
            set(gca, 'yscale', 'log')
            axis tight
            
            legend(hPlot,'LU ($m = 2s$)','LU ($m = s \log(s)$)','CC ($m = 2s$)','CC ($m = s \log(s)$)','LC ($m = 2s$)','LC ($m = s \log(s)$)','location','northeast');
            xlabel('$s$', 'interpreter', 'latex')
            ylabel('Relative $L^{\infty}(\mathcal{U})$-error', 'interpreter', 'latex')
                        
            set_axis_param
            set_fonts
        
            saveas(fig,['../../figs/chpt5/',fig_name],'epsc');
            
        end
    end
end
