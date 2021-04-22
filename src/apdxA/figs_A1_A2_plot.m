%--- Description ---%
%
% Filename: figs_A1_A2_plot.m 
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates the plots for Figures A.1 and A.2

clear all; close all; clc;
addpath(genpath('../../utils'))

%%% Plotting defaults %%%

[ms, lw, fs, colors, markers, AlphaLevel] = get_fig_param();

%%% Loop over all subfigures %%%

i_curve = 3;

for fig_number = 1:2  
    for col_number = 1:2
        
        switch fig_number
            case 1
                plot_type = 'mean_std_log10';
                plot_name = 'shaded plot (this book''s style)';
                figure_name = ['fig_A1_',num2str(col_number)]; % filenames
                
            case 2
                plot_type = 'mean_std_eps';
                plot_name = 'shaded plot (without $\log_{10}$ transformation)';
                figure_name = ['fig_A2_',num2str(col_number)]; % filenames       
        end
             
        switch col_number
            case 1
                
                load('../../data/chpt5/fig_53_1_data')
                
                fig = figure;
                
                hPlot = plot_book_style(s_values, cond_num_data(:,:,i_curve), 'shaded', plot_type);
                hold on
                hPlotRaw = plot(s_values, cond_num_data(:,:,i_curve),'*k');
                
                set(gca, 'yscale', 'log')
                axis tight
                
                ylim([1, 10^13])  
                
                legend([hPlot, hPlotRaw(1)], plot_name, 'raw data', 'location','northwest')
                xlabel('$s$','interpreter', 'latex')
                ylabel('$\mathrm{cond}(${\boldmath$A$}$)$','interpreter', 'latex')
                
    
            case 2
                
                load('../../data/chpt5/fig_54_1_data')
                
                fig = figure;
                
                hPlot = plot_book_style(s_values, L2_error_data(:,:,i_curve), 'shaded', plot_type);
                hold on
                hPlotRaw = plot(s_values, L2_error_data(:,:,i_curve),'*k', 'markersize',ms);
                
                set(gca, 'yscale', 'log')
                
                axis tight
                
                ylim([eps, 10])
                
                legend([hPlot, hPlotRaw(1)], plot_name, 'raw data', 'location','northeast')
                xlabel('$s$','interpreter', 'latex')
                ylabel('Relative $L^2_\varrho(\mathcal{U})$-error','interpreter', 'latex')
                
        end
        
        set_axis_param
        set_fonts

        
        saveas(fig,['../../figs/apdxA/',figure_name],'epsc');
    end
end
