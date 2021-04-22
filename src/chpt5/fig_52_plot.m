%--- Description ---%
%
% Filename: fig_52_plot.m 
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates the plots for Figure 5.2

clear all; close all; clc;
addpath(genpath('../../utils'))

%%% Plotting defaults %%%

[ms, lw, fs, colors, markers, AlphaLevel] = get_fig_param();

%%% Loop over all subfigures %%%

for col_num = 1:2
    
    fig_name = ['fig_52_',num2str(col_num)];
    load(['../../data/chpt5/',fig_name,'_data'])
    
    fig = figure(col_num);
    
    s_min = min(s_values_data(:)); s_max = max(s_values_data(:));

    for l = 1:num_d
        loglog(s_values_data(:,l),kappa_values_data(:,l),markers{l},'markersize',ms,'MarkerFaceColor',colors{l},'MarkerEdgeColor',colors{l},'LineWidth',lw,'Color',colors{l});
        hold on
    end
    loglog([s_min s_max],[s_min^2 s_max^2],'k--','LineWidth',lw);
    loglog([s_min s_max],[s_min s_max],'k--','LineWidth',lw);
    
    hold off  
    
    xlabel('$|S|$', 'interpreter','latex')
    ylabel('$\kappa(\mathcal{P}_S)$', 'interpreter','latex')

    h = legend('$d = 2$','$d = 4$','$d = 8$','$d = 16$','location','northwest');
    set(h,'Interpreter','latex');
    
    axis tight
    
    set_axis_param
    set_fonts
    
    saveas(fig,['../../figs/chpt5/',fig_name],'epsc');
    
end
