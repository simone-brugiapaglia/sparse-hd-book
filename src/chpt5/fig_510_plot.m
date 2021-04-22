%--- Description ---%
%
% Filename: fig_510_plot.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates the plots for Figure 5.10

clear all; close all; clc;
addpath(genpath('../../utils'))

%%% Create a grid for plotting %%%

m = 1000; y_grid = (linspace(-1,1,m))';

%%% Plotting defaults %%%

[ms, lw, fs, colors, markers, AlphaLevel] = get_fig_param();

s_values = [4 8 16];
s_max = max(s_values);
num_s = length(s_values);

%%% Chebyshev density %%%

cheb_density = @(y) 1./(pi*sqrt(1-y.^2));

% generate matrix whose columns are the Legendre polynomials on the grid
I = 0:(s_max-1); A = sqrt(m)*generate_measurement_matrix('legendre',I,y_grid);

for l = 1:num_s
    
    s = s_values(l);
    K = sum(abs(A(:,1:s)).^2,2)/(2*s); % compute K(y)/(2*s)
    
    fig = figure(l);
    plot(y_grid,K,'LineWidth',lw,'Color',colors{1});
    hold on
    plot(y_grid,cheb_density(y_grid),'--','LineWidth',lw,'Color',[0,0,0]);
    hold off
    
    ylim([0.2,2]);
    
    set_axis_param
    set_fonts
    
    saveas(fig,['../../figs/chpt5/','fig_510_',num2str(l)],'epsc');
    
end






