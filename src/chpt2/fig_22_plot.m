%--- Description ---%
%
% Filename: fig_22_plot.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates Figure 2.2

clear all; close all; clc;
addpath ../../utils

%%% create a grid for plotting %%%

m = 5000; y_grid = (linspace(-1,1,m))'; 

%%% Plot defaults %%%

[ms, lw, fs, colors, markers, AlphaLevel] = get_fig_param();

%%% Legendre polynomials %%%

fig = figure(1);

% generate matrix whose columns are the Legendre polynomials on the grid
I = 0:6; A = sqrt(m)*generate_measurement_matrix('legendre',I,y_grid);

figure(1);
for i = 1:length(I)
   plot(y_grid,A(:,i),'LineWidth',lw,'Color',colors{i});
   hold on
end

hold off

h = legend('$\Psi_0$','$\Psi_1$','$\Psi_2$','$\Psi_3$','$\Psi_4$','$\Psi_5$','$\Psi_6$','location','southeast');
set(h,'Interpreter','latex');

set_axis_param
set_fonts

saveas(fig,['../../figs/chpt2/','fig_22_1'],'epsc');

%%% Chebyshev polynomials %%%

fig = figure(2);

% generate matrix whose columns are the Chebyshev polynomials on the grid
A = sqrt(m)*generate_measurement_matrix('chebyshev',I,y_grid);

for i = 1:length(I)
   plot(y_grid,A(:,i),'LineWidth',lw,'Color',colors{i});
   hold on
end

hold off

h = legend('$\Psi_0$','$\Psi_1$','$\Psi_2$','$\Psi_3$','$\Psi_4$','$\Psi_5$','$\Psi_6$','location','southeast');
set(h,'Interpreter','latex');

set_axis_param
set_fonts

saveas(fig,['../../figs/chpt2/','fig_22_2'],'epsc');