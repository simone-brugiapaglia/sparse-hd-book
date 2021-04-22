%--- Description ---%
%
% Filename: fig_57_plot.m 
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates the plots for Figure 5.7

clear all; close all; clc;
addpath(genpath('../../utils'))

%%% Create a grid for plotting %%%

m = 5000; y_grid = (linspace(-1,1,m))'; 

%%% Plotting defaults %%%

[ms, lw, fs, colors, markers, AlphaLevel] = get_fig_param();

%%% Legendre polynomials and envelope bound %%%

env = @(y) 2./(sqrt(pi)*(1-y.^2).^(1/4));

% generate matrix whose columns are the Legendre polynomials on the grid
I = 0:9; A = sqrt(m)*generate_measurement_matrix('legendre',I,y_grid);

fig = figure(1);

for i = 1:length(I)
   plot(y_grid,A(:,i),'LineWidth',lw,'Color',colors{1});
   hold on
end

plot(y_grid,env(y_grid),'LineWidth',lw*1.5,'Color',colors{2});
plot(y_grid,-env(y_grid),'LineWidth',lw*1.5,'Color',colors{2});

hold off

set_axis_param
set_fonts

saveas(fig,['../../figs/chpt5/','fig_57_1'],'epsc');

%%% Preconditioned Legendre polynomials %%%

% generate matrix whose columns are the Legendre polynomials on the grid
I = 0:9; A = sqrt(m)*generate_measurement_matrix('preconditioned',I,y_grid);

fig = figure(2);

for i = 1:length(I)
   plot(y_grid,A(:,i),'LineWidth',lw,'Color',colors{1});
   hold on
end

plot(y_grid,2*ones(size(y_grid)),'LineWidth',lw*1.5,'Color',colors{2});
plot(y_grid,-2*ones(size(y_grid)),'LineWidth',lw*1.5,'Color',colors{2});
plot(y_grid,sqrt(pi)*ones(size(y_grid)),'--','LineWidth',lw*1.5,'Color',colors{2});

hold off

set_axis_param
set_fonts

ylim([-2.5,2.5])

saveas(fig,['../../figs/chpt5/','fig_57_2'],'epsc');
