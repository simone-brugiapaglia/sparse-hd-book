%--- Description ---%
%
% Filename: set_fonts.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: loads a set of graphical parameters related to fonts

[ms, lw, fs, colors, markers] = get_fig_param();

set(gca, 'TickLabelInterpreter', 'latex')
set(gca, 'FontSize', fs)

hLegend = findobj(gcf, 'Type', 'Legend');
set(hLegend, 'interpreter', 'latex', 'fontsize',fs)

hXLabel = get(findobj(gcf,'type','axe'), 'xlabel');
set(hXLabel, 'interpreter', 'latex', 'fontsize',fs)

hYLabel = get(findobj(gcf,'type','axe'), 'ylabel');
set(hYLabel, 'interpreter', 'latex', 'fontsize',fs)
