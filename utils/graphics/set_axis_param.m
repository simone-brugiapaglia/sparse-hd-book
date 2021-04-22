%--- Description ---%
%
% Filename: set_axis_param.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: loads a set of standard axis parameters 

[ms, lw, fs, colors, markers] = get_fig_param();

ax = gca;
ax.YMinorTick = 'on';
ax.YMinorGrid = 'on';
ax.FontSize = fs;
ax.LineWidth = lw;
ax.YGrid = 'on';
ax.XGrid = 'on';
ax.XMinorTick = 'on';
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.YMinorTick = 'on';