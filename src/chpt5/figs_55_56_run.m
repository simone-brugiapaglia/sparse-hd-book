%--- Description ---%
%
% Filename: figs_55_56_run.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates all data for Figures 5.5 and 5.6

clear all; close all; clc;
addpath ../../utils

for fig_num = 5:6
    for row_num = 1:3
        disp(' ');
        disp(['------------------------------------------------------------------------']);
        disp(['Running Figure 5.',num2str(fig_num),'_',num2str(row_num)]);
        disp(['------------------------------------------------------------------------']);
        disp(' ');
        figs_55_56_data(fig_num,row_num);
    end
end

