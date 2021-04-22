%--- Description ---%
%
% Filename: figs_61_62_run.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates all data for Figures 6.1 and 6.2

clear all; close all; clc;
addpath ../../utils

for fig_num = 1:2
    for row_num = 1:3
        for col_num = 1:2
            disp(' ');
            disp(['------------------------------------------------------------------------']);
            disp(['Running Figure 6.',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num)]);
            disp(['------------------------------------------------------------------------']);
            disp(' ');
            figs_61_62_data(fig_num,row_num,col_num);
        end
    end
end