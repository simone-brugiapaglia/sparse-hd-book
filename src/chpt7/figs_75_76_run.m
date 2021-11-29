%--- Description ---%
%
% Filename: figs_75_76_run.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates all data for Figures 7.5 and 7.6

clear all; close all; clc;
addpath ../../utils

for fig_num = 5:6
    for row_num = 1:2
        for col_num = 1:2
             
            disp(' ');
            disp(['------------------------------------------------------------------------']);
            disp(['Running Figure 7.',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num)]);
            disp(['------------------------------------------------------------------------']);
            disp(' ');
            figs_75_76_data(fig_num,row_num,col_num);
            
        end
    end
end