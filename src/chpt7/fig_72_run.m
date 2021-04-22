%--- Description ---%
%
% Filename: fig_72_run.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates all data for Figure 7.2

clear all; close all; clc;
addpath ../../utils

for row_num = 1:3
    for col_num = 1:2
        disp(' ');
        disp(['------------------------------------------------------------------------']);
        disp(['Running Figure 7.2_',num2str(row_num),'_',num2str(col_num)]);
        disp(['------------------------------------------------------------------------']);
        disp(' ');
        fig_72_data(row_num,col_num);
    end
end
