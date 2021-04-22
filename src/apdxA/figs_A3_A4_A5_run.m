%--- Description ---%
%
% Filename: figs_A3_A4_A5_run.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates all data for Figures A.3, A.4 and A.5

clear all; close all; clc;
addpath ../../utils

for fig_num = 3:5
    for row_num = 1:2
        for col_num = 1:2
            disp(' ');
            disp(['------------------------------------------------------------------------']);
            disp(['Running Figure A.',num2str(fig_num),'_',num2str(row_num),'_',num2str(col_num)]);
            disp(['------------------------------------------------------------------------']);
            disp(' ');
            figs_A3_A4_A5_data(fig_num,row_num,col_num);
        end
    end
end
