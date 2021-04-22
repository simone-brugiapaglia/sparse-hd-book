%--- Description ---%
%
% Filename: fig_52_run.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates all data for Figure 5.2

clear all; close all; clc;
addpath ../../utils

for col_num = 1:2
    disp(' ');
    disp(['------------------------------------------------------------------------']);
    disp(['Running Figure 5.2','_',num2str(col_num)]);
    disp(['------------------------------------------------------------------------']);
    disp(' ');
    fig_52_data(col_num);
end
