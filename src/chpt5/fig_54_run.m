%--- Description ---%
%
% Filename: fig_54_run.m
% Authors: Ben Adcock, Simone Brugiapaglia and Clayton Webster
% Part of the book "Sparse Polynomial Approximation of High-Dimensional
% Functions", SIAM, 2021
%
% Description: generates all data for Figure 5.4

clear all; close all; clc;
addpath ../../utils

for col_num = 1:2
    disp(' ');
    disp(['------------------------------------------------------------------------']);
    disp(['Running Figure 5.4','_',num2str(col_num)]);
    disp(['------------------------------------------------------------------------']);
    disp(' ');
    fig_54_data(col_num);
end
