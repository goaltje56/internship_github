clear all;
close all;
clc;

iter = 664;
aP = [3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3]';
aW = [0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]';
aE = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0]';
b  = [210 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 1010]';

T = 150*ones(20,1);
n = length(aP);

[T r] = GS_solve(n, T, aW, aE, aP, b, iter);
