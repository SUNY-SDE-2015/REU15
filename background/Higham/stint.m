%STINT  Approximate stochastic integrals
%
% Ito and Stratonovich integrals of W dW
%
% Adapted from stint.m as given in the paper
%
% Higham., D.,
% An Algorithmic Introduction to Numerical Simulation of Stochastic Differential Equations
% SIAM Review
% Volume 43
% Number 3
% pages 525-546
% 2001
% doi 10.1137/S0036144500378302
% http://dx.doi.org/10.1137/S0036144500378302
%


randn('state',100)                      % set the state of randn
T = 1; N = 500; dt = T/N;

dW = sqrt(dt)*randn(1,N);               % increments
W = cumsum(dW);                         % cumulative sum

ito = sum([0,W(1:end-1)].*dW)    
strat = sum((0.5*([0,W(1:end-1)]+W) + 0.5*sqrt(dt)*randn(1,N)).*dW) 

itoerr = abs(ito - 0.5*(W(end)^2-T))
straterr = abs(strat - 0.5*W(end)^2)
