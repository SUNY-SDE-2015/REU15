%BPATH2  Brownian path simulation: vectorized 
% Adapted from bpath2.m as given in the paper
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


randn('state',100)          % set the state of randn
T = 1; N = 500; dt = T/N;
tic
dW = sqrt(dt)*randn(1,N);   % increments
W = cumsum(dW);            % cumulative sum
toc
W

plot([0:dt:T],[0,W],'r-')   % plot W against t
xlabel('t','FontSize',16)
ylabel('W(t)','FontSize',16,'Rotation',0)