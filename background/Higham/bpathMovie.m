%BPATH Brownian path simulation: 
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

% Script to generate the path for a Brownian motion and plot it.

%randn('state',100)          % set the state of randn
T=1; N=250; dt=T/N;          % Set up the bounds
skip = 5;                    % Plot every "skip" time steps.
colours = ['r-';'b-';'g-'];   % Define the colours to plot in
timeValues = 0:dt:T;         % Set up the time values

% Set up the plotting window so it looks all nice and proper like.
clf;
axis([0 T -4 4])
xlabel('t','FontSize',16)
ylabel('W(t)','FontSize',16,'Rotation',0)
title('Random Walks Examples')
hold on
clear F;

frameNumber = 1;
for runNumber = 1:size(colours,1),
    
    currentRun = colours(runNumber,:);
    %printf('Printing using colour %s\n',currentRun)
    
    % Generate a sample Brownian motion.
    dW=sqrt(dt)*randn(1,N);     % increments
    W=cumsum(dW);               % cumulative sum

    % Plot it in time.
    for k=2:skip:size(W,2),
        plot(timeValues(1:k+1),[0,W(1:k)],currentRun)   % plot W against t
        drawnow
        F(frameNumber) = getframe(gcf);
        frameNumber = frameNumber + 1;
    end
    
end

%movie(F)

%writerObj = VideoWriter('randomWalk.avi','Uncompressed AVI');
writerObj = VideoWriter('randomWalk.mp4','MPEG-4');
open(writerObj);
writeVideo(writerObj,F)
