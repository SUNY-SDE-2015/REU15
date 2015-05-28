%EM  Euler-Maruyama vs. Milstein
%
% SDE is  dX = lambda*X dt + mu*X dW,   X(0) = Xzero,
%      where lambda = 2, mu = 1 and Xzero = 1.
%
% Discretized Brownian path over [0,1] has dt = 2^(-8).
%% EM
% Euler-Maruyama uses timestep R*dt.

s=5000; %
emerror=zeros(s,2); %
%approx=zeros(s,1);

%randn('state',100)
lambda = 2; mu = 1; Xzero = 1;    % problem parameters
T = 1; N = 2^9; dt = T/N;    
for i=1:s %
    dW = sqrt(dt)*randn(1,N);         % Brownian increments
    W = cumsum(dW);                   % discretized Brownian path 

    Xtrue = Xzero*exp((lambda-0.5*mu^2)*([dt:dt:T])+mu*W); 
    plot([0:dt:T],[Xzero,Xtrue],'m-'), hold on

    R = 4; Dt = R*dt; L = N/R;        % L EM steps of size Dt = R*dt
    Xem = zeros(1,L);                 % preallocate for efficiency
    Xtemp = Xzero;
    for j = 1:L
       Winc = sum(dW(R*(j-1)+1:R*j)); 
       Xtemp = Xtemp + Dt*lambda*Xtemp + mu*Xtemp*Winc;
       Xem(j) = Xtemp;
    end

    plot([0:Dt:T],[Xzero,Xem],'r--*'), hold off
    xlabel('t','FontSize',12)
    ylabel('X','FontSize',16,'Rotation',0,'HorizontalAlignment','right')

    emerr = abs(Xem(end)-Xtrue(end));

    emerror(i,1)=emerr; %
    emerror(i,2)=Xem(end); %
    i=i+1; %

end %
%% Milstein

s=5000; %
milerror=zeros(s,2); %
%approx=zeros(s,1);

%randn('state',100)
lambda = 2; mu = 1; Xzero = 1;    % problem parameters
T = 1; N = 2^9; dt = T/N;    
for i=1:s %
    dW = sqrt(dt)*randn(1,N);         % Brownian increments
    W = cumsum(dW);                   % discretized Brownian path 

    Xtrue = Xzero*exp((lambda-0.5*mu^2)*([dt:dt:T])+mu*W); 
    plot([0:dt:T],[Xzero,Xtrue],'m-'), hold on

    R = 4; Dt = R*dt; L = N/R;        % L EM steps of size Dt = R*dt
    X = zeros(1,L);                 % preallocate for efficiency
    Xtemp = Xzero;
    for j = 1:L
       Winc = sum(dW(R*(j-1)+1:R*j)); 
       Xtemp = Xtemp + Dt*lambda*Xtemp + mu*Xtemp*Winc...
           +.5*mu*mu*Xtemp*((Winc*Winc) - Dt);
       Xmil(j) = Xtemp;
    end

    plot([0:Dt:T],[Xzero,Xmil],'r--*'), hold off
    xlabel('t','FontSize',12)
    ylabel('X','FontSize',16,'Rotation',0,'HorizontalAlignment','right')

    milerr = abs(Xmil(end)-Xtrue(end));

    milerror(i,1)=milerr; %
    milerror(i,2)=Xmil(end); %
    i=i+1; %

end %