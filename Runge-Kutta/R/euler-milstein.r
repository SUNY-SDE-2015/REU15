# Euler/Milstein

r <- 2; K <- 1; beta <- 0.25; Xzero <- 0.5;   % problem parameters
T <- 1; N <- 2^(11); dt <- T/N;              %             
M <- 500;                                  % number of paths sampled
R <- [1; 16; 32; 64; 128];                 % Milstein stepsizes are R*dt

dW <- sqrt(dt)* runif(1, M, N);
Xmil <- double(M, 5);
for p = 1:5
    Dt = R[p]*dt;
    L = N/R[p];
    Xtemp = Xzero*rep(M,1);
    for j = 1:L
        Winc = sum(dW(:, R[p]*(j-1)+1:R[p]*j),2);
        Xtemp = Xzero + Dt*r*Xtemp.*(K-Xtemp) + beta*Xtemp.*Winc + .5*beta^2*Xtemp.*(Winc.^2 - Dt);
    end
    