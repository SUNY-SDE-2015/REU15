# Euler/Milstein

r <- 2; K <- 1; beta <- 0.25; Xzero <- 0.5;   # problem parameters
T <- 1; N <- 2^(11); dt <- T/N;               #             
M <- 500;                                     # number of paths sampled
R <- [1; 16; 32; 64; 128];                    # Milstein stepsizes are R*dt

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
    Xmil(:,p) <- Xtemp; #store MIlstein solution at t=1
end

Xref = Xmil(:,1);                            # Reference solution
Xerr = abs(Xmil(:,2:5) - repmat (Xref,1,4)); # Error in each path
mean(Xerr);                                  # Mean pathwise errors
Dtvals = dt*R(2:5);                          # Milstein timesteps used

subplot(224)
plot(Dtvals, mean(Xerr), log = "xy", type = "l");
lines(Dtvals, Dtvals, log = "xy");
axis([1e-3 1e-1 1e-4 1])
xlabel('\Delta t')
ylabel('Sample average of | X(T) - X_L|')
main('milstein', 'FontSize', 10)