%RungeKutta Approximation and Heun's Method
r = -5;
a = 1/2;
b = 1/2;
alpha = 1;
beta = 1;


final = 100;
steps = 10000;
dt = final/steps;


n = 0;
y = 1;
clf
semilogx(0,y,'rx')
%loglog(0,y,'rx')
%semilogy(0,y,'rx')
while n<steps
    k1 = dt*r*y;
    k2 = dt*r*(y+beta*k1);
    y = y+a*k1+b*k2;
    hold on;
    %plot(dt*(n+1),y-exp(r*n),'rx')
    %loglog(dt*(n+1),y-exp(r*n),'bx')
    semilogx(dt*(n+1),y-exp(r*n),'gx')
    %semilogy(dt*(n+1),y-exp(r*n),'kx')
    n=n+1;
end
xlabel('Error');
ylabel('Time');
title('SemiLogX with r = -5 and dt = 10^-2');
grid on