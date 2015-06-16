
fp=fopen('y_e-y.csv','w');
fprintf(fp,"Time,Approximation,Real,Error\n");

n=0;
xi=0;
y=1;

final=8;
steps=20;
dt=final/steps;
t=xi:dt:final;

clf
hold on
plot(xi,y,'rx')
plot(xi,e^(r*xi),'bx')
title('Another approximation e^-y')
xlabel('t')
ylabel('y(t)')

fprintf(fp,"%f,%f,%f,%f\n",xi,y,ln(xi),ln(xi)-y);

while n<steps
	y=y+exp(-n*dt)*dt+(1-0.5*dt);
	n=n+1;
	s=ln(xi+n*dt);
	fprintf(fp,"%f,%f,%f,%f\n",xi+n*dt,y,s,abs(s-y));
	
	hold on
	plot(xi+n*dt,y,'rx');
	plot(xi+n*dt,s,'bx');
	end 
	
print appvsreal00.jpg
fclose(fp);