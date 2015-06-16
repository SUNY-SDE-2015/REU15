
fp=fopen('y_e-y.csv','w');
fprintf(fp,"Time,Approximation,Real,Error\n");

n=0;
xi=1;
y=0;

final=8;
steps=20;
dt=final/steps;
t=xi:dt:final;

clf
hold on
plot(xi,y,'rx')
plot(xi,e^(r*xi),'bx')
title('Another approximation ln y')
xlabel('t')
ylabel('y(t)')

fprintf(fp,"%f,%f,%f,%f\n",xi,y,log(xi),log(xi)-y);

while n<steps
	y=y+exp(-n*dt)*dt+(1-0.5*dt);
	n=n+1;
	s=log(xi+n*dt);
	fprintf(fp,"%f,%f,%f,%f\n",xi+n*dt,y,s,abs(s-y));
	
	hold on
	plot(xi+n*dt,y,'rx');
	plot(xi+n*dt,s,'bx');
	end 
	
print appvsreal00.jpg
fclose(fp);