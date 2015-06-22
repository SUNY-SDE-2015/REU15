
fp=fopen('y_y00.csv','w');
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
title('Another approximation')
xlabel('t')
ylabel('y(t)')

fprintf(fp,"%f,%f,%f,%f\n",xi,y,e^(xi),e(xi)-y);

while n<steps
	y=y*(1+dt+0.5*dt^2);
	n=n+1;
	s=e^(xi+n*dt);
	fprintf(fp,"%f,%f,%f,%f\n",xi+n*dt,y,s,abs(s-y));

	hold on
	plot(xi+n*dt,y,'rx');
	plot(xi+n*dt,s,'bx');
    end 
	
print appvsreal00.jpg
fclose(fp);