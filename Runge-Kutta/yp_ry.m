fp=fopen('RK_y_ry_02.csv','w');
fprintf(fp,"Approximation,r,steps\n");
r=-3;
a=1/2;
b=1/2;
bet=1;
alfa=1;

final=1;
steps=10;
dt=final/steps;
fprintf(fp,"1,%d,%d\n",r,steps);

n=0;
y=1;

while n<steps
	k1=dt*r*y;
	k2=dt*r*(y+bet*k1);
	y=y+a*k1+b*k2;
	n=n+1;
	fprintf(fp,"%f\n",y);
	end 
y
s=e^(r)
s-y
fclose(fp)