function [dt error]=aproximacion_y_y(w,x)
	n=0;
	xi=0;
	y=1;
	
	final=w;
	steps=x;
	dt=final/steps;
		
	while n<steps
		y=y*(1+dt+0.5*dt^2);
		n=n+1;
		s=e^((xi+n*dt));
		error=abs(s-y);
		end
	endfunction
clf
fp=fopen('error_dt_f1.csv','w');
fprintf(fp,"dt,Error\n");
for i=10:100
	[r1 r2]=aproximacion_y_y(1,i);
	fprintf(fp,"%f,%f\n",r1,r2);
	i=i+10;
	
	hold on 
	plot(r1,r2,'.');
	endfor
fclose(fp);
title('Error in y(1)')
xlabel('dt')
ylabel('error')
hold off
print error_y1.jpg



