fp=fopen('Trials.csv','w');
fprintf(fp,'dt,beta,g,tau,trial,theta,initMacro,initCoral,initTurf,macroalgae,coral');

TAU_START = 0.2;
TAU_END = 0.6;
NUMBER_TAU = 2;

G_START = 0.2;
G_END = 0.6;
NUMBER_G = 3;

BETA_START = 0.2;
BETA_END = 1.2;
NUMBER_BETA = 6;

THETA_START = 0.0;
THETA_END = pi*0.5;
NUMBER_THETA = 15;
NUMBER_TRIALS = 5;
FINAL_TIME = 5.0;
BASE_DT = 0.01;
RADIUS = 0.06;

%Define the constants.
%These are the parameters that are used in the operator for the differential equations.
        
a = 0.1;
gamma = 0.8;
r = 1.0;
d = 0.44;
%tau; %= 0.5;
g  = 0.6;
beta = 0.4;
dt = BASE_DT; %BASE_DT; Set the initial time step
theta = 0;


steps = 0;    %The number of steps to take in a single simulation.
reefState = zeros(1,2);    %The state of the system, coral and macroalgae
macroalgaePast = 0;    %Past values for the macroalgae density.
coralPast = 0;         %Past values for the coral density.

        
for tauStep = 0:1:(NUMBER_TAU-1)
    
    tau = TAU_START+(TAU_END - TAU_START)/(NUMBER_TAU - 1)*(tauStep);
    
    %Determine the number of time steps required to move back to the delay in time.
    %The number of cells needed for the delay (changes with dt)
    n = 1;
    if tau > 0.0
       n = int32(tau/BASE_DT+0.5);
    end
    
    %Allocate the space for the states of the system
    macroalgaePast = zeros(n,1);		%macroalgae density in the past
    coralPast      = zeros(n,1);		%coral density in the past

    steps = FINAL_TIME/dt;

    for theta = [10, 12] * pi/80
        %Make an approximation for different initial conditions.
        %Make an arc through 0 to pi/2 radians from the origin.

        for k = 1:NUMBER_TRIALS
            macroalgaePast(1) = RADIUS*cos(theta); %initial Macroalgae level
            coralPast(1)      = RADIUS*sin(theta); %initial Coral level
            icMacro = macroalgaePast(1);
            icCoral = coralPast(1);
            
            for l = 1:n %//fills in the past times for y, z, v, and w
                macroalgaePast(l) = macroalgaePast(1);
                coralPast(l) = coralPast(1);
            end
            
            m = n;
            p = 1;
            
            %Step through every time step.
            for c = 1:steps
                %Update the macroalgae term
                B = normrnd(0, sqrt(dt));
                reefState(1) =  macroalgaePast(m) + ...
                    (macroalgaePast(m)*(gamma-gamma*macroalgaePast(m)+(a-gamma)*coralPast(m))-(g*macroalgaePast(p)/(1.0-coralPast(p))))*dt ...
                    + beta*macroalgaePast(m)*(1.0-macroalgaePast(m))*B...
                    + 0.5*beta*(1.0-2.0*macroalgaePast(m))*beta*macroalgaePast(m)*(1-macroalgaePast(m))*(B*B-dt);
                %Computes the deterministic component for Coral
                reefState(2) =  coralPast(m)+(coralPast(m)*(r-d-(a+r)*macroalgaePast(m)-r*coralPast(m)))*dt;

%               ****************************************************************
%                   Account for extinction and overgrowing!!
%               ****************************************************************
                for i = 1:2
                    if reefState(i) < 0.0
                       reefState(i) = 0.0;
                    elseif reefState(i) > 1.0
                       reefState(i) = 1.0;
                    end
                end
            
            %Updates delay and cell index
            m = 1 + mod(m,n);
            p = 1 + mod(p,n);
            macroalgaePast(m) = reefState(1);
            coralPast(m) = reefState(2);
            
            fprintf(fp, '%f,%f,%f,%f,%d,%f,%f,%f,%f,%f,%f\n',...
            dt, beta, g, tau, k, theta, icMacro, icCoral, 1-icMacro-icCoral, reefState(1), reefState(2));
            
            dt, beta, g, tau, k,c
            input('enter');
            end %for(c<steps)                  
        end %for(k<trials)
    end %for(theta<pi/2

end %tau = {9,10,12}

fp.close();                     

