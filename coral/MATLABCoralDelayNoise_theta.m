fp=fopen('Trials.csv','w');
fprintf(fp,'dt,beta,g,tau,trial,theta,initMacro,initCoral,initTurf,macroalgae,coral,turf,lgMacro,lgCoral,lgTurf');

TAU_START = 0.2;
TAU_END = 0.6;
NUMBER_TAU = 3;

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
FINAL_TIME = 35.0;
BASE_DT = 0.00001;
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

randomNumbers = normrnd(0,sqrt(dt));
steps = 0;    %The number of steps to take in a single simulation.
reefState = zeros(1,2);    %The state of the system, coral and macroalgae
macroalgaePast = 0;    %Past values for the macroalgae density.
coralPast = 0;         %Past values for the coral density.

        
for tauStep = 0.2:0.2:NUMBER_TAU
    
    tau = TAU_START+(TAU_END - TAU_START)/((NUMBER_TAU - 1)*(tauStep));
    
    %Determine the number of time steps required to move back to the delay in time.
    %The number of cells needed for the delay (changes with dt)
    
    if tau > 0.0
       n = int32(tau/(BASE_DT+0.5));
    end
    
    %Allocate the space for the states of the system
    macroalgaePast = calloc(n,sizeof(double));		%macroalgae density in the past
    coralPast      = calloc(n,sizeof(double));		%coral density in the past

    if(isempty(macroalgaePast) || isempty(coralPast))
       message = 'Error - unable to allocate necessary memory.'
       free(macroalgaePast);
       free(coralPast);
    end
    
    steps = FINAL_TIME/dt;

    for theta = {9 , 10 ,12} 
        %Make an approximation for different initial conditions.
        %Make an arc through 0 to pi/2 radians from the origin.

        for k = 0:NUMBER_TRIALS
            macroalgaePast(1,1) = RADIUS*cos(theta); %initial Macroalgae level
            coralPast(1,1)      = RADIUS*sin(theta); %initial Coral level
            
            for l = 1:n %//fills in the past times for y, z, v, and w
                macroalgaePast(1,2) = macroalgaePast(1,1);
                coralPast(1,2) = coralPast(1,1);
            end
            
            m = n-1;
            p = 0;
            calcRandom = 0;
            
            %Step through every time step.
            for c =0:steps
                %Update the macroalgae term
                reefState(1,1) =  macroalgaePast(1,m) + (macroalgaePast(1,m)*(gamma-gamma*macroalgaePast(1,m)+(a-gamma)*coralPast(1,m))-(g*macroalgaePast(1,p)/(1.0-coralPast(1,p))))*dt + beta*macroalgaePast(1,m)*(1.0-macroalgaePast(1,m))*B(1,calcRandom) + 0.5*beta*(1.0-2.0*macroalgaePast(1,m))*beta*macroalgaePast(1,m)*(1-macroalgaePast(1,m))*(B(1,calcRandom)*B(1,calcRandom)-dt);
                %Computes the deterministic component for Coral
                reefState(1,1) =  coralPast(1,m)+(coralPast(1,m)*(r-d-(a+r)*macroalgaePast(1,m)-r*coralPast(1,m)))*dt;

%               ****************************************************************
%                   Account for extinction and overgrowing!!
%               ****************************************************************
                for i = 0:2
                    if reefState(1,i) < 0.0
                       reefState(1,i) = 0.0;
                    elseif reefState(1,i) > 1.0
                       reefState(1,i) = 1.0;
                    end
                end

            %Updates delay and cell index
            m = mod((m+1),n);
            p = mod((p+1),n);
            macroalgaePast(1,m) = reefState(1,1);
            coralPast(1,m) = reefState(1,1);
            calcRandom = mod((calcRandom+1),2); %update which random number to use.
            end %for(c<steps)                  
        end %for(k<trials)
    end %for(theta<pi/2

    %Free up the allocated memory.
    free(macroalgaePast);
    free(coralPast);

    fp.close();


end %tau = {9,10,12}

                     

