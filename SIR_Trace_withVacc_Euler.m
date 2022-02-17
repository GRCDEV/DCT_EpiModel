function [ S, I, R, V, Q_S, Q_I, Q_T, Qa_T, T, tEndInfection ] = SIR_Trace_withVacc_Euler(N, I0, R0, K, B, CT, CA, tau_Q, delta, gamma, OMEGA, v, tau_T, MAX_T, dt)
% Resolution of the SIR model considering vaccination, using three types of
% quarantines; solved numerical using the Euler method.
% PARAMETERS
%   N: Population 
%   I0: initial infected population
%   R0: initial recovered population (inmunised)
%   K: contact rate - number of contact per unit time (i.e. Lambda) -
%   depending on time
%   B: probability of transmitting the disease (NOTE that beta = k*b) -
%   depending on time
%   CT: True traced contacts ratio
%   CA: False alerts ratio
%   tau_Q: 1/tau_Q is the average time spent in isolation
%   delta: rate of detecting and isolating infected individuals
%   gamma: removal or recovery rate
%   omega: vaccination rate (per day)
%   v: weighted efficay of the vaccines
%   tau_T: time of tracing.
%   MAX_T: Simulation time (s)-> If MAX_T is negative then only returns the 
%       values when there is infection (I(t) > 1) the simulation
% RETURNS:
%   S: Susceptible population 
%   I: Infected population 
%   R: Recovered population
%   V: Vaccinated population
%   Q_S: Susceptible individuals in quarantine
%   Q_I: Infected individuals in quarantine and trace
%   Q_T: Infected quarantined being traced.
%   Qa_T: Accumulated individuals forced to be quarantined by trace
%   T: Time vector 
%   tEndInfection: Time of end of infection.

n = ceil(MAX_T/dt);
S = zeros(1,n+1);
I = zeros(1,n+1);
R = zeros(1,n+1);
V = zeros(1,n+1);
Q_S = zeros(1,n+1);
Q_I = zeros(1,n+1);
Q_T = zeros(1,n+1);
Qa_T = zeros(1,n+1);

I(1) = I0; R(1) = R0;  S(1) = N - I0 - R0;
End_RQ= false;
% Normalise q
CT = CT*tau_T;
CA = CA*tau_T;
tau_Qi = 1/(1/tau_Q-1/tau_T);


T = 0:dt:MAX_T;
% Define the ODEs of the model and solve numerically by Euler's method:
tEndInfection = MAX_T;
for i = 1:n-1
    k = K(floor(i*dt)+1); b = B(floor(i*dt)+1); omega = OMEGA(floor(i*dt)+1);
    cT = CT(floor(i*dt)+1); cA = CA(floor(i*dt)+1);
    S(i+1) = S(i) + dt*( -(1-cT)*(k*b*I(i)*S(i))/N -cT*(k*b*I(i)*S(i))/N + ...     
             -cA*(k*(1-b)*Q_I(i)*S(i))/N - omega*S(i) + tau_Q*Q_S(i)); 
    I(i+1) = I(i) + dt*((1-cT)*(k*b*I(i)*S(i))/N + (1-v)*(k*b*I(i)*V(i))/N + ...
              - delta*I(i) - gamma*I(i));
    R(i+1) = R(i) + dt*(gamma*I(i)+tau_Q*Q_I(i));    
    V(i+1) = V(i) + dt*(omega*S(i) - (1-v)*(k*b*I(i)*V(i))/N);
    Q_S(i+1) = Q_S(i) + dt*(cA*(k*(1-b)*Q_I(i)*S(i))/N - tau_Q*Q_S(i));
    Q_I(i+1) = Q_I(i) + dt*(tau_T*Q_T(i) - tau_Qi*Q_I(i));
    Q_T(i+1) = Q_T(i) + dt*(delta*I(i) + cT*(k*b*I(i)*S(i))/N - tau_T*Q_T(i));
    Qa_T(i+1) = Qa_T(i) + dt*(cT*(k*b*I(i)*S(i))/N +cA*(k*(1-b)*I(i)*S(i))/N);
    
    if (I(i+1)+Q_I(i+1)+Q_T(i+1)) < 1
        tEndInfection = i*dt;
        break;
    end    
    iEnd = i;

end

S = S(1:iEnd); 
I = I(1:iEnd);
R = R(1:iEnd); 
V = V(1:iEnd);
Q_S = Q_S(1:iEnd);
Q_I = Q_I(1:iEnd);
Q_T = Q_T(1:iEnd);
Qa_T = Qa_T(1:iEnd); 
T = T(1:iEnd);

end
