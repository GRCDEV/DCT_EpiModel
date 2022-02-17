% This files contains the code to generate all the plots of the paper
% using the SIR_Trace_withVacc_Euler function.

% Time unit = 1day
clear;

Y1 = 1;
R1 = 0;
N = 50e6;
beta = 0.52;
gamma = 1/10;       % 1/gamma = days of recovery
tau_Q = 1/14;       % 1/tau_Q = days of quarantine
tau_T = 1;          % 1/tau_T = Trace time (days)

TPR = 0.69;         % true positive ratio for Bluetooth
FPR = 0.45;         % False positive ratio for Bluetooth

IFR = 0.01;         % Infection fatality rate

MAX_T = 365*2-60; %  TWO year

% Measures
R_0 =  3;           % Reproductive ratio
                    % Beta = k*b
k = 8;
b = R_0*gamma/k;    % We obtain b from the reproductive ratio (see eq. in book R0 = k*b/gamma;)

% NO MEASURES
K = k*ones(1,MAX_T);
Re = R_0*ones(1,MAX_T); 

delta = 0.01;       % rate of detecting and isolating infected individuals

% Vaccination
omega = 0.004;
OMEGA = omega*ones(1,MAX_T);
OMEGA(1:360) = zeros(1,360); 
v = 0.9;


% Measure periods (default)
m_periods = [ 75, 160; 183, 250; 274, 305; 330, 400];


% TYPE OF GRAPH
NO_DCT = 1;  % FIGURE 3
DCT = 2;   % FIGURES 4 AND 6,7

%--------------------
% TYPE OF FIGURES
%--------------------

%%%%
FIGURE = 4;   % ---> FIGURE TO BE GENERATED.
%%%%

if FIGURE == 3
    TYPE_OF_GRAPH = NO_DCT; 
    AR = zeros(1,MAX_T);            % Utilisation
    TC = 0.0;            % Ratio of users checking
end
if FIGURE == 4   % GOOGLE/APPLE API
    % Figures 4a,4b
    TYPE_OF_GRAPH = DCT; 
    tau_T = 1/2;         % 1/tau_T = Trace time (days)
    AR = zeros(1,MAX_T); 
    UTILISATION =  0.28;
    AR((8*30+1):12*30) = UTILISATION*ones(1,12*30-8*30);           % Utilisation
    TC = 0.8;            % Ratio of users checking
    TPR = 0.69;          % true positive ratio for Bluetooth
    FPR = 0.45;          % False positive ratio for Bluetooth
end
if FIGURE == 43   % GOOGLE/APPLE API
    % Figure 4c
    TYPE_OF_GRAPH = DCT; 
    m_periods = [ 75, 160; 183, 240; 390, 400];
    tau_T = 1/2;         % 1/tau_T = Trace time (days)
    AR = zeros(1,MAX_T); 
    UTILISATION = 0.28;
    AR((8*30+1):12*30) = UTILISATION*ones(1,12*30-8*30);           % Utilisation
    TC = 0.8;            % Ratio of users checking
    TPR = 0.69;          % true positive ratio for Bluetooth
    FPR = 0.45;          % False positive ratio for Bluetooth
end
if FIGURE == 5    % CHINESE/KOREAN APP
    % Figures 5a,5b
    TYPE_OF_GRAPH = DCT; 
    tau_T = 1;         % 1/tau_T = Trace time (days)
    AR = zeros(1,MAX_T); 
    UTILISATION =  0.8;
    AR((8*30+1):12*30) = UTILISATION*ones(1,12*30-8*30);           % Utilisation
    TC = 1;            % Ratio of users checking
    TPR = 0.5;          % true positive ratio for Bluetooth
    FPR = 0.4;          % False positive ratio for Bluetooth
end
if FIGURE == 53   % CHINESE/KOREAN APP
    % Figure 5c
    TYPE_OF_GRAPH = DCT; 
    m_periods = [ 75, 160; 183, 240; 390, 400];
    tau_T = 1;         % 1/tau_T = Trace time (days)
    AR = zeros(1,MAX_T); 
    UTILISATION = 0.8;
    AR((8*30+1):12*30) = UTILISATION*ones(1,12*30-8*30);           % Utilisation
    TC = 1;            % Ratio of users checking
    TPR = 0.5; % true positive ratio for Bluetooth
    FPR = 0.4;          % False positive ratio for Bluetooth
end


if FIGURE == 7
    % Figures 7a,7b
    TYPE_OF_GRAPH = DCT; 
    tau_T = 1/2;         % 1/tau_T = Trace time (days)
    AR = zeros(1,MAX_T); 
    UTILISATION = 0.7;
    AR((8*30+1):12*30) = UTILISATION*ones(1,12*30-8*30);           % Utilisation
    TC = 0.8;            % Ratio of users checking
    TPR = 0.69;          % true positive ratio for Bluetooth
    FPR = 0.45;          % False positive ratio for Bluetooth
end

if FIGURE == 73
    % Figure 7c
    TYPE_OF_GRAPH = DCT; 
    tau_T = 1/2;         % 1/tau_T = Trace time (days)
    AR = zeros(1,MAX_T); 
    UTILISATION = 0.7;
    AR((8*30+1):12*30) = UTILISATION*ones(1,12*30-8*30);           % Utilisation
    TC = 0.8;            % Ratio of users checking
    TPR = 0.8;          % true positive ratio for Bluetooth
    FPR = 0.1;          % False positive ratio for Bluetooth
end

if FIGURE == 8
    % Figures 8a y 8b
    TYPE_OF_GRAPH = DCT; 
    m_periods = [ 75, 160; 183, 240];
    tau_T = 1;         % 1/tau_T = Trace time (days)
    AR = zeros(1,MAX_T); 
    UTILISATION = 0.7;
    AR((8*30+1):12*30) = UTILISATION*ones(1,12*30-8*30);           % Utilisation
    TC = 0.8;            % Ratio of users checking
    TPR = 0.8;          % true positive ratio for Bluetooth
    FPR = 0.1;          % False positive ratio for Bluetooth
end
    
    

% WITH MEASURES
Re(m_periods(1,1):end) = 0.75*Re(m_periods(1,1):end);
measures = 2*ones(1,MAX_T);
measures(1:m_periods(1,1))= zeros(1,m_periods(1,1));
T_m = (1:MAX_T)./30;
for i = 1:size(m_periods, 1)
    m_i = m_periods(i,1);
    m_e = m_periods(i,2);
    K(m_i:m_e) = K(m_i:m_e)*0.4; % 0.1;
    Re(m_i:m_e) = Re(m_i:m_e)*0.40; % 0.05;
    measures(m_i:m_e) = 2*measures(m_i:m_e);
end
%Re(m_periods(end,2):end) = 0.85*Re(m_periods(end,2):end);
B = Re.*gamma./K;


cT = TPR*AR.^2*TC;
cA = FPR*AR.^2*TC;
dt = 0.1;

[ S, I, R, V, Q_S, Q_I, Q_T, Qa_T, T, tEndInfection ] = SIR_Trace_withVacc_Euler(N, Y1, R1, K, B, cT, cA, tau_Q, delta, gamma,  OMEGA, v, tau_T, MAX_T, dt);


p_i = 8*30*10;
p_e = 12*30*10;
fprintf('Test end %7.2f  End of infection %.2f (days), Recovered %.2f, Deaths %.2f, Qa_T = %7.2f\n', S(end)+I(end)+V(end)+Q_S(end)+Q_I(end)+Q_T(end)+R(end), tEndInfection, R(end), R(end)*IFR, Qa_T(end));
fprintf('PERIOD_ 9-12 MONTHS Recovered %.2f, Deaths %.2f, Qa_T = %7.2f\n', R(p_e)-R(p_i), (R(p_e)-R(p_i))*IFR, Qa_T(p_e)-Qa_T(p_i));


Pop_Unit = 1000000;

T = T/30;
tEndInfection = tEndInfection/30;


if TYPE_OF_GRAPH == NO_DCT
% Show vaccination

T = T + 0.5;
T_m = T_m + 0.5;

figure;
plot(T,S/Pop_Unit,'k','LineWidth',2);
hold on;
plot(T,I/Pop_Unit,'r','LineWidth',2);
hold on;
plot(T,R/Pop_Unit,'b','LineWidth',2);
hold on;
plot(T,V/Pop_Unit-0.1,'color',[0.4660 0.6740 0.1880],'LineWidth',2); % Resto 0.1 para que no salga antes 
hold on;
% plot(T,Q_S/Pop_Unit,'k--','LineWidth',2);
% hold on;
% plot(T,Q_I/Pop_Unit,'r--','LineWidth',2);
% hold on;
% plot(T,Q_T/Pop_Unit,'m--','LineWidth',2);
% hold on;
plot(T_m,measures*3,':','color',[0.4940 0.1840 0.5560],'LineWidth',2);

xlabel('months');
ylabel('Population (millions)');
%legend('S','I','R','V','Q_S','Q_I','Q_T','Meas.','Location','northeast','NumColumns',2);
legend('S','I','R','V','Meas.','Location','northeast','NumColumns',2);

set(gca,'FontSize',22);
xlim([0+0.5 tEndInfection+0.5]);
ylim([0 N/Pop_Unit]);
xticks(1:2:23);

end

if TYPE_OF_GRAPH == DCT

T = T + 0.5;
T_m = T_m + 0.5;
    
figure;
plot(T,I/Pop_Unit,'r','LineWidth',2);
hold on;
plot(T,R/Pop_Unit,'b','LineWidth',2);
hold on;
plot(T_m,measures,':','color',[0.4940 0.1840 0.5560],'LineWidth',2);

xlabel('months');
ylabel('Population (millions)');
legend('I','R','Meas.','Location','northwest','NumColumns',1);

set(gca,'FontSize',22);
xlim([8.90 12.99]);
%ylim([0 12]);
xlim([8.5 12.5]);
xticks(9:12);

figure;
plot(T,Q_S/Pop_Unit,'k','LineWidth',2);
hold on;
plot(T,Q_I/Pop_Unit,'r','LineWidth',2);
hold on;
plot(T,Q_T/Pop_Unit,'m','LineWidth',2);


xlabel('months');
ylabel('Population (millions)');
legend('Q_S','Q_I','Q_T','Location','northwest','NumColumns',1);

set(gca,'FontSize',22);
xlim([8.5 12.5]);
xticks(9:12);

end






