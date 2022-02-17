% This script contains the code for generating the plot threshold 
clear;
% COVID-19 DATA
% Time unit = 1day

beta = 0.52;
gamma = 1/10;       % 1/gamma = days of recovery
tau_Q = 1/14;       % 1/tau_Q = days of quarantine
TPR = 0.69;          % true positive ratio for Bluetooth
FPR = 0.45;          % False positive ratio for Bluetooth
k = 8;
R_0 = 3;        % Reproductive ratio
b = R_0*gamma/k;  % We obtain b from the reproductive ratio (see eq. in book R0 = k*b/gamma;)

delta1 = 0.05;      % rate of detecting and isolating infected individuals
tau_T = 1/1;      % 1/tau_T = Trace time (days)

figure;



FN=0.5;


N = 1e6;
S = N;
V = 0;
v = 0.8;

S_N = 1;
cT = 0.05:0.05:1;
SRN = (1-cT)*S_N;
delta = 0.01;
Re1 = (delta+gamma)./(gamma.*SRN);

S_N = 1;
cT = 0.05:0.05:1;
SRN = (1-cT)*S_N;
delta = 0.05;
Re2 = (delta+gamma)./(gamma.*SRN);

S_N = 0.75;
cT = 0.05:0.05:1;
SRN = (1-cT)*S_N;
delta = 0.01;
Re3 = (delta+gamma)./(gamma.*SRN);

S_N = 0.75;
cT = 0.05:0.05:1;
SRN = (1-cT)*S_N;
delta = 0.05;
Re4 = (delta+gamma)./(gamma.*SRN);

S_N = 0.5;
cT = 0.05:0.05:1;
SRN = (1-cT)*S_N;
delta = 0.01;
Re5 = (delta+gamma)./(gamma.*SRN);

S_N = 0.5;
cT = 0.05:0.05:1;
SRN = (1-cT)*S_N;
delta = 0.05;
Re6 = (delta+gamma)./(gamma.*SRN);


plot(Re1, cT, 'k','LineWidth',2);
hold on;
plot(Re2, cT, 'r','LineWidth',2);
hold on;
plot(Re3, cT, 'k--','LineWidth',2);
hold on;
plot(Re4, cT, 'r--','LineWidth',2);
hold on;
plot(Re5, cT, 'k-.','LineWidth',2);
hold on;
plot(Re6, cT, 'r-.','LineWidth',2);
hold on;
% plot(U*100,RC4*100,'r:s','LineWidth',2);
% hold on;
% plot(U*100,RC5*100,'b-^','LineWidth',2);
% hold on;
% plot(U*100,RC6*100,'r:v','LineWidth',2);

% text(0.09,0.25,'S/N=1','FontSize',16)
% text(0.063,0.2,'S/N=0.8','FontSize',16)
% text(0.038,0.15,'S/N=0.6','FontSize',16)
% text(0.01,0.1,'S/N=0.4','FontSize',16)


ylabel('Ratio of contacts traced (c_T)');
xlabel('Reproductive ratio (R_e)'); 
set(gca,'FontSize',20);
legend('Im=0,\delta=0.1', 'Im=0,\delta=0.5', 'Im=0.25,\delta=0.1','Im=0.25,\delta=0.5', 'Im=0.5,\delta=0.1','Im=0.5,\delta=0.5');

xlim([0 16]);
xticks([0:2:16]);
ylim([0 1]);
