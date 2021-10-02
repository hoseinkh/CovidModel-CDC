% Pandemic
%
%
%
close all, clear all, clc;
% clear all, clc;
%% Initial Values
S0 = 0.99;
E0 = 0.01;
I0 = 0.00;
R0 = 0;
n0 = 0.8;
%
x0 = 0.7;
%
%
Beta0 = 0.8;
alpha = 0.2; 
%
%
%
MaxTime = 1000;
tspan=[0, MaxTime];
%
%
recovRate = 0;
epsil = 0.5;
Lambda = 0.0;
mu = 1*Lambda;
eta = 0.5;
nu = 0; % percentage die from disease
rhoo = 0.3; % has no role
N0 = S0+E0+I0+R0;
%
outResultsCell = cell(1, 7)
%
for ttt = [1:7]
    %% Parameters    
    %
    %% Oscillating TOC 45-90 degree
    A_fr = [0.5  , 1;
            1.00 , 1.25]*1;   % corresponding to n=1. Cooperative defector is NE
    %
    A_fd = [3.80 , 1.00;
            2.00 , 0.20]*1; % corresponding to n=0. Cooperative dominance is NE
    %
    %
    %
    parameters = cell(1,11);
    parameters{1} = A_fr;parameters{2} = A_fd;parameters{3} = Lambda;
    parameters{4} = mu; parameters{5} = eta;parameters{6} = alpha;parameters{7} = nu;
    parameters{8} = rhoo; parameters{9} = N0; parameters{10} = epsil; 
    parameters{11} = recovRate; parameters{12} = Beta0;
    %
    y0 = [S0;E0;I0;R0;x0;n0];
    %
    %
    %
    %% Solving set of differential equations
    [t,y] = ode23(@(t,y) call_dstate_9(t,y,parameters),tspan,y0); 
    S=y(:,1); E=y(:,2); I=y(:,3); R=y(:,4); x=y(:,5); n=y(:,6);
    %
    %
    HerdImmunIndex = zeros(length(S),1);
    for hhh = [1:length(S)]
        HerdImmunIndex(hhh) = (Beta0/alpha)*S(hhh);
    end
    %
    %% Saving results
    outMatrix = zeros(7,length(t));
    outMatrix(1,:) = t';
    outMatrix(2,:) = E';
    outMatrix(3,:) = I';
    outMatrix(4,:) = R';
    outMatrix(5,:) = x';
    outMatrix(6,:) = n';
    outMatrix(7,:) = HerdImmunIndex;
    %
    %
    outResultsCell{1,ttt} = outMatrix;
    %
end
%
hold off
%
%% Plots
%
%
figure(1)
ttt = 2; %% We are plotting for Oscillating TOC 45-90 degree
outMatrix = outResultsCell{1,ttt};
t = outMatrix(1,:);
n = outMatrix(6,:);
x = outMatrix(5,:);
I = outMatrix(3,:);
R = outMatrix(4,:);
HerdImmunIndex = outMatrix(7,:);
fprintf("max of index = %2.4f \n", max(HerdImmunIndex))
%
TimeHerdImmun = -1;
TimeIndexFound = 0;
for hhh = [2:length(HerdImmunIndex)]
    if (HerdImmunIndex(hhh) < 1) & (HerdImmunIndex(hhh-1) >= 1) & (TimeIndexFound == 0)
        TimeIndexFound = 1;
        TimeHerdImmun = hhh;
        break
    end
end
ylim([0 0.8])
Beta = Beta0 * (1-x);
plot(t,I,'k','LineWidth',3); hold on;
plot(t,x,'--k','LineWidth',3); hold on;
axis([0,MaxTime,0,1.05*max([max(I),max(x)])]); hold off;
xlabel("time",'FontSize',24)
ylabel("fraction",'FontSize',24)
legend("i","x", 'FontSize',24)
title(append("epsilon = ", num2str(epsil), ", R_0 = A_{11}^{FD} = ", num2str(A_fd(1,1))));
ax = gca;
ax.FontSize = 22; 
%
HH = 5+6;
%% Function
function dydt = call_dstate_9(t,y,parameters)
    %
%     parameters = [A_fr,A_fd,b,muu,muuP,alpha1,alpha2,delta,psi,omega,sigma,gammaS,gammaA,etaS,etaA];
    A_fr=parameters{1}*1;A_fd=parameters{2}*1;Lambda=parameters{3};
    mu=parameters{4};eta=parameters{5};alpha=parameters{6};
    nu=parameters{7};rhoo=parameters{8}; N0=parameters{9};
    epsil=parameters{10};recovRate=parameters{11};Beta0=parameters{12};
    %
    %
    S=y(1); E = y(2); I=y(3); R=y(4); x=y(5); n=y(6);
    %
    A_n = (1-n)*A_fd + n*A_fr;
    avg_payoff_C = A_n(1,1)*x + A_n(1,2)*(1-x);
    avg_payoff_D = A_n(2,1)*x + A_n(2,2)*(1-x);
    %
    %
    dSdt = Lambda - 1 * Beta0 * (1 - x) * S * I + recovRate * R - mu * S;
    dEdt = 1* Beta0 * (1 - x) * S * I - (eta + mu) * E;
    dIdt = eta * E - (alpha + mu + nu) * I;
    dRdt = alpha * I - recovRate * R - mu * R;
    dxdt = x * (1-x) * (avg_payoff_C - avg_payoff_D);
    dydt = [
        %% %% %% %% %% %% %%     d S / d t     %% %% %% %% %%
        dSdt;
        %% %% %% %% %% %% %%     d E / d t     %% %% %% %% %%
        dEdt;
        %% %% %% %% %% %% %%     d I / d t     %% %% %% %% %%
        dIdt;
        %% %% %% %% %% %% %%     d R / d t     %% %% %% %% %%
        dRdt;
        %% %% %% %% %% %% %%     d x / d t     %% %% %% %% %%
        dxdt;
        %% %% %% %% %% %% %%     d n / d t     %% %% %% %% %%
        -epsil*dIdt;
    ];
end
%
%

    
    
    
    
    
    