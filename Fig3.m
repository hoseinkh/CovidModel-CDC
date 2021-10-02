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
Beta0 = 1;
alpha = 0.2;
%
%
recovRate = 0;
epsil = 2;
Lambda = 0;
mu = 1*Lambda;
eta = 1;
nu = 0; % percentage die from disease
rhoo = 0.3; % has no role
N0 = S0+E0+I0+R0;
%
for ttt = [1:7]
    ttt;
    %% Parameters    
    if ttt == 1
        %% Averted 0-45 degree
        A_fr = [0.5  , 1;
                1.00 , 1.25]*1;   % corresponding to n=1. Cooperative defector is NE
        %
        A_fd = [2.50 , 1.00;
                  2.00 , 0.25]*1; % corresponding to n=0. Cooperative dominance is NE
%
    elseif ttt == 2
        %% Oscillating TOC 45-90 degree
        A_fr = [0.5  , 1;
                1.00 , 1.25]*1;   % corresponding to n=1. Cooperative defector is NE
        %
        A_fd = [3.90 , 1.00;
                  2.00 , 0.25]*1; % corresponding to n=0. Cooperative dominance is NE
%
    elseif ttt == 3
        %% Averted 315-360 degree
        A_fr = [0.5  , 1;
                1.00 , 1.25]*1;   % corresponding to n=1. Cooperative defector is NE
        %
        A_fd = [3.00 , 1.00;
                  3.20 , 0.25]*1; % corresponding to n=0. Cooperative dominance is NE
%
    elseif ttt == 4
        %% TOC 270-315 degree
        A_fr = [0.5  , 1;
                1.00 , 1.25]*1;   % corresponding to n=1. Cooperative defector is NE
        %
        A_fd = [2.50 , 1.00;
                  4.50 , 0.25]*1; % corresponding to n=0. Cooperative dominance is NE
%
    elseif ttt == 5
        %% TOC 180-270 degree
        A_fr = [0.5  , 1;
                1.00 , 1.25]*1;   % corresponding to n=1. Cooperative defector is NE
        %
        A_fd = [3.00 , 1.00;
                  3.20 , 1.25]*1; % corresponding to n=0. Cooperative dominance is NE
%
    elseif ttt == 6
        %% TOC 135-180 degree
        A_fr = [0.5  , 1;
                1.00 , 1.25]*1;   % corresponding to n=1. Cooperative defector is NE
        %
        A_fd = [3.00 , 1.00;
                  2.80 , 1.75]*1; % corresponding to n=0. Cooperative dominance is NE
%
    elseif ttt == 7
    %% TOC 90-135 degree
        A_fr = [0.5  , 1;
                1.00 , 1.25]*1;   % corresponding to n=1. Cooperative defector is NE
            %
            A_fd = [2.70 , 1.00;
                    2.00 , 1.25]*1; % corresponding to n=0. Cooperative dominance is NE
    end
    %
    %
    Lambda = 0.;
    mu = 1*Lambda;
    eta = 1;
    alpha = 0.2; % 0.50444;
    nu = 0; % percentage die from disease
    rhoo = 0.3;
    N0 = S0+E0+I0+R0;
    %
    parameters = cell(1,11);
    parameters{1} = A_fr;parameters{2} = A_fd;parameters{3} = Lambda;
    parameters{4} = mu; parameters{5} = eta;parameters{6} = alpha;parameters{7} = nu;
    parameters{8} = rhoo; parameters{9} = N0; parameters{10} = epsil; 
    parameters{11} = recovRate; parameters{12} = Beta0;
    %
    y0 = [S0;E0;I0;R0;x0;n0];
    %
    MaxTime = 500;
    tspan=[0, MaxTime];
    %
    %% Solving set of differential equations
    [t,y] = ode23(@(t,y) call_dstate_9(t,y,parameters),tspan,y0); 
    S=y(:,1); E=y(:,2); I=y(:,3); R=y(:,4); x=y(:,5); n=y(:,6);
    %
    %
%     %% Find herd immunity
%     IndexHerdImmunity = -1;
%     for gg = [1:length(Beta)]
%         if (Beta(gg)*S(gg)*I(gg)/alpha < 1) && (IndexHerdImmunity == -1)
%             IndexHerdImmunity = gg;
%             BetaHerdImmunity = Beta(gg);
%             break
%         end
%     end
    %% Plots
    figure(ttt)
    plot(x,I, 'b','LineWidth',2)
    ax = gca;
    ax.FontSize = 18; 
    xlabel("x",'FontSize',24,'FontWeight','bold')
    ylabel("i",'FontSize',24,'FontWeight','bold')
%
    axis([min(x)*0.9,max(x)*1.1,min(I)*0.9,max(I)*1.1])
    if ttt == 1
        title("Averted 0-45 degree")
        fprintf("Averted 0-45 degree")
        fprintf("R inf = %2.2f \n",R(end))
    elseif ttt == 2
        title("Oscillating TOC 45-90 degree")
        fprintf("Oscillating TOC 45-90 degree")
        fprintf("R inf = %2.2f \n",R(end))
    elseif ttt == 3
        title("Averted 315-360 degree")
        fprintf("Averted 315-360 degree")
        fprintf("R inf = %2.2f \n",R(end))
    elseif ttt == 4
        title("TOC 270-315 degree")
        fprintf("TOC 270-315 degree")
        fprintf("R inf = %2.2f \n",R(end))
    elseif ttt == 5
        title("TOC 180-270 degree")
        fprintf("TOC 180-270 degree")
        fprintf("R inf = %2.2f \n",R(end))
    elseif ttt == 6
        title("TOC 135-180 degree")
        fprintf("TOC 135-180 degree")
        fprintf("R inf = %2.2f \n",R(end))
    elseif ttt == 7
        title("TOC 90-135 degree")
        fprintf("TOC 90-135 degree")
        fprintf("R inf = %2.2f \n",R(end))
    end
    %
end
%
%
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

    
    
    
    
    
    