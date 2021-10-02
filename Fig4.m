% Pandemic
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
Beta0 = 0.8;
alpha = 0.2;
%
recovRate = 0;
Lambda = 0;
mu = 1*Lambda;
eta = 1;
nu = 0; % percentage die from disease
rhoo = 0.3; % has no role
N0 = S0+E0+I0+R0;
%
List_epsil = [0.1:0.02:2.3];
outResultsCell = cell(length(List_epsil), 7);
%
for epsil = List_epsil
    for ttt = [1:7]
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
        parameters = cell(1,12);
        parameters{1} = A_fr;parameters{2} = A_fd;parameters{3} = Lambda;
        parameters{4} = mu; parameters{5} = eta;parameters{6} = alpha;parameters{7} = nu;
        parameters{8} = rhoo; parameters{9} = N0; parameters{10} = epsil; 
        parameters{11} = recovRate; parameters{12} = Beta0;
        %
        y0 = [S0;E0;I0;R0;x0;n0];
        %
        MaxTime = 2000;
        tspan=[0, MaxTime];
        %
        %% Solving set of differential equations
        [t,y] = ode23(@(t,y) call_dstate_9(t,y,parameters),tspan,y0); 
        S=y(:,1); E=y(:,2); I=y(:,3); R=y(:,4); x=y(:,5); n=y(:,6);
        %
        %
%         %% Find herd immunity
%         IndexHerdImmunity = -1;
%         for gg = [1:length(Beta)]
%             if (Beta0*S(gg)/alpha < 1) && (IndexHerdImmunity == -1) && I(gg) >= 0.01
%                 IndexHerdImmunity = gg;
%                 BetaHerdImmunity = Beta(gg);
%                 break
%             end
%         end
        %% Saving results
        outMatrix = zeros(8,length(t));
        outMatrix(1,:) = t';
        outMatrix(2,:) = E';
        outMatrix(3,:) = I';
        outMatrix(4,:) = R';
        outMatrix(5,:) = x';
        outMatrix(6,:) = n';
        %
        currIndex = find(List_epsil == epsil);
        outResultsCell{currIndex,ttt} = outMatrix;
        %
    end
    %
    hold off
end
%
%% Plots
counter = 0;
counter = counter + 1;
figure(counter); hold on;
set(gcf, 'Position',  [100, 100, 1200, 800])
currSubplotIndex = 0;
for ttt = [1:3]
    for u = [1]
        Rinf4epsil = [];
        for epsil = List_epsil
            currIndex = find(List_epsil == epsil);
            outMatrix = outResultsCell{currIndex,ttt};
            R = outMatrix(4,:);
            Rinf4epsil = [Rinf4epsil, R(end)];
        end
        if ttt == 1
            plot(List_epsil, Rinf4epsil, 'k', 'LineWidth',3)
            hold on
        elseif ttt == 2
            plot(List_epsil, Rinf4epsil, '-.k', 'LineWidth',3)
            hold on
        elseif ttt == 3
            plot(List_epsil, Rinf4epsil, ':k', 'LineWidth',3)
            hold on
        end
    end
    %
end
ax = gca;
ax.FontSize = 18; 
xlabel("\epsilon",'FontSize',26)
ylabel("r (\infty)",'FontSize',26)
legend("Dimin OSC", "OSC TOC", "Anti-Coord Game",'FontSize',16)
hold off
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

    
    
    
    
    
    