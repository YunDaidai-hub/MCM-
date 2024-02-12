% Consider a more complex ecosystem, with four species:
% Parasites-Competitors 1+ Lamprey, three levels of the food chain
clear;
clc;
close all;

% In the initial population, lamprey and competitor 2 can only choose one
init_value = [100, 50, 50, 80, 500, 0, 0]; % 初始值
% Parasite;Male-Lamprey;Female-Lamprey;Competitor 1;Producer;Competitor 2

% time range
Time = linspace(0, 500, 1000);

% solving equations
[T1, Y1] = ode45(@question_4_2_B, Time, init_value);

% visualization results
figure;
plot( ...
    T1, Y1(:, 1), "m-.",...
    T1, Y1(:, 2), "r-", ...
    T1, Y1(:, 3), "b-", ...
    T1, Y1(:, 4), "c-" ...
);
hold on;
plot(T1, Y1(:, 6), "LineStyle", "-", "Color", [0.9290 0.6940 0.1250]);
plot( ...
    T1, Y1(:, 5), "g-", ...
    T1, Y1(:, 2)+Y1(:, 3), "k-" ...
);
xlabel("Time");
ylabel("Population");
legend( ...
   "Parasites", ...
   "Lamprey-Male", ...
   "Lamprey-Female", ...
   "Competitor1", ...
   "Competitor1", ...
   "Producer", ...
   "Lamprey-Total" ...
);

function eq = question_4_2_B(T, Y) 
    %{
    init_value = [100, 50, 50, 80, 500, 0, 0]; % 初始值
    Under this parameter set, although the number of each species does not 
    converge stably, it presents periodic stability, and there are alternate 
    changes in the proportion of population with large and small periodic competitors.
    1. Keeping other parameters constant,[a14a41] = [-2e-3 -2e-3] can make 
       the model converge slowly, and decreasing [a14a41] can accelerate convergence 
       and shorten the period.
    2. Based on this model, either the heptagill or competitor species are 
       deleted, and the same amount removed is added to the initial number 
       of the other. Looking at the model, it is found that the competitor 
       has difficulty maintaining the survival of the parasite, while the 
       heptagill can provide a relatively superior living environment for the parasite.
    %}

    eq = zeros(7,1);
    % Subfunctions and Parameter Definitions
    A = (Y(5)+Y(6))./(1e-50+Y(2)+Y(3)+Y(4)); % 1e-50 Prevent denominator from zero
    % R = alpha * ln(A + beta) + gamma; (5:3, 78%) (1:5, 50%) (1:20, 40%)
    p_alpha = -0.1514;
    p_beta = 0.4585;
    p_gamma = 0.7569;
    R = p_alpha .* log(A + p_beta) + p_gamma; % Sex determines the proportion of males
    R_r = Y(2)./(1e-50+Y(2)+Y(3)); % The true male ratio of lamprey populations at a given time
    
    % Lamprey birth rate
    p_sigma = 0.72;
    birth_rate = p_sigma.*R_r.*(1-R_r);
    % Lamprey mortality
    death_rate = 0.20;
    
    f_m = [ % Transmission efficiency
        % Since the parasite natural growth rate is negative, f<0 should be set to help it survive.
        0, -0.8, -0.8, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0;
    ];
    
    a_m = [
        0, 2.5e-3, 2.5e-3, 2.5e-3, 0, 0, 0;
        -2.5e-3, 0, 0, 0, (6e-4).*(1+Y(3)./(Y(2)+1e-50)), 0, 0; % M
        -2.5e-3, 0, 0, 0, (6e-4).*(1+Y(2)./(Y(3)+1e-50)), 0, 0;
        -2.5e-3, 0, 0, 0, 1e-3, 0, 0;
        0, -6e-4, -6e-4, -1e-3, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0;
    ];

    c_m = [
        0, 0.55, 0.55, 0.55, 0, 0, 0;
        0.2, 0, 0, 0, R.*(0.4.*Y(2) + 0.38.*Y(3))./(Y(2)+Y(3)+1e-50), 0, 0;  % M
        0.2, 0, 0, 0, (1-R).*(0.4.*Y(2) + 0.38.*Y(3))./(Y(2)+Y(3)+1e-50), 0, 0;
        0.2, 0, 0, 0, 0.39, 0, 0;
        0, 1, 1, 1, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0;
        0, 0, 0, 0, 0, 0, 0;
    ];
    
    p_v = [
        Y(1);
        Y(2)+Y(3);
        Y(2)+Y(3);
        Y(4);
        Y(5);
        Y(6);
        Y(7);
    ];

    r_v = [
       -0.16;
       (1+Y(3)./(Y(2)+1e-50)).*R.*(birth_rate) - death_rate;
       (1+Y(2)./(Y(3)+1e-50)).*(1-R).*(birth_rate) - death_rate;
       -0.1;
       0.1;
       0;
       0;
    ];  % Natural Growth Rate
    
    % k_m = [1e6, 5e3, 5e3, 1e4, 1e5, 1e5, 1e8]'; 
    k_v = [ % k2 == k3
        5e3;    
        5e3;    % M
        5e3;    % F 
        4e3;
        5e3;    % Producer
        1e-9;
        1e-9;
    ];  % Environmental carrying capacity

    eq(1) = r_v(1).*Y(1).*(1-(p_v(1)-f_m(1,:)*Y)./k_v(1)) + Y(1).*(a_m(1,:).*c_m(1,:))*Y;
    eq(2) = r_v(2).*Y(2).*(1-(p_v(2)-f_m(2,:)*Y)./k_v(2)) + Y(2).*(a_m(2,:).*c_m(2,:))*Y;
    eq(3) = r_v(3).*Y(3).*(1-(p_v(3)-f_m(3,:)*Y)./k_v(3)) + Y(3).*(a_m(3,:).*c_m(3,:))*Y;
    eq(4) = r_v(4).*Y(4).*(1-(p_v(4)-f_m(4,:)*Y)./k_v(4)) + Y(4).*(a_m(4,:).*c_m(4,:))*Y;
    eq(5) = r_v(5).*Y(5).*(1-(p_v(5)-f_m(5,:)*Y)./k_v(5)) + Y(5).*(a_m(5,:).*c_m(5,:))*Y;
    eq(6) = r_v(6).*Y(6).*(1-(p_v(6)-f_m(6,:)*Y)./k_v(6)) + Y(6).*(a_m(6,:).*c_m(6,:))*Y;
    eq(7) = r_v(7).*Y(7).*(1-(p_v(7)-f_m(7,:)*Y)./k_v(7)) + Y(7).*(a_m(7,:).*c_m(7,:))*Y;
end