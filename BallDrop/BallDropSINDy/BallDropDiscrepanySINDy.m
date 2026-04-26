clear;

BallDrop;
%% Constants
v_wDiscrepancy(1) = 0; % [m/s]
s_wDiscrepancy(1) = 0; % [m]
Rho = 1.225; % [kg/m^3]
mu = 1.81e-5; % [Pa*s]
D = 0.021963 * 2; % [m] 
A = D^2 * pi() / 4; % [m^2]
Re(1) = 0;
Cd(1) = 0;
a_wDiscrepancy(1) = (m * g) / m;
discrepancy = 0.05;

%% Equations with an additional v added to F
% a = -Fd / m + (m * g) / m + discrepancy * v

% Fd = (1/2) * Rho * v^2 * Cd(Re) * A
% Re = (Rho * v * D) / mu

% a = -(1/2 * Rho * v^2 * Cd(Re) * A) / m + (m * g) / m + discrepancy * v

for i = 1 : length(t) - 1
    v_wDiscrepancy(i+1) = v_wDiscrepancy(i) + a_wDiscrepancy(i) * dt;
    s_wDiscrepancy(i+1) = s_wDiscrepancy(i) + v_wDiscrepancy(i) * dt;
    
    Re(i+1) = (Rho * v_wDiscrepancy(i+1) * D) / mu;
    Cd(i+1) = 24/Re(i+1) ...
    + (2.6.*(Re(i+1)/5.0)) / (1 + (Re(i+1)/5.0).^1.52) ...
    + (0.411.*(Re(i+1)/(2.63e5))^(-7.94)) / (1 + (Re(i+1)/(2.63e5))^(-8.00)) ...
    + (0.25.*(Re(i+1)/10e6)) / (1 + (Re(i+1)/10e6));
    
    a_wDiscrepancy(i+1) = -(1/2 * Rho * v_wDiscrepancy(i+1)^2 * Cd(i+1) * A) / m + (m * g) / m + discrepancy * v_wDiscrepancy(i+1);
end

%% Next step is to identify the discrepany between this model and the ideal model with SINDy



for i = 1 : length(a_ideal)
    Theta(i,1:10) = [1, v_wDiscrepancy(i), s_wDiscrepancy(i),...
    v_wDiscrepancy(i)*s_wDiscrepancy(i), v_wDiscrepancy(i)^2,...
    s_wDiscrepancy(i)^2, s_wDiscrepancy(i)^2*v_wDiscrepancy(i),...
    s_wDiscrepancy(i)*v_wDiscrepancy(i)^2, s_wDiscrepancy(i)^3, v_wDiscrepancy(i)^3];
end

for i = 1 : length(s_ideal)
    x(i, 1) = s_wDiscrepancy(i) - s_ideal(i);
    x(i, 2) = v_wDiscrepancy(i) - v_ideal(i);
end

for i = 1 : length(v_ideal)
    x_der(i, 1) = v_wDiscrepancy(i) - v_ideal(i);
    x_der(i, 2) = a_wDiscrepancy(i) - a_ideal(i);
end

Xi = (Theta'*Theta)^-1 * Theta' * x_der(:, 2);

correction = Theta * Xi;
a_DiscrepancySINDy = a_ideal(:) + correction;
v_DiscrepancySINDy = 0;
s_DiscrepancySINDy = 0;
for i = 1 : length(t) - 1
    v_DiscrepancySINDy(i+1) = v_DiscrepancySINDy(i) + a_DiscrepancySINDy(i) * dt;
    s_DiscrepancySINDy(i+1) = s_DiscrepancySINDy(i) + v_DiscrepancySINDy(i) * dt;
end

%% Plot
clf;
subplot(3,1,1);
title("Acceleration [m/s^2]")
hold on;
grid on;
grid minor;
plot(t, a_DiscrepancySINDy,'b-', t, a_wDiscrepancy,'r--', t, a_ideal, 'k-');
legend("Discrepancy SINDy", "Meas", "Ideal");
subplot(3,1,2);
title("Velocity [m/s]")
hold on;
grid on;
grid minor;
plot(t, v_DiscrepancySINDy,'b-', t, v_wDiscrepancy,'r--', t, v_ideal, 'k-');
legend("Discrepancy SINDy", "Meas", "Ideal");
subplot(3,1,3);
title("Distance [m]")
hold on;
grid on;
grid minor;
plot(t, s_DiscrepancySINDy,'b-', t, s_wDiscrepancy,'r--', t, s_ideal, 'k-');
legend("Discrepancy SINDy", "Meas", "Ideal");
