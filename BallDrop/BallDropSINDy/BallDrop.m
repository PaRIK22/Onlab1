clear;
clf;
%% Constants
g = 9.81; % [m/s^2]
m = 0.045359; % [kg]
dt = 1/15; % [s]
t = 0:dt:30; % [s]
v_ideal(1) = 0; % [m/s]
s_ideal(1) = 0; % [m]
Rho = 1.225; % [kg/m^3]
mu = 1.81e-5; % [Pa*s]
D = 0.021963 * 2; % [m] 
A = D^2 * pi() / 4; % [m^2]
Re(1) = 0;
Cd(1) = 0;
a_ideal(1) = (m * g) / m;

%% Equations
% F = m * a
% - Fd + (m * g) = m * a
% a = -Fd / m + (m * g) / m

% Fd = (1/2) * Rho * v^2 * Cd(Re) * A
% Re = (Rho * v * D) / mu

% a = -(1/2 * Rho * v^2 * Cd(Re) * A) / m + (m * g) / m

for i = 1 : length(t) - 1
    v_ideal(i+1) = v_ideal(i) + a_ideal(i) * dt;
    s_ideal(i+1) = s_ideal(i) + v_ideal(i) * dt;
    
    Re(i+1) = (Rho * v_ideal(i+1) * D) / mu;
    Cd(i+1) = 24/Re(i+1) ...
    + (2.6.*(Re(i+1)/5.0)) / (1 + (Re(i+1)/5.0).^1.52) ...
    + (0.411.*(Re(i+1)/(2.63e5))^(-7.94)) / (1 + (Re(i+1)/(2.63e5))^(-8.00)) ...
    + (0.25.*(Re(i+1)/10e6)) / (1 + (Re(i+1)/10e6));
    
    a_ideal(i+1) = -(1/2 * Rho * v_ideal(i+1)^2 * Cd(i+1) * A) / m + (m * g) / m;
end

%% SINDy

s_meas = s_ideal;

noise_level = 0;    
gauss = noise_level * randn(size(s_meas));

s_meas = s_meas + gauss;

v_meas = diff(s_meas);
v_meas(length(v_meas) + 1) = v_meas(length(v_meas));
v_meas = v_meas / dt;

a_meas = diff(v_meas);
a_meas(length(a_meas) + 1) = a_meas(length(a_meas));
a_meas = a_meas / dt;


for i = 1 : length(a_ideal)
    Theta(i,1:10) = [1, v_meas(i), s_meas(i), v_meas(i)*s_meas(i), v_meas(i)^2, ...
    s_meas(i)^2, s_meas(i)^2*v_meas(i), s_meas(i)*v_meas(i)^2, s_meas(i)^3, v_meas(i)^3];
end

for i = 1 : length(s_ideal)
    x(i, 1) = s_meas(i);
    x(i, 2) = v_meas(i);
end

for i = 1 : length(v_ideal)
    x_der(i, 1) = v_meas(i);
    x_der(i, 2) = a_meas(i);
end

Xi= inv(Theta'*Theta) * Theta' * x_der(:, 2);

a_IdealSINDy = Theta * Xi;
v_IdealSINDy = 0;
s_IdealSINDy = 0;
for i = 1 : length(t) - 1
    v_IdealSINDy(i+1) = v_IdealSINDy(i) + a_IdealSINDy(i) * dt;
    s_IdealSINDy(i+1) = s_IdealSINDy(i) + v_IdealSINDy(i) * dt;
end

%% Plot
subplot(3,1,1);
plot(t, a_IdealSINDy, 'b-', t, a_ideal, 'r--');
title("Acceleration [m/s^2]")
grid on;
grid minor;
legend("SINDy", "Ideal");
subplot(3,1,2);
plot(t, v_IdealSINDy, 'b-', t, v_ideal, 'r--');
title("Velocity [m/s]")
grid on;
grid minor;
legend("SINDy", "Ideal");
subplot(3,1,3);
plot(t, s_IdealSINDy, 'b-', t, s_ideal, 'r--');
title("Distance [m]")
grid on;
grid minor;
legend("SINDy", "Ideal");

