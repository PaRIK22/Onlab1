clear;
%% Constants
g = 9.81; % [m/s^2]
m = 0.045359; % [kg]
dt = 1/15; % [s]
t = 0:dt:30; % [s]
v(1) = 0; % [m/s]
s(1) = 0; % [m]
Rho = 1.225; % [kg/m^3]
mu = 1.81e-5; % [Pa*s]
D = 0.021963 * 2; % [m] 
A = D^2 * pi() / 4; % [m^2]
Re(1) = 0;
Cd(1) = 0;

%% Equations
% F = m * a
% - Fd + (m * g) = m * a
% a = -Fd / m + (m * g) / m

% Fd = (1/2) * Rho * v^2 * Cd(Re) * A
% Re = (Rho * v * D) / mu

% a = -(1/2 * Rho * v^2 * Cd(Re) * A) / m + (m * g) / m

a(1) = (m * g) / m;


for i = 1 : length(t) - 1
    v(i+1) = v(i) + a(i) * dt;
    s(i+1) = s(i) + v(i) * dt;
    
    Re(i+1) = (Rho * v(i+1) * D) / mu;
    Cd(i+1) = 24/Re(i+1) ...
    + (2.6.*(Re(i+1)/5.0)) / (1 + (Re(i+1)/5.0).^1.52) ...
    + (0.411.*(Re(i+1)/(2.63e5))^(-7.94)) / (1 + (Re(i+1)/(2.63e5))^(-8.00)) ...
    + (0.25.*(Re(i+1)/10e6)) / (1 + (Re(i+1)/10e6));
    a(i+1) = -(1/2 * Rho * v(i+1)^2 * Cd(i+1) * A) / m + (m * g) / m;
end

%% SINDy

s_sim = s;

noise_level = 0.1;    
gauss = noise_level * randn(size(s_sim));

s_sim = s_sim + gauss;

v_sim = gradient(s_sim,t);
a_sim = gradient(v_sim,t);


for i = 1 : length(a)
    Theta(i,1:10) = [1, v_sim(i), s_sim(i), v_sim(i)*s_sim(i), v_sim(i)^2,...
    s_sim(i)^2, s_sim(i)^2*v_sim(i), s_sim(i)*v_sim(i)^2, s_sim(i)^3, v_sim(i)^3];
end

for i = 1 : length(s)
    x(i, 1) = s_sim(i);
    x(i, 2) = v_sim(i);
end

for i = 1 : length(v)
    x_der(i, 1) = v_sim(i);
    x_der(i, 2) = a_sim(i);
end

Xi= (Theta'*Theta)^-1 * Theta' * x_der(:, 2);

a_sindy = Theta * Xi;
v_sindy = 0;
s_sindy = 0;
for i = 1 : length(t) - 1
    v_sindy(i+1) = v_sindy(i) + a_sindy(i) * dt;
    s_sindy(i+1) = s_sindy(i) + v_sindy(i) * dt;
end

%% Plot
subplot(3,1,1);
title("Acceleration [m/s^2]")
hold on;
grid on;
grid minor;
plot(t, a_sindy, t, a);
legend("SINDy", "Ideal");
subplot(3,1,2);
title("Velocity [m/s]")
hold on;
grid on;
grid minor;
plot(t, v_sindy, t, v);
legend("SINDy", "Ideal");
subplot(3,1,3);
title("Distance [m]")
hold on;
grid on;
grid minor;
plot(t, s_sindy, t, s);
legend("SINDy", "Ideal");

