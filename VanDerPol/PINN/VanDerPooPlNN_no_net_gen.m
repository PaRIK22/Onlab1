clear; clc; close all;

load('D:\Egyetem\!Msc\1.Félév\Önlab\Matlab\VanDerPol\PINN\VanDerPol_PINN_network.mat', 'net');

% 2. Új bemeneti adatok generálása
dt = 0.01;
t_uj = (0:dt:20)'; 

% Adatok
ode_func = @(t, x) [x(2); (1.5 - 0.05*t) * (1 - x(1)^2) * x(2) - x(1)];
[~, x_adat] = ode45(ode_func, t_uj, [0.2; 0]); % Más kezdőfeltétel

x_uj = x_adat(:, 1);
v_uj = x_adat(:, 2);

% Adatok előkészítése 
bemenetek = [x_uj'; v_uj'; t_uj'];
X_dl = dlarray(bemenetek, 'CB');

% Pred
a_pred_dl = predict(net, X_dl);

% Kicsomagoljuk a számokat
a_nn_pred = extractdata(a_pred_dl)';

% 5. Eredmény megjelenítése
figure('Name', 'PINN Inferencia');
plot(t_uj, a_nn_pred, 'r', 'LineWidth', 2);
grid on;
title('A betöltött PINN hálózat jóslata új adatokra');
xlabel('Idő [s]');
ylabel('Gyorsulás [a]');