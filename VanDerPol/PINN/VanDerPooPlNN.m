%% VAN DER POL OSZCILLÁTOR - NEURÁLIS HÁLÓZAT (MLP)
clear; clc; close all;

mu = 1.5; 
dt = 0.01;
t_meres = (0:dt:20)';

ode_func = @(t, x) [x(2); (1.5 - 0.05*t) * (1 - x(1)^2) * x(2) - x(1)];
[~, x_adat] = ode45(ode_func, t_meres, [0.1; 0]); 
x_mert = x_adat(:, 1); 

zaj_szint = 0.01; 
x_zajos = x_mert + zaj_szint * randn(size(x_mert));
ablak_meret = 50; 
x_szurt = smoothdata(x_zajos, 'sgolay', ablak_meret);

% Kétszer deriválunk, hogy meglegyen a gyorsulás (a célpont)
v_szamitott = diff(x_szurt) / dt; 
a_szamitott = diff(v_szamitott) / dt;

x_bemenet = x_szurt(1 : end-2);
x_mert = x_mert(1 : end-2);
v_bemenet = v_szamitott(1 : end-1);
t_bemenet = t_meres(1 : end-2);

bemenetek = [x_bemenet'; v_bemenet'; t_bemenet'];

celok = a_szamitott';

%feedforward 2x10
net = fitnet([10 10]);

%Tanítási paraméterek
net.trainParam.showWindow = true;

net = train(net, bemenetek, celok);

a_nn_predikcio = net(bemenetek);
a_nn_plot = a_nn_predikcio'; 

% ==========================================================

clc;
fprintf('--- Van der Pol Oszcillátor (Neurális Hálóval) ---\n\n');

% --- PLOT ---
figure('Name', 'Van der Pol - Neurális Hálózat', 'Position', [100, 100, 1200, 400]);

subplot(1, 3, 1);
plot(t_bemenet, x_bemenet, 'b', 'LineWidth', 1.5); hold on;
plot(t_bemenet, v_bemenet, 'r', 'LineWidth', 1.5);
title('1. Időtartomány (Bemenetek)');
xlabel('Idő [s]'); ylabel('Amplitúdó');
legend('Feszültség (x)', 'Sebesség (v)');
grid on;

subplot(1, 3, 2);
plot(x_bemenet, v_bemenet, 'k', 'LineWidth', 1.5); hold on;
plot(x_bemenet(1), v_bemenet(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
title('2. Fázissík');
xlabel('Feszültség (x)'); ylabel('Sebesség (v)');
legend('Rendszer pályája', 'Kezdőpont', 'Location', 'southeast');
grid on;

subplot(1, 3, 3);
plot(t_bemenet, a_szamitott, 'k', 'LineWidth', 0.5); hold on;
plot(t_bemenet, a_nn_plot, 'r--', 'LineWidth', 2);
plot(t_bemenet, x_mert, 'g--', 'LineWidth', 2);
title('3. Gyorsulás: Valóság vs. Hálózat');
xlabel('Idő [s]'); ylabel('Gyorsulás (a)');
legend('Kiszámolt (Cél)', 'Neurális Háló Jóslata');
grid on;