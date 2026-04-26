clear; clc; close all;

mu = 1.5; 
dt = 0.01;
t_meres = (0:dt:20)';

ode_func = @(t, x) [x(2); (1.5 - 0.05*t) * (1 - x(1)^2) * x(2) - x(1)];
[~, x_adat] = ode45(ode_func, t_meres, [0.1; 0]); 

x_mert = x_adat(:, 1); 

% zaj_szint = 0.01; 
% x_mert = x_mert + zaj_szint * randn(size(x_mert));
% ablak_meret = 50; 
% x_mert = smoothdata(x_mert, 'sgolay', ablak_meret);

v_szamitott = diff(x_mert) / dt; 
a_szamitott = diff(v_szamitott) / dt;

x_konyvtar = x_mert(1 : end-2);
t_konyvtar = t_meres(1 : end-2);
v_konyvtar = v_szamitott(1 : end-1);

Theta = [ones(size(x_konyvtar)), x_konyvtar, v_konyvtar, t_konyvtar, t_konyvtar.*v_konyvtar, ...
         x_konyvtar.*v_konyvtar, t_konyvtar.*x_konyvtar, x_konyvtar.^2, v_konyvtar.^2,...
         (x_konyvtar.^2).*v_konyvtar, x_konyvtar.*(v_konyvtar.^2), ...
         t_konyvtar.*(x_konyvtar.^2).*v_konyvtar, ... 
         x_konyvtar.^3, v_konyvtar.^3]; 

konyvtar_nevek = {'1', 'x', 'v', 't', 't*v', 'x*v', 't*x', 'x^2', ...
                  'v^2', 'x^2*v', 'x*v^2', 't*x^2*v', 'x^3', 'v^3'};

lambda = 0.03;
Xi = Theta \ a_szamitott; 

for k = 1:10
    kis_indexek = (abs(Xi) < lambda);
    Xi(kis_indexek) = 0;              
    nagy_indexek = ~kis_indexek;      
    
    if any(nagy_indexek)
        Xi(nagy_indexek) = Theta(:, nagy_indexek) \ a_szamitott;
    end
end

clc;
fprintf('--- Van der Pol Oszcillátor (1%%-os zajjal) ---\n\n');
fprintf('A keresett (elméleti) együtthatók:\n');
fprintf('x     : -1.0000\n');
fprintf('v     :  1.5000\n');
fprintf('x^2*v : -1.5000\n\n');

fprintf('--- SINDy Által Talált Egyenlet ---\n');
for i = 1:length(Xi)
    if Xi(i) ~= 0
        fprintf('%6s tag együtthatója: %8.4f\n', konyvtar_nevek{i}, Xi(i));
    end
end

a_sindy = Theta * Xi;
t_plot = t_meres(1:end-2);

figure('Name', 'Van der Pol Oszcillátor - SINDy Zajos Eredmények', 'Position', [100, 100, 1200, 400]);

subplot(1, 3, 1);
plot(t_plot, x_konyvtar, 'b', 'LineWidth', 1.5); 
hold on;
plot(t_plot, v_konyvtar, 'r', 'LineWidth', 1.5);
title('1. Időtartomány (Jelek)');
xlabel('Idő [s]'); ylabel('Amplitúdó');
legend('Feszültség (x)', 'Sebesség (v)');
grid on;

subplot(1, 3, 2);
plot(x_konyvtar, v_konyvtar, 'k', 'LineWidth', 1.5); hold on;
plot(x_konyvtar(1), v_konyvtar(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
title('2. Fázissík (Határciklus)');
xlabel('Feszültség (x)'); ylabel('Sebesség (v)');
legend('Rendszer pályája', 'Kezdőpont', 'Location', 'southeast');
grid on;

subplot(1, 3, 3);
plot(t_plot, a_szamitott, 'k', 'LineWidth', 2); hold on;
plot(t_plot, a_sindy, 'r--', 'LineWidth', 1.5);
title('3. Modell Teszt (Gyorsulás)');
xlabel('Idő [s]'); ylabel('Gyorsulás (a)');
legend('Számított (diff)', 'SINDy Modell (Theta * Xi)');
grid on;