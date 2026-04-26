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

% Csak EGYSZER deriválunk 
v_szamitott = diff(x_szurt) / dt; 

% Mivel csak egyszer deriváltunk, a vektorok 1-gyel rövidebbek
x_konyvtar = x_szurt(1 : end-1);
t_konyvtar = t_meres(1 : end-1);

v_cel = v_szamitott - v_szamitott(1);

Theta = [ones(size(x_konyvtar)), x_konyvtar, v_szamitott, t_konyvtar, t_konyvtar.*v_szamitott, ...
         x_konyvtar.*v_szamitott, t_konyvtar.*x_konyvtar, x_konyvtar.^2, v_szamitott.^2,...
         (x_konyvtar.^2).*v_szamitott, x_konyvtar.*(v_szamitott.^2), ...
         t_konyvtar.*(x_konyvtar.^2).*v_szamitott, ... 
         x_konyvtar.^3, v_szamitott.^3]; 

konyvtar_nevek = {'1', 'x', 'v', 't', 't*v', 'x*v', 't*x', 'x^2', ...
                  'v^2', 'x^2*v', 'x*v^2', 't*x^2*v', 'x^3', 'v^3'};
%Integrált könyvtár
Theta_int = cumtrapz(t_konyvtar, Theta);

lambda = 0.02;
Xi = Theta_int \ v_cel; 

for k = 1:10
    kis_indexek = (abs(Xi) < lambda);
    Xi(kis_indexek) = 0;              
    nagy_indexek = ~kis_indexek;      
    
    if any(nagy_indexek)
        Xi(nagy_indexek) = Theta_int(:, nagy_indexek) \ v_cel;
    end
end

clc;
fprintf('--- Van der Pol Oszcillátor (Integrál SINDy 1%% zajjal) ---\n\n');
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

% --- PLOT ---
a_sindy = Theta * Xi;

%Csak a plot miatt
a_zajos_plot = diff(v_szamitott) / dt;
t_plot_a = t_konyvtar(1:end-1);

figure('Name', 'Van der Pol - Integrál SINDy (Zajos)', 'Position', [100, 100, 1200, 400]);

subplot(1, 3, 1);
plot(t_konyvtar, x_konyvtar, 'b', 'LineWidth', 1.5); hold on;
plot(t_konyvtar, v_szamitott, 'r', 'LineWidth', 1.5);
title('1. Időtartomány');
xlabel('Idő [s]'); ylabel('Amplitúdó');
legend('Feszültség (x)', 'Sebesség (v)');
grid on;

subplot(1, 3, 2);
plot(x_konyvtar, v_szamitott, 'k', 'LineWidth', 1.5); hold on;
plot(x_konyvtar(1), v_szamitott(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
title('2. Fázissík');
xlabel('Feszültség (x)'); ylabel('Sebesség (v)');
legend('Rendszer pályája', 'Kezdőpont', 'Location', 'southeast');
grid on;

subplot(1, 3, 3);
plot(t_plot_a, a_zajos_plot, 'k', 'LineWidth', 0.5); hold on;
plot(t_konyvtar, a_sindy, 'r--', 'LineWidth', 2);
title('3. Gyorsulás');
xlabel('Idő [s]'); ylabel('Gyorsulás (a)');
legend('Zajos diff(v)', 'Integrál SINDy');
grid on;