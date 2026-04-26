%% DINAMIKUS C_D (REYNOLDS-FÜGGŐ) ÉS SINDy MODELL VISSZASZIMULÁLÁSA
clear; clc; close all;

% --- 1. Fizikai Paraméterek ---
m = 0.5; g = 9.81; rho = 1.225; D = 0.05; A = (pi * D^2) / 4; mu = 1.81e-5;

%White -féle
Re_func = @(v) (rho * max(abs(v), 1e-6) * D) / mu;
Cd_func = @(v) 24 ./ Re_func(v) + 6 ./ (1 + sqrt(Re_func(v))) + 0.4;

% Valós ODE
ode_func_valos = @(t, x) [x(2); g - (Cd_func(x(2)) * rho * A * x(2)^2) / (2 * m)];

% Fix időlépések az egyszerűbb összehasonlításért (5 mp, 100 lépés)
t_span = linspace(0, 5, 100)'; 
[~, x_valos] = ode45(ode_func_valos, t_span, [0; 0]);

s_valos = x_valos(:, 1);
v_valos = x_valos(:, 2);

% Valós gyorsulás kiszámítása
a_valos = zeros(size(v_valos));
for i = 1:length(v_valos)
    a_valos(i) = g - (Cd_func(v_valos(i)) * rho * A * v_valos(i)^2) / (2 * m);
end

%% --- 2. SINDy ALGORITMUS ---
a_ismert = g .* ones(size(t_span));
a_maradek = a_valos - a_ismert;

% Könyvtár építése
Theta = [ones(size(s_valos)), s_valos, v_valos, s_valos.*v_valos, ...
         s_valos.^2, v_valos.^2, (s_valos.^2).*v_valos, s_valos.*(v_valos.^2), ...
         s_valos.^3, v_valos.^3]; 

konyvtar_nevek = {'1', 's', 'v', 's*v', 's^2', 'v^2', 's^2*v', 's*v^2', 's^3', 'v^3'};

% SINDy optimalizálás
lambda = 0.0001; % Elég kicsi, hogy a SINDy próbálkozhasson a polinomokkal
Xi = Theta \ a_maradek; 
for k = 1:10
    kis_indexek = (abs(Xi) < lambda);
    Xi(kis_indexek) = 0;
    nagy_indexek = ~kis_indexek;
    if any(nagy_indexek)
        Xi(nagy_indexek) = Theta(:, nagy_indexek) \ a_maradek;
    end
end

%% --- 3. A SINDy MODELL VISSZASZIMULÁLÁSA (A VÁRT ÉS SZIMULÁLT ÖSSZEVETÉSE) ---
% Létrehozunk egy új ODE függvényt a SINDy által talált paraméterekkel
% a_sindy = g + Theta(pillanatnyi) * Xi
ode_func_sindy = @(t, x) [
    x(2); 
    g + [1, x(1), x(2), x(1)*x(2), x(1)^2, x(2)^2, (x(1)^2)*x(2), x(1)*(x(2)^2), x(1)^3, x(2)^3] * Xi
];

% Végigfuttatjuk a szimulációt a kitalált modellel is
[~, x_sindy] = ode45(ode_func_sindy, t_span, [0; 0]);
s_sindy = x_sindy(:, 1);
v_sindy = x_sindy(:, 2);

% A becsült gyorsulás kiszámítása
a_sindy = zeros(size(t_span));
for i = 1:length(t_span)
    % Az adott pillanatban érvényes állapotokból felépítjük a Theta sort
    Theta_aktualis = [1, s_sindy(i), v_sindy(i), s_sindy(i)*v_sindy(i), ...
                      s_sindy(i)^2, v_sindy(i)^2, (s_sindy(i)^2)*v_sindy(i), s_sindy(i)*(v_sindy(i)^2), ...
                      s_sindy(i)^3, v_sindy(i)^3];
    a_sindy(i) = g + Theta_aktualis * Xi;
end

%% --- 4. EREDMÉNYEK KIÍRÁSA ÉS ÁBRÁZOLÁSA ---
clc;
fprintf('--- SINDy Által Talált Egyenlet ---\n');
talalt = false;
for i = 1:length(Xi)
    if Xi(i) ~= 0
        fprintf('%6s tag együtthatója: %8.4f\n', konyvtar_nevek{i}, Xi(i));
        talalt = true;
    end
end
if ~talalt
    fprintf('A SINDy csak a konstans gravitációt hagyta meg.\n');
end

% Ábrázolás
figure('Name', 'Összehasonlítás: Valós (Gyökös C_D) vs. SINDy (Polinomos)', 'Position', [100, 100, 1000, 350]);

subplot(1,3,1);
plot(t_span, s_valos, 'k', 'LineWidth', 2); hold on;
plot(t_span, s_sindy, 'r--', 'LineWidth', 2);
title('Út (s)'); xlabel('Idő [s]'); ylabel('[m]'); 
legend('Valós (Bonyolult C_D)', 'SINDy (Becsült)'); grid on;

subplot(1,3,2);
plot(t_span, v_valos, 'k', 'LineWidth', 2); hold on;
plot(t_span, v_sindy, 'r--', 'LineWidth', 2);
title('Sebesség (v)'); xlabel('Idő [s]'); ylabel('[m/s]'); 
legend('Valós', 'SINDy'); grid on;

subplot(1,3,3);
plot(t_span, a_valos, 'k', 'LineWidth', 2); hold on;
plot(t_span, a_sindy, 'r--', 'LineWidth', 2);
title('Gyorsulás (a)'); xlabel('Idő [s]'); ylabel('[m/s^2]'); 
legend('Valós', 'SINDy'); grid on;