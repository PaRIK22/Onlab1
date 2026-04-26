%% 1. FIZIKAI SZIMULÁCIÓ ÉS KAMERA (15 fps)
clear; clc; close all;

% Paraméterek
m = 0.5; g = 9.81; rho = 1.225; A = 0.05; Cd = 0.47;
k_valos = (Cd * rho * A) / (2 * m); 

% --- Kamera beállítások ---
fps = 15; % 15 fps esetén borzasztó
dt = 1/fps;
t_cam = (0:dt:2)'; % 2 másodperces videofelvétel (oszlopvektor)

% Valós (zajmentes) adatok generálása a kamera időpillanataiban
ode_func = @(t, x) [x(2); g - k_valos * x(2)^2];
[~, x_data] = ode45(ode_func, t_cam, [0; 0]);
s_valos = x_data(:, 1);

% --- Zaj hozzáadása (Kamera és képelemzés pontatlansága) ---
% Tegyük fel, hogy a képpontfelbontás miatt átlagosan 1 cm (0.01 m) a hiba
zaj_szint = 0.01; 
s_zajos = s_valos + zaj_szint * randn(size(s_valos)); % Gauss-zaj

%% 2. SEBESSÉG ÉS GYORSULÁS SZÁMÍTÁSA NUMERIKUSAN
% Itt jön a nehézség: a zajos pozícióból kell deriválnunk!
% A MATLAB 'gradient' függvénye központi differenciát használ
v_szamitott = gradient(s_zajos, dt);
a_szamitott = gradient(v_szamitott, dt);

% Ábrázolás, hogy lássuk a zaj drasztikus hatását
figure('Name', 'Kamera Adatok (15 fps + Zaj)', 'Position', [100, 100, 1000, 300]);
subplot(1,3,1); plot(t_cam, s_zajos, 'b.-', t_cam, s_valos, 'k--'); 
title('Út (s)'); legend('Mért (zajos)', 'Valós'); xlabel('Idő [s]'); grid on;

subplot(1,3,2); plot(t_cam, v_szamitott, 'r.-'); 
title('Számított Sebesség (v)'); xlabel('Idő [s]'); grid on;

subplot(1,3,3); plot(t_cam, a_szamitott, 'g.-'); 
title('Számított Gyorsulás (a)'); xlabel('Idő [s]'); grid on;

%% 3. SINDy ALGORITMUS A ZAJOS ADATOKON
% A hibrid modell: a_maradek = a_szamitott - g
a_ismert = g .* ones(size(t_cam));
a_maradek = a_szamitott - a_ismert;

% Könyvtár építése a ZAJOS sebességből
Theta = [ones(size(v_szamitott)), v_szamitott, v_szamitott.^2, v_szamitott.^3]; 
konyvtar_nevek = {'1', 'v', 'v^2', 'v^3'};

% SINDy paraméterek
lambda = 0.01; 
Xi = Theta \ a_maradek; 

for k = 1:10
    kis_indexek = (abs(Xi) < lambda);
    Xi(kis_indexek) = 0;
    nagy_indexek = ~kis_indexek;
    if any(nagy_indexek)
        Xi(nagy_indexek) = Theta(:, nagy_indexek) \ a_maradek;
    end
end

%% Eredmények kiírása
clc;
fprintf('--- FIZIKA ---\n');
fprintf('A valós közegellenállási együttható (amit keresünk): k = %.4f\n\n', -k_valos);

fprintf('--- SINDy EREDMÉNYEK ZAJOS ADATOKBÓL ---\n');
for i = 1:length(Xi)
    fprintf('%5s tag együtthatója: %8.4f\n', konyvtar_nevek{i}, Xi(i));
end