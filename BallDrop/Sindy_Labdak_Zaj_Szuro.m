%% 3. SAVITZKY-GOLAY SZŰRŐ ÉS SINDy A ZAJOS ADATOKON
clear; clc; close all;

% --- 1. Fizikai Paraméterek és Valós Adatok ---
m = 0.5; g = 9.81; rho = 1.225; A = 0.05; Cd = 0.47;
k_valos = (Cd * rho * A) / (2 * m); 

fps = 15;
dt = 1/fps;
t_cam = (0:dt:2)'; % 2 másodperc videó

ode_func = @(t, x) [x(2); g - k_valos * x(2)^2];
[~, x_data] = ode45(ode_func, t_cam, [0; 0]);
s_valos = x_data(:, 1);

% --- 2. Kamera Zaj Hozzáadása ---
rng(42); % Fixáljuk a véletlenszám-generátort, hogy mindig ugyanazt a zajt kapjuk
zaj_szint = 0.01; % 1 cm átlagos hiba
s_zajos = s_valos + zaj_szint * randn(size(s_valos));

% --- 3. ADATSZŰRÉS: Savitzky-Golay ---
% A cikkben is ezt használták a zaj minimalizálására deriválás előtt
polinom_fok = 3;   % 3. fokú lokális polinomokat illesztünk
ablak_meret = 7;   % 7 adatpontos (képkockás) csúszóablakban (páratlan kell legyen!)

s_szurt = sgolayfilt(s_zajos, polinom_fok, ablak_meret);

% --- 4. Deriválás a SZŰRT adatokból ---
v_szurt = gradient(s_szurt, dt);
a_szurt = gradient(v_szurt, dt);

% Ábrázolás: Nézzük meg, mit csinált a szűrő!
figure('Name', 'Szűrt vs. Zajos Gyorsulás', 'Position', [100, 100, 800, 400]);
a_valos = g - k_valos * x_data(:,2).^2; % Az elméleti, tökéletes gyorsulás
a_zajos_nyers = gradient(gradient(s_zajos, dt), dt); % Szűrés nélküli, vad gyorsulás

plot(t_cam, a_valos, 'k', 'LineWidth', 2); hold on;
plot(t_cam, a_zajos_nyers, 'r:', 'LineWidth', 1);
plot(t_cam, a_szurt, 'b.-', 'LineWidth', 1.5);
title('A gyorsulás ($a$) kiszámítása'); 
legend('Valós fizika', 'Nyers derivált (Szűrés nélkül)', 'Szűrt derivált (Savitzky-Golay)');
xlabel('Idő [s]'); ylabel('Gyorsulás [m/s^2]'); grid on;

% --- 5. SINDy ALGORITMUS ---
a_ismert = g .* ones(size(t_cam));
a_maradek = a_szurt - a_ismert;

Theta = [ones(size(v_szurt)), v_szurt, v_szurt.^2, v_szurt.^3]; 
konyvtar_nevek = {'1', 'v', 'v^2', 'v^3'};

% SINDy hiperparaméter
lambda = 0.1; % Magasabb küszöböt állítunk be, hogy a maradék zajt kiszűrjük
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
fprintf('--- FIZIKA (Ideális) ---\n');
fprintf('A valós együttható a v^2 taghoz: %.4f\n\n', -k_valos);

fprintf('--- SINDy EREDMÉNYEK (Savitzky-Golay szűrt adatokból) ---\n');
talalt_tag = false;
for i = 1:length(Xi)
    if Xi(i) ~= 0
        fprintf('%5s tag együtthatója: %8.4f\n', konyvtar_nevek{i}, Xi(i));
        talalt_tag = true;
    end
end
if ~talalt_tag
    fprintf('A SINDy a zaj miatt minden tagot eldobott (Túl nagy lambda vagy túl kis jel).\n');
end