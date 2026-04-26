%% 1. FIZIKAI SZIMULÁCIÓ (Valós adatgenerálás)
clear; clc; close all;

% Paraméterek (SI mértékegységekben)
m = 0.9;        % Labda tömege (kg)
g = 9.81;       % Gravitációs gyorsulás (m/s^2)
rho = 1.225;    % Levegő sűrűsége (kg/m^3)
A = 0.5;       % Keresztmetszet (m^2)
Cd = 0.47;      % Közegellenállási tényező (Gömb)

% Közegellenállás konstansa: F_d = (1/2) * C_d* rho * A * v^2 = k * v^2 
% m * a = m * g - k * v^2
% a = g - (Cd * rho * A / 2m) * v^2
k_valos = (Cd * rho * A) / (2 * m); 

fprintf('A valós közegellenállási együttható: k = %.4f\n', k_valos);

% Differenciálegyenlet megoldása (ode45)
% Állapotvektor: x = [s; v] (út, sebesség)
% A pozitív irányt lefelé definiáljuk!
t_span = [0 20]; % 20 másodperc szimuláció
x0 = [0; 0];    % Kezdeti út 0, kezdeti sebesség 0

% ODE függvény definiálása: dx/dt = [v; g - k*v^2]
ode_func = @(t, x) [x(2); g - k_valos * x(2)^2];

[t_data, x_data] = ode45(ode_func, t_span, x0);

s_data = x_data(:, 1); % Út (magasság változása lefelé)
v_data = x_data(:, 2); % Sebesség

% Gyorsulás kiszámítása a valós adatokból
a_valos = g - k_valos * v_data.^2;

%% Eredmények ábrázolása
figure('Name', 'Valós Fizikai Szimuláció', 'Position', [100, 100, 1000, 300]);
subplot(1,3,1); plot(t_data, s_data, 'b', 'LineWidth', 2); title('Út (s)'); xlabel('Idő [s]'); ylabel('[m]'); grid on;
subplot(1,3,2); plot(t_data, v_data, 'r', 'LineWidth', 2); title('Sebesség (v)'); xlabel('Idő [s]'); ylabel('[m/s]'); grid on;
subplot(1,3,3); plot(t_data, a_valos, 'k', 'LineWidth', 2); title('Gyorsulás (a)'); xlabel('Idő [s]'); ylabel('[m/s^2]'); grid on;


%% 2. SINDy ALGORITMUS TESZTELÉSE
%A teljes rendszer = M1 (ismert) + NN/SINDy (ismeretlen)
% Kiszámoljuk azt a gyorsulást, amit a hiányos modellünk meg tud magyarázni:
a_ismert = g .* ones(size(t_data)); 

% A maradék, amit a SINDy-nek kell megmagyaráznia (ez a Fd/m hatása):
a_maradek = a_valos - a_ismert; 

% --- SINDy Könyvtár (Library) építése ---
% Létrehozunk egy Theta mátrixot a jelöltekből
Theta = [ones(size(v_data)), v_data, v_data.^2, v_data.^3]; 
konyvtar_nevek = {'1', 'v', 'v^2', 'v^3'};

% --- SINDy Megoldó: Sequential Thresholded Least Squares (STLSQ) ---
% Kezdeti becslés p_hat = (A^T * A)^-1 * A^T * y
% Matlabban \ pont ezt csinálja
Xi = Theta \ a_maradek; 

% SINDy Sparsitás (Ritkítás) paramétere
lambda = 0.01; % Ez alatti együtthatókat nullának tekintjük -> ha túl nagy kiszűri a légellenállást
iteraciok = 10;

for k = 1:iteraciok
    % Megkeressük a kis együtthatókat
    kis_indexek = (abs(Xi) < lambda);
    
    % Nullázzuk a kis együtthatókat (ez a SINDy lényege, kiszűri a zajt/felesleget)
    Xi(kis_indexek) = 0;
    
    % Újraszámoljuk a legkisebb négyzeteket csak a megmaradt elemekre
    nagy_indexek = ~kis_indexek;
    
    if any(nagy_indexek) % Ha maradt egyáltalán valami
        % Csak a kiválasztott oszlopokkal számolunk
        Xi(nagy_indexek) = Theta(:, nagy_indexek) \ a_maradek;
    end
end

%% Eredmények kiírása
fprintf('\n--- SINDy Eredmények ---\n');
fprintf('A SINDy által talált egyenlet a maradék gyorsulásra:\n a_maradek = ');

talalt_tagok = 0;
for i = 1:length(Xi)
    if Xi(i) ~= 0
        if talalt_tagok > 0
            fprintf(' + ');
        end
        fprintf('(%.4f) * %s', Xi(i), konyvtar_nevek{i});
        talalt_tagok = talalt_tagok + 1;
    end
end

if talalt_tagok == 0
    fprintf('0 (Nem talált összefüggést)');
end
fprintf('\n\n');
fprintf('Összehasonlítás:\n');
fprintf('Elvárt együttható a v^2 taghoz: %.4f\n', -k_valos);
fprintf('Talált együttható a v^2 taghoz: %.4f\n', Xi(3));