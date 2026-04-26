clear; clc; close all;

mu = 1.5; 
dt = 0.01;
t_meres = (0:dt:20)';

%x(1) a feszültség, x(2) a sebesség
ode_func = @(t, x) [x(2); (1.5 - 0.05*t) * (1 - x(1)^2) * x(2) - x(1)];
[~, x_adat] = ode45(ode_func, t_meres, [0.1; 0]); 
x_mert = x_adat(:, 1); 

zaj_szint = 0.001; 
x_zajos = x_mert + zaj_szint * randn(size(x_mert));

ablak_meret = 20; 
x_szurt = smoothdata(x_zajos, 'sgolay', ablak_meret);

v_szamitott = diff(x_szurt) / dt; 
a_szamitott = diff(v_szamitott) / dt;

x_bemenet = x_szurt(1 : end-2);
v_bemenet = v_szamitott(1 : end-1);
t_bemenet = t_meres(1 : end-2);

bemenetek = [x_bemenet'; v_bemenet'; t_bemenet'];
celok = a_szamitott';

% Konvertálás 'dlarray' (Deep Learning Array) formátumba a gradiens
% számításhoz -> diff helyett autodiff
X_dl = dlarray(bemenetek, 'CB'); %dlarray -> háttérben eltárolja a korábbi értékeket is
T_dl = dlarray(celok, 'CB'); %autodiff miatt a dlgradient  vissza tudja számolni melyik változás mennyivel csökkenti a hibát

layers = [
    featureInputLayer(3, 'Normalization', 'none') % 3 bemenet: x, v, t
    fullyConnectedLayer(30)                       % 1. rejtett réteg
    tanhLayer                                     % Aktivációs függvény
    fullyConnectedLayer(30)                       % 2. rejtett réteg
    tanhLayer                                     % Aktivációs függvény
    fullyConnectedLayer(1)                        % Kimenet: a (gyorsulás)
];
net = dlnetwork(layers);

numEpochs = 5000;
learnRate = 0.01;
velocity = []; % SGDM optimalizáló belső memóriája

for epoch = 1:numEpochs
    % A hiba és a gradiensek kiszámolása
    [loss, gradients] = dlfeval(@pinnLoss, net, X_dl, T_dl); %X_dl, T_dl a gradiens számításhoz
    
    % A hálózat súlyainak frissítése a hibák alapján (SGDM algoritmus)
    [net, velocity] = sgdmupdate(net, gradients, velocity, learnRate, 0.9); %90%-ban megtatja az előző sebesságet és irányt -> nem ragad be lok. minimumokba
    
    % if mod(epoch, 100) == 0
    %     fprintf('Epoch %d/%d - Teljes Hiba: %.4f\n', epoch, numEpochs, extractdata(loss));
    % end
end

a_nn_predikcio = extractdata(predict(net, X_dl)); %dlarrayról visszatérés, predikció kész hálóval
a_nn_plot = a_nn_predikcio'; 
% save('D:\Egyetem\!Msc\1.Félév\Önlab\Matlab\VanDerPol\PINN\VanDerPol_PINN_network.mat', 'net');

x_tiszta = x_adat(1:end-2, 1);
v_tiszta = x_adat(1:end-2, 2);
a_tiszta = (1.5 - 0.05*t_bemenet) .* (1 - x_tiszta.^2) .* v_tiszta - x_tiszta;

figure('Name', 'Van der Pol - Physics-Informed Neural Network', 'Position', [100, 100, 800, 400]);

plot(t_bemenet, a_szamitott, 'k', 'LineWidth', 0.5); hold on;
plot(t_bemenet, a_tiszta, 'b', 'LineWidth', 2);
plot(t_bemenet, a_nn_plot, 'r--', 'LineWidth', 2);

title('Gyorsulás: Zajos Érték vs. Várt Érték vs. PINN');
xlabel('Idő [s]'); ylabel('Gyorsulás [v^2])');
legend('Zajos jel', 'Várt kimenet', 'PINN');
grid on;

function [loss, gradients] = pinnLoss(net, X, T)
    % Háló predikció
    Y_pred = forward(net, X);
    %y_pred_diff = gradient(y_pred, x)
    %y_pred_diff_diff = gradient(y_pred_diff, x)
    
    % 2. Adat-Hiba -> Mekkora az eltérés
    loss_data = mse(Y_pred, T);
    
    % Fizikai paraméterek
    x = X(1,:);
    v = X(2,:);
    t = X(3,:);
    
    % Számított érték
    a_fizika = (1.5 - 0.05*t) .* (1 - x.^2) .* v - x;
    
    % Büntetés, ha eltér a predikció a valós fizikától
    loss_physics = mse(Y_pred, a_fizika);
    
    % Adat + Fizika hiba
    lambda_phys = 1; % 0 = Sima háló - Nagyobb érték -> Mekkora hatása legyen az egyenletnek
    loss = loss_data * 0 + lambda_phys * loss_physics;
    
    % Gradiens backpropagationhöz
    gradients = dlgradient(loss, net.Learnables); %Minden neuron súlyához meredekség -> dlnetwork miatt
end