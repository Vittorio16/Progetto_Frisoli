clear; clc; close all;

% --- 1. PARAMETRI DEL SISTEMA ---
m = [100; 100; 50; 50; 50; 50; 50; 50];           % Masse dei link [kg]
L = [3; 2; 2; 2; 2; 2; 2; 2];            % Lunghezze dei link [m]
g = 9.81;                % Gravità [m/s^2]

% Parametri base oscillante: x0(t) = A * sin(w * t)
A = 2;                   
w = 2 * pi * 0.1;        

n = length(m);
q0 = deg2rad([90; 90; 90; 90; 90; 90; 90; 90]);   % Angoli iniziali
dq0 = zeros(n, 1);       % Velocità iniziali
initial_state = [q0; dq0];

% Pre-calcolo della "massa cumulativa" 
mu = flipud(cumsum(flipud(m(:))));

% --- 2. INTEGRAZIONE NUMERICA ---
disp(['Avvio simulazione numerica per ', num2str(n), '-Pendolo...']);
t_span = [0, 20];

options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

tic;
[t_out, state_out] = ode113(@(t, state) n_pendulum_ode(t, state, mu, L, g, A, w, n), ...
                            linspace(t_span(1), t_span(2), 5000), ...
                            initial_state, options);
toc;

% --- 3. CALCOLO ENERGIA E NORMALIZZAZIONE ---
disp('Calcolo delle energie del sistema...');
theta_out = state_out(:, 1:n)';
omega_out = state_out(:, n+1:end)';

T_tot = zeros(1, length(t_out));
V_tot = zeros(1, length(t_out));
E_tot = zeros(1, length(t_out));

for k = 1:length(t_out)
    t = t_out(k);
    q = theta_out(:, k);
    dq = omega_out(:, k);
    
    v_x0 = A * w * cos(w * t);
    
    vx = zeros(n, 1);
    vy = zeros(n, 1);
    y  = zeros(n, 1);
    
    vx(1) = v_x0 + L(1) * dq(1) * cos(q(1));
    vy(1) = L(1) * dq(1) * sin(q(1));
    y(1)  = -L(1) * cos(q(1));
    
    for i = 2:n
        vx(i) = vx(i-1) + L(i) * dq(i) * cos(q(i));
        vy(i) = vy(i-1) + L(i) * dq(i) * sin(q(i));
        y(i)  = y(i-1) - L(i) * cos(q(i));
    end
    
    T_tot(k) = 0.5 * sum(m .* (vx.^2 + vy.^2));
    V_tot(k) = sum(m .* g .* y);
    E_tot(k) = T_tot(k) + V_tot(k);
end

% --- NORMALIZZAZIONE ---
% Calcolo dell'Energia Potenziale Caratteristica (Pendolo orizzontale)
E_rif = g * sum(mu .* L);

T_norm = T_tot / E_rif;
V_norm = V_tot / E_rif;
E_norm = E_tot / E_rif;
E_initial_norm = E_norm(1);

% Plot Energia Normalizzata
figure('Name', 'Energia Normalizzata', 'Color', 'w');
hold on; grid on;
plot(t_out, T_norm, 'b', 'LineWidth', 1.2, 'DisplayName', 'T / E_{rif}');
plot(t_out, V_norm, 'r', 'LineWidth', 1.2, 'DisplayName', 'V / E_{rif}');
plot(t_out, E_norm, 'k', 'LineWidth', 2, 'DisplayName', 'E_{tot} / E_{rif}');
yline(E_initial_norm, 'k--', 'LineWidth', 1.5, 'DisplayName', 'E_0 / E_{rif}');

xlabel('Tempo [s]', 'FontWeight', 'bold');
ylabel('Energia Adimensionale [-]', 'FontWeight', 'bold');
title('Fluttuazione dell''Energia Normalizzata');
legend('Location', 'best');
xlim([0, t_span(2)]);

% --- 4. POST-PROCESSING E ANIMAZIONE ---
disp('Avvio animazione...');

x_anim = zeros(n + 1, length(t_out));
y_anim = zeros(n + 1, length(t_out));

x_anim(1, :) = A * sin(w * t_out');
y_anim(1, :) = 0;

for i = 1:n
    x_anim(i+1, :) = x_anim(i, :) + L(i) * sin(theta_out(i, :));
    y_anim(i+1, :) = y_anim(i, :) - L(i) * cos(theta_out(i, :));
end

figure('Name', 'Animazione n-Pendolo su Base Oscillante', 'Color', 'w');
ax = axes; hold on; axis equal; grid on;
lim_x = sum(L) + A + 0.5;
lim_y = sum(L) + 0.5;
xlim(ax, [-lim_x, lim_x]);
ylim(ax, [-lim_y, 1]);
xlabel('X [m]'); ylabel('Y [m]');
title(sprintf('N-Pendolo Dinamico (n=%d)', n));

plot([-A-0.5, A+0.5], [0, 0], 'k--', 'LineWidth', 1);
h_line = plot(ax, x_anim(:, 1), y_anim(:, 1), 'o-k', 'LineWidth', 2.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'r');

disp('Premi un tasto qualsiasi per avviare l''animazione...');
pause; % Attende che tu sia pronto a guardare

% Loop animazione sincronizzato
tic; % Avvia il cronometro reale
for k = 1:15:length(t_out)
    set(h_line, 'XData', x_anim(:, k), 'YData', y_anim(:, k));
    drawnow; % Aggiorna la grafica forzatamente
    
    % Sincronizza il tempo di simulazione con il tempo reale
    tempo_simulato = t_out(k);
    while toc < tempo_simulato
        pause(0.002); % Pausa leggerissima per non bloccare la CPU
    end
end
disp('Processo terminato.');

% =========================================================================
% MOTORE FISICO MATRICIALE 
% =========================================================================
function dstate = n_pendulum_ode(t, state, mu, L, g, A, w, n)
    q = state(1:n);
    dq = state(n+1:end);
    
    ddx0 = -A * w^2 * sin(w * t);
    
    M = zeros(n, n); 
    C = zeros(n, 1); 
    G = zeros(n, 1); 
    F = zeros(n, 1); 
    
    for i = 1:n
        G(i) = mu(i) * g * L(i) * sin(q(i));
        F(i) = -mu(i) * L(i) * cos(q(i)) * ddx0;
        
        for j = 1:n
            idx = max(i, j);
            M(i, j) = mu(idx) * L(i) * L(j) * cos(q(i) - q(j));
            C(i) = C(i) + mu(idx) * L(i) * L(j) * sin(q(i) - q(j)) * dq(j)^2;
        end
    end
    
    ddq = M \ (F - C - G);
    dstate = [dq; ddq];
end