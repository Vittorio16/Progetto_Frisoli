clear; clc; close all;

% --- 1. PARAMETRI DEL SISTEMA ---
m = [5; 4; 3];           % Masse dei link [kg]
L = [3; 2; 1];           % Lunghezze dei link [m]
g = 9.81;                % Gravità [m/s^2]

% Parametri base oscillante: x0(t) = A * sin(w * t)
A = 2;                   
w = 5.0 * pi * 0.1;        

n = length(m);
q0 = deg2rad([90; 90; 90]);   % Angoli iniziali
dq0 = zeros(n, 1);       % Velocità iniziali
initial_state = [q0; dq0];

% massa cumulativa 
mu = flipud(cumsum(flipud(m(:))));

% --- 2. INTEGRAZIONE NUMERICA ---
disp(['Avvio simulazione numerica per ', num2str(n), '-Pendolo...']);
t_span = [0, 60];

options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

tic;
[t_out, state_out] = ode15s(@(t, state) n_pendulum_ode(t, state, mu, L, g, A, w, n), ...
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
P_base = zeros(1, length(t_out)); % Potenza

for k = 1:length(t_out)
    t = t_out(k);
    q = theta_out(:, k);
    dq = omega_out(:, k);
    
    % Velocità e accelerazione base
    v_x0 = A * w * cos(w * t);
    ddx0 = -A * w^2 * sin(w * t);
    
    % Richiama l'ODE per ottenere le accelerazioni angolari (ddq)
    state_k = state_out(k, :)';
    dstate = n_pendulum_ode(t, state_k, mu, L, g, A, w, n);
    ddq = dstate(n+1:end);
    
    vx = zeros(n, 1);
    vy = zeros(n, 1);
    y  = zeros(n, 1);
    ax = zeros(n, 1);
    
    % Cinematica link 1
    vx(1) = v_x0 + L(1) * dq(1) * cos(q(1));
    vy(1) = L(1) * dq(1) * sin(q(1));
    y(1)  = -L(1) * cos(q(1));
    ax(1) = ddx0 + L(1) * (ddq(1) * cos(q(1)) - dq(1)^2 * sin(q(1)));
    
    % Cinematica link successivi
    for i = 2:n
        vx(i) = vx(i-1) + L(i) * dq(i) * cos(q(i));
        vy(i) = vy(i-1) + L(i) * dq(i) * sin(q(i));
        y(i)  = y(i-1) - L(i) * cos(q(i));
        ax(i) = ax(i-1) + L(i) * (ddq(i) * cos(q(i)) - dq(i)^2 * sin(q(i)));
    end
    
    % Calcolo Energia
    T_tot(k) = 0.5 * sum(m .* (vx.^2 + vy.^2));
    V_tot(k) = sum(m .* g .* y);
    E_tot(k) = T_tot(k) + V_tot(k);
    
    % Calcolo Potenza istantanea
    F_base_su_catena = sum(m .* ax);
    P_base(k) = F_base_su_catena * v_x0;
end

% Calcolo Lavoro ed Energia Attesa
W_base = cumtrapz(t_out, P_base);
E_attesa = E_tot(1) + W_base;

% --- NORMALIZZAZIONE ---
% Calcolo dell'Energia del pendolo orizzontale
E_rif = g * sum(mu .* L);

E_norm = E_tot / E_rif;
E_attesa_norm = E_attesa / E_rif;
E_initial_norm = E_norm(1);

% --- PLOT ENERGIA ---
figure('Name', 'Energia Normalizzata', 'Color', 'w', 'Position', [100, 100, 800, 450]);
hold on; grid on; grid minor;
set(gca, 'FontSize', 12, 'LineWidth', 1.2, 'TickLabelInterpreter', 'latex', ...
    'GridAlpha', 0.15, 'MinorGridAlpha', 0.1);

c_E = [0.200, 0.200, 0.200]; % Grigio Scuro (Simulata)
c_A = [0.466, 0.674, 0.188]; % Verde (Attesa)
c_0 = [0.000, 0.447, 0.741]; % Blu (Iniziale)

plot(t_out, E_norm, 'Color', c_E, 'LineWidth', 2.5, 'DisplayName', '$\tilde{E}_{tot}$ (Energia Simulata)');
plot(t_out, E_attesa_norm, '--', 'Color', c_A, 'LineWidth', 2.5, 'DisplayName', '$\tilde{E}_{attesa}$ (Energia Teorica)');
yline(E_initial_norm, ':', 'Color', c_0, 'LineWidth', 1.5, 'DisplayName', '$\tilde{E}_0$ (Iniziale)');

xlabel('Tempo $t$ [s]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Energia Adimensionale [-]', 'Interpreter', 'latex', 'FontSize', 14);
title('\textbf{Confronto Energia Simulata vs Attesa}', 'Interpreter', 'latex', 'FontSize', 16);

lgd = legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);
lgd.Box = 'off';
xlim([0, t_span(2)]);

% Evita che il grafico esca fuori vista
min_val = min([E_norm, E_attesa_norm, E_initial_norm]);
max_val = max([E_norm, E_attesa_norm, E_initial_norm]);
margin = max(0.1 * abs(max_val - min_val), 0.1); % 10% margin
ylim([min_val - margin, max_val + margin]);


% --- 4. ANIMAZIONE ---
disp('Avvio animazione...');

x_anim = zeros(n + 1, length(t_out));
y_anim = zeros(n + 1, length(t_out));

x_anim(1, :) = A * sin(w * t_out');
y_anim(1, :) = 0;

for i = 1:n
    x_anim(i+1, :) = x_anim(i, :) + L(i) * sin(theta_out(i, :));
    y_anim(i+1, :) = y_anim(i, :) - L(i) * cos(theta_out(i, :));
end

figure('Name', 'Animazione n-Pendolo su Base Oscillante', 'Color', 'w', 'Position', [950, 100, 600, 600]);
ax = axes; hold on; axis equal; grid on;

lim_x = sum(L) + A + 0.5;
lim_y = sum(L) + 0.5;
xlim(ax, [-lim_x, lim_x]);
ylim(ax, [-lim_y, lim_y]);

xlabel('X [m]'); ylabel('Y [m]');
title(sprintf('N-Pendolo Dinamico (n=%d)', n));

plot([-A-0.5, A+0.5], [0, 0], 'k--', 'LineWidth', 1);
h_line = plot(ax, x_anim(:, 1), y_anim(:, 1), 'o-k', 'LineWidth', 2.5, ...
    'MarkerSize', 8, 'MarkerFaceColor', 'r');

disp('Premi un tasto qualsiasi per avviare l''animazione...');
pause;

% Loop animazione sincronizzato
tic;
for k = 1:15:length(t_out)
    set(h_line, 'XData', x_anim(:, k), 'YData', y_anim(:, k));
    drawnow;
    
    % Sincronizza il tempo di simulazione con il tempo reale
    tempo_simulato = t_out(k);
    while toc < tempo_simulato
        pause(0.002);
    end
end
disp('Processo terminato.');

% --- 5. MOTORE FISICO ---
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