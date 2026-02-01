clear all; clc; close all;

%% ========================================================================
%  1. PARAMETRI DEL SISTEMA E CONTROLLO (MATLAB PURO)
% ========================================================================
P.Mm      = 1.028;        % Inerzia lato motore [kg m^2]
P.M       = 2.021;        % Inerzia lato link [kg m^2]
P.K       = 1000.0;       % Rigidezza elastica [Nm/rad]
Ts        = 0.001;        % Passo di campionamento per l'analisi FFT

% Parametri per l'estrazione finale (alias per coerenza con le tue formule)
Mm = P.Mm;
M  = P.M;
K_molla = P.K;

%questi sono i guadagni su cui funziona anche la stima normale 
%P.Kp_theta = 7.07;     
%P.Kd_theta = 18.93;      
%P.Kp_tau   = 8.41;     
%P.Kd_tau   = 0.4;  


P.Kp_theta = 70.07;     
P.Kd_theta = 46.93;      
P.Kp_tau   = 6.41;     
P.Kd_tau   = 0.9;  



% Per confronto finale
K_p_theta = P.Kp_theta;
K_d_theta = P.Kd_theta;
K_p_tau   = P.Kp_tau;
K_d_tau   = P.Kd_tau;

%% ========================================================================
%  FASE 2: SIMULAZIONE SEQUENZIALE (GENERAZIONE DATI CON ODE45)
% ======================================================================== 
disp('--- Inizio Scansione Frequenze (MATLAB PURO) ---');

% Scansione in frequenza
frequncy_to_test = logspace(log10(0.1), log10(100), 50); 
N_test = length(frequncy_to_test);

Freqs_rad = zeros(N_test, 1);  
H_complex = zeros(N_test, 1);  

% Limiti fisici per il calcolo dell'ampiezza dinamica
MAX_POS_AMP = 0.8;  
MAX_TAU_AMP = 10.0; 

for k = 1:N_test
    w_test = frequncy_to_test(k); 
    
    % Calcolo ampiezza sicura
    Amp_q_limit_kinematic = MAX_POS_AMP;
    Amp_q_limit_dynamic   = MAX_TAU_AMP / (P.M * w_test^2);
    Amp_q = min(Amp_q_limit_kinematic, Amp_q_limit_dynamic);
    
    T_period = 2*pi / w_test;
    
    if w_test < 1
        Num_Periodi_Transitorio = 3; Num_Periodi_Regime = 3;  
    elseif w_test < 15.0
        Num_Periodi_Transitorio = 3; Num_Periodi_Regime = 10;
    else
        Num_Periodi_Transitorio = 50; Num_Periodi_Regime = 100; 
    end
    
    T_transitorio = Num_Periodi_Transitorio * T_period;
    T_regime      = Num_Periodi_Regime * T_period;
    Tf = T_transitorio + T_regime;
    t_vec = 0:Ts:Tf;
    
    % --- SIMULAZIONE ODE45 ---
    x0 = [0; 0; 0; 0];
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    [t_ode, x_ode] = ode45(@(t,x) system_dynamics(t, x, P, Amp_q, w_test), [0 Tf], x0, opts);
    
    % --- RICOSTRUZIONE SEGNALI ---
    % 1. Output: Coppia Misurata (Tau_J)
    % Interpoliamo gli stati sulla griglia fissa t_vec
    theta_interp = interp1(t_ode, x_ode(:,1), t_vec, 'spline');
    q_interp     = interp1(t_ode, x_ode(:,3), t_vec, 'spline');
    tau_meas     = P.K * (theta_interp - q_interp); % Coppia misurata
    
    % 2. Input: Coppia di Riferimento (Tau_ref = M * q_ddot_des)
    % Calcolata analiticamente per precisione
    acc_vals_exact = Amp_q * w_test^2 * sin(w_test * t_vec);
    tau_ref        = P.M * acc_vals_exact; 
    
    % --- CALCOLO FFT ---
    mask_steady = t_vec >= T_transitorio;
    
    % Input/Output per la stima (H = Misurata / Riferimento)
    u_steady = tau_ref(mask_steady);  
    y_steady = tau_meas(mask_steady); 
    
    % Gestione lunghezze
    len = min(length(u_steady), length(y_steady));
    u_steady = u_steady(1:len);
    y_steady = y_steady(1:len);

    FFT_u = fft(u_steady);
    FFT_y = fft(y_steady);
    N_samples = length(y_steady);
    
    [~, idx_peak] = max(abs(FFT_u(1:floor(N_samples/2)+1)));
    
    U_phasor = FFT_u(idx_peak);
    Y_phasor = FFT_y(idx_peak); 
    
    H_complex(k) = Y_phasor / U_phasor;
    Freqs_rad(k) = w_test;
    
    % Feedback progresso
    if mod(k,10)==0, fprintf('Simulazione freq %.1f rad/s completata.\n', w_test); end
end

%% ========================================================================
%  FASE 5: FITTING FDT CON LEAST SQUARES (WLS) - PESI AGGRESSIVI
% ========================================================================
disp('--- Fitting Funzione di Trasferimento (High Precision Low-Freq) ---');
na = 4; % Poli
nb = 2; % Zeri
N = length(Freqs_rad);
num_vars = (nb + 1) + na; 
mag_data = abs(H_complex);


A_mat = zeros(2 * N, num_vars);
B_vec = zeros(2 * N, 1);
for k = 1:N
    s_val = 1j * Freqs_rad(k);
    h_val = H_complex(k);
    
    % Numeratore: [1, s, s^2]
    reg_num = zeros(1, nb+1);
    for i = 0:nb, reg_num(i+1) = s_val^i; end
    
    % Denominatore: [-h*s, -h*s^2, -h*s^3, -h*s^4]
    reg_den = zeros(1, na);
    for i = 1:na, reg_den(i) = -h_val * (s_val^i); end
    
    row_complex = [reg_num, reg_den];
    rhs_complex = h_val;
    
    A_mat(2*k-1, :) = real(row_complex);
    B_vec(2*k-1)    = real(rhs_complex);
    A_mat(2*k, :)   = imag(row_complex);
    B_vec(2*k)      = imag(rhs_complex);
end

% Risoluzione
cond_num = cond(A_mat);
fprintf('Condizionamento Matrice: %.2e\n', cond_num);

Theta = (A_mat ) \ (B_vec );

% Creazione Modello
b_ls = Theta(1 : nb+1);       
a_ls = [1; Theta(nb+2 : end)]; 
sys_est = tf(flip(b_ls).', flip(a_ls).');
disp('Fitting Completato.');

%% ========================================================================
%  FASE 6: VALIDAZIONE GRAFICA
% ========================================================================
w_ver = logspace(log10(min(Freqs_rad)), log10(max(Freqs_rad)), 500);
[mag, phase] = bode(sys_est, w_ver);

figure('Name','Fit LS Vincolato (Da Sinusoidi Sequenziali)','Color','w');
subplot(2,1,1); 
semilogx(Freqs_rad, 20*log10(abs(H_complex)), 'bo', 'LineWidth', 1.5); hold on;
semilogx(w_ver, 20*log10(squeeze(mag)), 'r-', 'LineWidth', 2);
grid on; ylabel('Magnitudo [dB]'); title('Fit FDT');
legend('Dati Sperimentali (ODE45)', 'Modello Identificato');

% Linea teorica risonanza fisica (per vedere dove cade)
omega_n = sqrt(K_molla/M);
xline(omega_n, 'g--', 'Label', 'Risonanza Fisica');

subplot(2,1,2); 
ph_data = rad2deg(unwrap(angle(H_complex)));
ph_model = rad2deg(unwrap(deg2rad(squeeze(phase))));
off = round((ph_data(1)-ph_model(1))/360)*360;
semilogx(Freqs_rad, ph_data, 'bo', 'LineWidth', 1.5); hold on;
semilogx(w_ver, ph_model + off, 'r-', 'LineWidth', 2);
grid on; ylabel('Fase [deg]'); xlabel('Rad/s');

%% ========================================================================
%  FASE 7: ESTRAZIONE DEI GUADAGNI DEL CONTROLLORE
% ========================================================================
disp('--- Estrazione Parametri Controllore ---');

% 1. Recupero coefficienti dal modello identificato
[~, den_est] = tfdata(sys_est, 'v'); 
den_est = den_est / den_est(1); % Normalizzazione s^4
% Mappatura: den = [1, a3, a2, a1, a0]
a3 = den_est(2);
a2 = den_est(3);
a1 = den_est(4);
a0 = den_est(5);

% 2. Calcolo Inverso dei Guadagni
% K_p_theta (da a0)
Kp_theta_est = (a0 * M * Mm) / K_molla;

% K_d_theta (da a1)
Kd_theta_est = (a1 * M * Mm) / K_molla;

% K_d_tau (da a3) -> a3 * Mm = K * Kd_tau + Kd_theta
Kd_tau_est = (a3 * Mm - Kd_theta_est) / K_molla;

% K_p_tau (da a2)
term_mech = K_molla * (M + Mm) / M;
Kp_tau_est = (a2 * Mm - Kp_theta_est - term_mech) / K_molla;

% 3. Visualizzazione Risultati
fprintf('\n--- RISULTATI IDENTIFICAZIONE GUADAGNI ---\n');
fprintf('K_p_theta: Reale = %.4f | Stimato = %.4f (Err: %.1f%%)\n', ...
    K_p_theta, Kp_theta_est, abs((Kp_theta_est-K_p_theta)/K_p_theta)*100);
fprintf('K_d_theta: Reale = %.4f | Stimato = %.4f (Err: %.1f%%)\n', ...
    K_d_theta, Kd_theta_est, abs((Kd_theta_est-K_d_theta)/K_d_theta)*100);
fprintf('K_p_tau:   Reale = %.4f | Stimato = %.4f (Err: %.1f%%)\n', ...
    K_p_tau, Kp_tau_est, abs((Kp_tau_est-K_p_tau)/K_p_tau)*100);
fprintf('K_d_tau:   Reale = %.4f | Stimato = %.4f (Err: %.1f%%)\n', ...
    K_d_tau, Kd_tau_est, abs((Kd_tau_est-K_d_tau)/K_d_tau)*100);
disp('--------------------------------------------');

%% ========================================================================
%  FUNZIONE DINAMICA DEL SISTEMA
% ========================================================================
function dxdt = system_dynamics(t, x, P, Amp, w)
    th      = x(1); th_dot  = x(2);
    q       = x(3); q_dot   = x(4);

    % Riferimenti Analitici
    sin_wt = sin(w*t); cos_wt = cos(w*t);
    
    q_d         = -Amp * sin_wt;
    q_d_dot     = -Amp * w * cos_wt;
    q_d_ddot    =  Amp * w^2 * sin_wt;
    q_d_dddot   =  Amp * w^3 * cos_wt;
    
    % Feedforward Coppia
    tau_Jd      = P.M * q_d_ddot;
    tau_Jd_dot  = P.M * q_d_dddot;
    
    % Riferimento Motore
    theta_d     = q_d + tau_Jd / P.K;
    theta_d_dot = q_d_dot + tau_Jd_dot / P.K;

    % Variabili Interne
    tau_J     = P.K * (th - q);
    tau_J_dot = P.K * (th_dot - q_dot); 

    % Controllo
    e_tau   = tau_Jd - tau_J;
    e_theta = theta_d - th;
    
    u_tau   = P.Kp_tau * e_tau     - P.Kd_tau * tau_J_dot;
    u_theta = P.Kp_theta * e_theta - P.Kd_theta * th_dot;
    
    tau_motor = u_tau + u_theta;

    % Dinamica
    th_ddot = (tau_motor - tau_J) / P.Mm;
    q_ddot  = tau_J / P.M;

    dxdt = [th_dot; th_ddot; q_dot; q_ddot];
end




%% ========================================================================
%  FASE DI VALIDAZIONE: TRAIETTORIA REST-TO-REST (Quintica)
% ========================================================================
if ~exist('sys_est', 'var')
    error('Attenzione: Devi prima eseguire l''identificazione per avere sys_est!');
end
disp('--- Inizio Validazione Rest-to-Rest (MATLAB PURO) ---');

% 1. Parametri della Traiettoria
Ampiezza_Movimento = 0.9;   % [rad] Spostamento totale
Durata_Movimento   = 3;   % [s]
Durata_Pausa       = 0.5;   % [s]
T_total            = Durata_Movimento + Durata_Pausa;
dt_val             = 0.001; 
t_val              = (0:dt_val:T_total)'; % <--- NOTA: L'apice ' forza Vettore Colonna subito

% 2. Simulazione del "SISTEMA REALE" (Ground Truth) con ODE45
x0 = [0; 0; 0; 0]; 
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

fprintf('Validazione con K fisica = %.1f\n', P.K);
[t_ode, x_ode] = ode45(@(t,x) dynamics_rest_to_rest(t, x, P, Ampiezza_Movimento, Durata_Movimento), ...
                       [0 T_total], x0, opts);

% 3. Ricostruzione dei Segnali (Post-Processing)
% Poiché t_val è colonna, interp1 restituisce colonne
theta_meas = interp1(t_ode, x_ode(:,1), t_val, 'spline');
q_meas     = interp1(t_ode, x_ode(:,3), t_val, 'spline');

% Ricostruiamo l'ingresso (Analitico)
tau_ref_val = zeros(size(t_val));
theta_d_val = zeros(size(t_val)); 

for i = 1:length(t_val)
    [~, sigs] = dynamics_rest_to_rest(t_val(i), [0;0;0;0], P, Ampiezza_Movimento, Durata_Movimento);
    tau_ref_val(i) = sigs.tau_Jd;
    theta_d_val(i) = sigs.theta_d;
end

% Ricostruiamo l'uscita "Reale"
% CORREZIONE: Ora sono tutti vettori colonna, non serve trasporre, o se serve, controlliamo
tau_meas_real = P.K * (theta_meas - q_meas); 

% 4. Simulazione del "MODELLO STIMATO" (lsim)
[y_est, t_est] = lsim(sys_est, tau_ref_val, t_val);

% 5. Calcolo del FIT (Best Fit Percentage)
% CORREZIONE QUI SOTTO: Rimosso il trasposto su y_est per evitare broadcasting
valid_idx = ~isnan(tau_meas_real) & ~isnan(y_est); 

y_real_clean = tau_meas_real(valid_idx);
y_est_clean  = y_est(valid_idx); % Rimosso trasposto inutile

numeratore   = norm(y_real_clean - y_est_clean);
denominatore = norm(y_real_clean - mean(y_real_clean));

if denominatore < 1e-6
    FIT_val = 0; 
else
    FIT_val = 100 * (1 - numeratore/denominatore);
end

fprintf('\n>>> RISULTATO VALIDAZIONE <<<\n');
fprintf('FIT Traiettoria: %.2f %%\n', FIT_val);

% 6. Plotting
figure('Name', 'Validazione Rest-To-Rest', 'Color', 'w');
subplot(3,1,1);
plot(t_val, theta_d_val, 'k--', 'LineWidth', 1); hold on;
plot(t_val, theta_meas, 'b', 'LineWidth', 1.5);
ylabel('Posizione [rad]'); title('Inseguimento Traiettoria (Reale)');
legend('Riferimento', 'Theta Motore'); grid on;

subplot(3,1,2);
plot(t_val, tau_ref_val, 'k', 'LineWidth', 1.5);
ylabel('Input [Nm]'); title('Coppia di Riferimento (Feedforward)');
grid on;

subplot(3,1,3);
plot(t_val, tau_meas_real, 'b', 'LineWidth', 2); hold on;
plot(t_val, y_est, 'r--', 'LineWidth', 2);
ylabel('Output [Nm]'); xlabel('Tempo [s]');
title(['Confronto Coppia Misurata vs Stimata (FIT: ' num2str(FIT_val, '%.2f') '%)']);
legend('Reale (ODE45)', 'Stimato (lsim)'); grid on;

%% --- FUNZIONE DINAMICA SPECIFICA PER POLINOMIO 5° ---
function [dxdt, signals] = dynamics_rest_to_rest(t, x, P, Amp, T_move)
    % Stati
    th = x(1); th_dot = x(2); q = x(3); q_dot = x(4);
    
    % --- GENERAZIONE TRAIETTORIA (Polinomio 5°) ---
    if t <= T_move
        tau = t / T_move; % Tempo normalizzato 0..1
        
        % Posizione: s(t) = A * (10t^3 - 15t^4 + 6t^5)
        poly_pos = 10*tau^3 - 15*tau^4 + 6*tau^5;
        
        % Velocità: v(t) = (A/T) * (30t^2 - 60t^3 + 30t^4)
        poly_vel = (1/T_move) * (30*tau^2 - 60*tau^3 + 30*tau^4);
        
        % Accelerazione: a(t) = (A/T^2) * (60t - 180t^2 + 120t^3)
        poly_acc = (1/T_move^2) * (60*tau - 180*tau^2 + 120*tau^3);
        
        % Jerk (Serve per la derivata della coppia): j(t) = (A/T^3) * (60 - 360t + 360t^2)
        poly_jerk = (1/T_move^3) * (60 - 360*tau + 360*tau^2);
        
        q_d       = Amp * poly_pos;
        q_d_dot   = Amp * poly_vel;
        q_d_ddot  = Amp * poly_acc;
        q_d_dddot = Amp * poly_jerk;
    else
        % Fine movimento: manteniamo la posizione finale
        q_d       = Amp;
        q_d_dot   = 0;
        q_d_ddot  = 0;
        q_d_dddot = 0;
    end
    
    % --- CALCOLO VARIABILI ---
    % Feedforward Coppia
    tau_Jd      = P.M * q_d_ddot;
    tau_Jd_dot  = P.M * q_d_dddot; 
    
    % Riferimento Motore
    theta_d     = q_d + tau_Jd / P.K;
    theta_d_dot = q_d_dot + tau_Jd_dot / P.K;
    
    % Coppia Misurata
    tau_J     = P.K * (th - q);
    tau_J_dot = P.K * (th_dot - q_dot);
    
    % --- CONTROLLORE ---
    e_tau   = tau_Jd - tau_J;
    e_theta = theta_d - th;
    
    % Controllo PD doppio anello
    u_tau   = P.Kp_tau * e_tau     - P.Kd_tau * tau_J_dot;
    u_theta = P.Kp_theta * e_theta - P.Kd_theta * th_dot;
    
    tau_motor = u_tau + u_theta;
    
    % --- DINAMICA ---
    th_ddot = (tau_motor - tau_J) / P.Mm;
    q_ddot  = tau_J / P.M;
    
    dxdt = [th_dot; th_ddot; q_dot; q_ddot];
    
    % Output segnali ausiliari
    signals.tau_Jd  = tau_Jd;
    signals.theta_d = theta_d;
end