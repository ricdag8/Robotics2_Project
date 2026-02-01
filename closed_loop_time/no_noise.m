clear; clc; close all;

%% ========================================================================
%  1. PARAMETRI DEL SISTEMA E CONTROLLO
% ========================================================================
Mm      = 1.0;        % Inerzia lato motore [kg m^2]
M       = 2.0;        % Inerzia lato link [kg m^2]
K       = 1000.0;     % Rigidezza elastica [Nm/rad]
Ts      = 0.001;      % Passo di campionamento [s]
    % Durata simulazione [s]
K_p_theta = 7.07;     
K_d_theta = 11.93;      
K_p_tau   = 3.11;     
K_d_tau   = 0.91;   
 

% Parametri Rumore (Definiti qui per chiarezza)
STD_NOISE_THETA = 0;  % Deviazione standard rumore theta
STD_NOISE_TAU   = 0.00;   % Deviazione standard rumore tau_J
VAR_NOISE_THETA = STD_NOISE_THETA^2;
VAR_NOISE_TAU   = STD_NOISE_TAU^2;


% Assegnazione al Workspace per Simulink
assignin('base', 'K', K);
assignin('base', 'Kp_theta', K_p_theta);
assignin('base', 'Kd_theta', K_d_theta);
assignin('base', 'Kp_tau',   K_p_tau);
assignin('base', 'Kd_tau',   K_d_tau);
assignin('base', 'VAR_NOISE_THETA', VAR_NOISE_THETA);
assignin('base', 'VAR_NOISE_TAU',   VAR_NOISE_TAU);
assignin('base', 'Ts', Ts);




Ts = 0.001;           
T_total = 40.0;        
t_vec   = (0:Ts:T_total)';
syms t 
pi_sym = sym(pi);


q_sym = (0.6) * sin((pi_sym/7) * t);


q_dot_sym = diff(q_sym, t);
q_ddot_sym = diff(q_dot_sym, t);


func_q      = matlabFunction(q_sym);
func_q_dot  = matlabFunction(q_dot_sym);
func_q_ddot = matlabFunction(q_ddot_sym);


q_vals      = func_q(t_vec);
q_dot_vals  = func_q_dot(t_vec);
q_ddot_vals = func_q_ddot(t_vec);


figure; 
subplot(3,1,1); plot(t_vec, q_vals); title('Posizione'); grid on;
subplot(3,1,2); plot(t_vec, q_dot_vals); title('Velocità'); grid on;
subplot(3,1,3); plot(t_vec, q_ddot_vals); title('Accelerazione'); grid on;



tau_Jd_vec = M * q_ddot_vals;



theta_d_vec = (tau_Jd_vec/K) + q_vals;


 
theta_d_input = [t_vec, theta_d_vec];
tau_jd_input  = [t_vec, tau_Jd_vec];

assignin('base','theta_d_input',theta_d_input);
assignin('base','tau_jd_input', tau_jd_input);



mdl = 'modello_controllore3';  

disp('Fase 3: Esecuzione della simulazione Simulink...');

try
    simOut = sim(mdl, 'StopTime', num2str(T_total), ...
                 'SaveOutput','on', 'SignalLogging','on');
    fprintf('Simulazione completata.\n');

    % Uscite dai blocchi To Workspace (timeseries)
    theta_ts = simOut.theta;
    tauJ_ts  = simOut.tau_J;

    theta_data = theta_ts.Data;
    t_theta    = theta_ts.Time;

    tau_J_data = tauJ_ts.Data;
    t_tau_J    = tauJ_ts.Time;

    % Controllo al motore (devi avere un To Workspace chiamato tau_ctrl)
    if isprop(simOut,'tau_ctrl')
        tau_ctrl_ts   = simOut.tau_ctrl;
    else
        error('Manca il segnale "tau_ctrl" in simOut. Aggiungi un blocco To Workspace con quel nome.');
    end

    tau_u_data = tau_ctrl_ts.Data;
    t_tau_u    = tau_ctrl_ts.Time;

catch ME
    fprintf('\n\n !!! SIMULAZIONE FALLITA !!! \n');
    disp(ME.message);
    error('Controlla nomi dei segnali "theta", "tau_J" e "tau_ctrl" nei blocchi To Workspace.');
end

t_meas      = t_theta(:);
theta_meas  = theta_data(:);
tau_J_meas  = tau_J_data(:);
tau_u_val  = tau_u_data(:);


%% ========================================================================
%  TAU_J AND THETA_DOT DERIVATIVES CALCULATION WITH PRE-FILTERING, THIS
%  COMPUTATION IS ONLY PERFORMED FOR PLOTS PURPOSES. A NEW COMPUTATION WILL
%  BE DISPLAYED NEXT, BUT THE STEPS AND RESULTS ARE THE SAME
% ========================================================================
dt = mean(diff(t_vec));
N = length(t_vec);

% since the main signal component is at 0.7 Hz, we can set a cutoff
% frequency slightly above that to preserve the signal while attenuating noise.
 

theta_clean =  theta_meas;



% 3. finite derivative taking a total of 7 points (3 before and 3 after)
% only after we have the clean signal we can procede to calculate the derivatives
theta_ddot_final = zeros(N, 1);
coeffs_d2 = [2; -27; 270; -490; 270; -27; 2];

for i = 4 : (N - 3)
    idx = (i-3) : (i+3);
    
    th_win  = theta_clean(idx); 
    theta_ddot_final(i) = (th_win' * coeffs_d2) / (180 * dt^2);
end

theta_ddot_calc = fillmissing(theta_ddot_final, 'nearest');




tau_J_meas_clean =  tau_J_meas;





%% ========================================================================
%  PLOT DI CONFRONTO
% ========================================================================
figure('Name', 'Metodo Pre-Filtering', 'Color', 'w');
subplot(2,1,1);
plot(t_vec, theta_meas, 'Color', [0.7 0.7 0.7]); hold on;
plot(t_vec, theta_clean, 'b', 'LineWidth', 1.5);
title('Posizione: Misurata (Rumorosa) vs Filtrata');
legend('Theta Meas', 'Theta Clean'); grid on;

subplot(2,1,2);
plot(t_vec, theta_ddot_final, 'b', 'LineWidth', 1.5); hold on;
plot(t_vec, q_ddot_vals, 'g--', 'LineWidth', 1.5); % Confronto con ideale
title('Accelerazione Calcolata (da segnale pre-filtrato)');
legend('Accelerazione Calcolata', 'Accelerazione Ideale');
grid on;
fc = 5;
fs = 1/dt;
[b, a] = butter(2, fc/(fs/2), 'low');
theta_ddot_filt = filtfilt(b, a, theta_ddot_calc);

q_reconstructed = theta_clean - (tau_J_meas_clean / K);

tau_elastic_term = K * (theta_clean - q_reconstructed);

tau_reconstructed = (Mm * theta_ddot_filt) + tau_elastic_term;

buffer = 50;
valid_idx = buffer : (N - buffer);

t_plot = t_vec(valid_idx);
tau_real = tau_u_val(valid_idx);
tau_est  = tau_reconstructed(valid_idx);

err_tau = tau_real - tau_est;
rmse_tau = sqrt(mean(err_tau.^2));
fit_percent = 100 * (1 - norm(err_tau)/norm(tau_real - mean(tau_real)));

figure('Name', 'Validazione Ricostruzione Coppia', 'Color', 'w');

subplot(2,1,1);
plot(t_plot, tau_real, 'b', 'LineWidth', 1.5); hold on;
plot(t_plot, tau_est, 'r--', 'LineWidth', 1.5);
title(['Confronto Coppia Motore: Reale vs Ricostruita (Fit: ' num2str(fit_percent, '%.2f') '%)']);
ylabel('Coppia [Nm]');
legend('Tau Controllo (Sensore)', 'Tau Ricostruita (Modello Inverso)');
grid on;

subplot(2,1,2);
plot(t_plot, err_tau, 'k');
title(['Errore di Ricostruzione (RMSE: ' num2str(rmse_tau, '%.4f') ' Nm)']);
xlabel('Tempo [s]'); ylabel('Errore [Nm]');
grid on;




%% ========================================================================
%  PLOT DI ANALISI DELLE ACCELERAZIONI (Raw vs Filtered)
% ========================================================================
figure('Name', 'Analisi Derivata Numerica', 'Color', 'w');
hold on; grid on;

% 1. Plot dell'accelerazione "Grezza" (Calcolata con Stencil 7-punti)
%    La facciamo in grigio chiaro per vedere quanto "sporca" è sotto.
h_raw = plot(t_vec(valid_idx), theta_ddot_calc(valid_idx), ...
             'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);

% 2. Plot dell'accelerazione "Filtrata" (Butterworth + Filtfilt)
%    Questa è quella usata per ricostruire la coppia.
h_filt = plot(t_vec(valid_idx), theta_ddot_filt(valid_idx), ...
              'b', 'LineWidth', 2);

% 3. Estetica e Legenda
xlabel('Tempo [s]', 'FontSize', 12);
ylabel('Accelerazione Angolare [rad/s^2]', 'FontSize', 12);
title('Effetto del Filtraggio sulla Derivata Numerica', 'FontSize', 14);
legend([h_raw, h_filt], ...
       'Accelerazione Grezza (Noise)', ...
       'Accelerazione Filtrata (Clean)', ...
       'Location', 'best');

% Opzionale: Zoom automatico su una porzione per vedere i dettagli
% xlim([0, 5]); % Decommenta se vuoi vedere solo i primi 5 secondi







% Creazione della matrice unica [Tempo, Theta, Tau_J]
% Questa è la matrice che userai per il filtraggio
matrice_dati_sensori = [t_meas, theta_meas, tau_J_meas];
theta_d_resampled = interp1(t_vec, theta_d_vec, t_theta, 'linear', 'extrap');

% Calcolo dell'errore di inseguimento
error_theta = theta_d_resampled - theta_data;

% --- FIGURA 1: Tracking di Posizione (Motore) ---
figure('Name', 'Tracking Posizione Motore', 'Color', 'w');
subplot(2,1,1);
    plot(t_vec, theta_d_vec, 'r--', 'LineWidth', 1.5); hold on;
    plot(t_theta, theta_data, 'b', 'LineWidth', 1.2);
    grid on; grid minor;
    ylabel('Posizione [rad]');
    legend('Theta Desiderato (Ref)', 'Theta Misurato (Sim)', 'Location', 'best');
    title('Inseguimento di Traiettoria: \theta_{des} vs \theta_{meas}');

subplot(2,1,2);
    plot(t_theta, error_theta, 'k', 'LineWidth', 1.0);
    grid on; grid minor;
    ylabel('Errore [rad]');
    xlabel('Tempo [s]');
    title(['Errore di Inseguimento (Max: ' num2str(max(abs(error_theta)), '%.4f') ' rad)']);

% --- FIGURA 2: Analisi delle Coppie (Sforzo di controllo) ---
figure('Name', 'Analisi Coppie', 'Color', 'w');
subplot(2,1,1);
    plot(t_tau_u, tau_u_data, 'b', 'LineWidth', 1.2); hold on;
    % Aggiungo limiti di saturazione (se ne hai, es. +/- 100Nm) per riferimento visivo
    yline(100, 'r--'); yline(-100, 'r--'); 
    grid on;
    ylabel('Coppia Motore [Nm]');
    title('Coppia di Controllo (\tau_{ctrl})');
    legend('Uscita PID', 'Limiti (ipotetici)');

subplot(2,1,2);
    plot(t_tau_J, tau_J_data, 'g', 'LineWidth', 1.2); hold on;
    plot(t_vec, tau_Jd_vec, 'm--', 'LineWidth', 1.5);
    grid on;
    ylabel('Coppia Giunto [Nm]');
    xlabel('Tempo [s]');
    legend('Tau Giunto (Misurata)', 'Tau Giunto (Ideale)', 'Location', 'best');
    title('Coppia Trasmessa attraverso la Molla');

% --- FIGURA 3: Effetto Elastico (Differenza Motore-Link) ---
% Nota: Questo grafico ha senso se estrai anche 'q' (lato link) da Simulink.
% Se non lo hai, visualizziamo la deformazione della molla stimata: (tau_J / K)
deformazione_molla = tau_J_data / K;

figure('Name', 'Effetti Elastici', 'Color', 'w');
plot(t_tau_J, deformazione_molla * 180/pi, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
grid on;
title('Deformazione angolare della molla (\theta - q)');
ylabel('Deflessione [gradi]');
xlabel('Tempo [s]');
subtitle('Quanto "flette" il giunto durante il movimento');




%% ========================================================================
%  ELABORAZIONE DATI: FILTRAGGIO -> DERIVAZIONE -> MATRICE
% ========================================================================

% prepare the data
t_proc      = t_theta(:);
theta_raw   = theta_data(:);
tau_J_raw   = tau_J_data(:);
tau_u_raw   = tau_u_data(:);

dt = mean(diff(t_proc)); 
N_samples = length(t_proc);



% we apply the filter
theta_filt =  theta_raw;
tau_J_filt =  tau_J_raw;


tau_u_filt = tau_u_raw; % not needed 



% coefficients for first and second order derivative
coeffs_d1 = [-1; 9; -45; 0; 45; -9; 1];       
coeffs_d2 = [2; -27; 270; -490; 270; -27; 2]; 


th_dot   = zeros(N_samples, 1);
th_ddot  = zeros(N_samples, 1);
tau_dot  = zeros(N_samples, 1);

for i = 4 : (N_samples - 3)
    idx = (i-3) : (i+3);
    
    % we take the filtered data and then we compute the derivative
    win_th  = theta_filt(idx);
    win_tau = tau_J_filt(idx);
    
   
    th_dot(i)  = (win_th' * coeffs_d1) / (60 * dt);
    th_ddot(i) = (win_th' * coeffs_d2) / (180 * dt^2);
    tau_dot(i) = (win_tau' * coeffs_d1) / (60 * dt);
end

%interpolate the missing values which do not have neighbor. basically the
%last values will take the value of the nearest neighbor.
th_dot  = fillmissing(th_dot, 'nearest');
th_ddot = fillmissing(th_ddot, 'nearest');
tau_dot = fillmissing(tau_dot, 'nearest');

%we can now start to build the regressor matrix
buffer = 50; 
valid_idx = buffer : (N_samples - buffer);

% Matrix columns: 
% [1:Time, 2:Theta_Filt, 3:Tau_J_Filt, 4:Th_dot, 5:Th_ddot, 6:Tau_dot, 7:Tau_u_Raw(o Filt)]
Matrix_LS_Ready = [ ...
    t_proc(valid_idx), ...      % 1. Tempo
    theta_filt(valid_idx), ...  % 2. Posizione (FILTRATA)
    tau_J_filt(valid_idx), ...  % 3. Coppia Giunto (FILTRATA)
    th_dot(valid_idx), ...      % 4. Velocità Motore (Calcolata su filtrato)
    th_ddot(valid_idx), ...     % 5. Accel Motore (Calcolata su filtrato)
    tau_dot(valid_idx), ...     % 6. Derivata Coppia (Calcolata su filtrato)
    tau_u_filt(valid_idx)       % 7. Coppia Controllo (Target per LS)
];

fprintf('Elaborazione completata.\n');
fprintf('Dati filtrati all''origine e poi derivati.\n');
fprintf('Matrice pronta: %d righe x 7 colonne.\n', size(Matrix_LS_Ready,1));



%% ========================================================================
%  5. LEAST SQUARES
% ========================================================================

% data estraction
% Matrix_LS_Ready = [t, theta, tau, th_dot, th_ddot, tau_dot, tau_ddot]
t_est       = Matrix_LS_Ready(:, 1);
theta_meas  = Matrix_LS_Ready(:, 2);
tau_meas    = Matrix_LS_Ready(:, 3);
th_dot_meas = Matrix_LS_Ready(:, 4);
th_ddot_meas= Matrix_LS_Ready(:, 5);
tau_dot_meas= Matrix_LS_Ready(:, 6);

% just for safety we rebuild the reference value  
syms t_sym_k
pi_sym = sym(pi);
q_sym_ref = (0.6) * sin((pi_sym/7) * t_sym_k);

f_q      = matlabFunction(q_sym_ref);
f_q_ddot = matlabFunction(diff(diff(q_sym_ref, t_sym_k), t_sym_k));

% Valori ideali agli istanti misurati
q_d      = f_q(t_est);
q_ddot_d = f_q_ddot(t_est);

% Inversione Dinamica (Calcolo Riferimenti Controllo)
tau_Jd_ref  = M * q_ddot_d;                 % Coppia Giunto Desiderata
theta_d_ref = q_d + (tau_Jd_ref / K); % Posizione Motore Desiderata

% 3. COSTRUZIONE VETTORE USCITA (Y)
% Poiché NON c'è Feedforward, la coppia totale del motore è generata SOLO dal PID.
% Y = Coppia Totale Motore = Mm * Accel + Tau_Carico
Y_target = (Mm * th_ddot_meas) + tau_meas;

% 4. COSTRUZIONE MATRICE REGRESSORI (Phi)
% Modello: Tau = Kp_tau*err_tau - Kd_tau*d_tau + Kp_th*err_th - Kd_th*d_th

% Col 1: Errore Coppia (Kp_tau)
col_Kp_tau   = tau_Jd_ref - tau_meas;

% Col 2: -Derivata Coppia (Kd_tau)
col_Kd_tau   = -tau_dot_meas;

% Col 3: Errore Posizione (Kp_theta)
col_Kp_theta = theta_d_ref - theta_meas;

% Col 4: -Velocità Motore (Kd_theta)
col_Kd_theta = -th_dot_meas;

% Assemblaggio Phi
Phi = [col_Kp_tau, col_Kd_tau, col_Kp_theta, col_Kd_theta];

% 5. SOLUZIONE LEAST SQUARES
Gains_est = pinv(Phi)*Y_target;

% 6. RISULTATI
Kpt_est  = Gains_est(1);
Kdt_est  = Gains_est(2);
Kpth_est = Gains_est(3);
Kdth_est = Gains_est(4);

fprintf('\n------------------------------------------------\n');
fprintf(' RISULTATI STIMA (SENZA FEEDFORWARD)\n');
fprintf('------------------------------------------------\n');
fprintf('Kp_tau   (Reale: %.2f) -> Stimato: %.4f\n', K_p_tau, Kpt_est);
fprintf('Kd_tau   (Reale: %.2f) -> Stimato: %.4f\n', K_d_tau, Kdt_est);
fprintf('Kp_theta (Reale: %.2f) -> Stimato: %.4f\n', K_p_theta, Kpth_est);
fprintf('Kd_theta (Reale: %.2f) -> Stimato: %.4f\n', K_d_theta, Kdth_est);
fprintf('------------------------------------------------\n');

% Validazione Grafica
Y_model = Phi * Gains_est;
fit_perc = 100 * (1 - norm(Y_target - Y_model)/norm(Y_target - mean(Y_target)));

figure('Name','Gains estimate','Color','w');
subplot(2,1,1);
plot(t_est, Y_target, 'b', 'LineWidth', 1.5); hold on;
plot(t_est, Y_model, 'r--', 'LineWidth', 1.5);
legend('Total motor torque (Y)', 'Reconstructed tau from the computed gains');
title(['Model fit (Fit: ' num2str(fit_perc,'%.2f') '%)']);
ylabel('Torque [Nm]'); grid on;




















