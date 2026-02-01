clear; clc; close all;

%% ========================================================================
%  1. SYSTEM AND CONTROLLER PARAMETERS
% ========================================================================
Mm      = 1.0036;       
M       = 2.010;       
K       = 1004.096;     
Ts      = 0.001;  



K_p_theta = 70.07;     
K_d_theta = 46.93;      
K_p_tau   = 6.11;     
K_d_tau   = 0.91;   
 
% Noise parameters
STD_NOISE_THETA = 1e-4;  
STD_NOISE_TAU   = 0.1;    
VAR_NOISE_THETA = STD_NOISE_THETA^2;
VAR_NOISE_TAU   = STD_NOISE_TAU^2;

assignin('base', 'K', K);
assignin('base', 'Mm', Mm);
assignin('base', 'M', M);
assignin('base', 'Kp_theta', K_p_theta);
assignin('base', 'Kd_theta', K_d_theta);
assignin('base', 'Kp_tau',   K_p_tau);
assignin('base', 'Kd_tau',   K_d_tau);
assignin('base', 'VAR_NOISE_THETA', VAR_NOISE_THETA);
assignin('base', 'VAR_NOISE_TAU',   VAR_NOISE_TAU);
assignin('base', 'Ts', Ts);

mdl = 'modello_controllore3'; 
load_system(mdl);


%% ========================================================================
%   SIMULATION 
% ======================================================================== 
disp('--- Scanning frequencies ---');
frequency_to_test = logspace(log10(0.1), log10(40), 50); 
N_test = length(frequency_to_test);

Freqs_rad = zeros(N_test, 1);  
H_complex = zeros(N_test, 1);  

MAX_POS_AMP = 1;  
MAX_TAU_AMP = 15.0; 
Amps_Vector      = zeros(N_test, 1);
Limit_Dyn_Vector = zeros(N_test, 1);
Limit_Kin_Vector = zeros(N_test, 1);

for k = 1:N_test
    w_test = frequency_to_test(k); 
    
    Limit_Tau_Dynamic = MAX_TAU_AMP;
    Limit_Tau_Kinematic = M * w_test^2 * MAX_POS_AMP;
    
    MIN_TORQUE_REQUIRED = 0.15; 
    A_tau = max(MIN_TORQUE_REQUIRED, min(Limit_Tau_Dynamic, Limit_Tau_Kinematic));
    

    Amps_Vector(k)      = A_tau;
    Limit_Dyn_Vector(k) = Limit_Tau_Dynamic;
    Limit_Kin_Vector(k) = Limit_Tau_Kinematic;


    T_period = 2*pi / w_test;
    if w_test < 1
        Num_Periodi_Regime = 3;  
    elseif w_test < 15.0
        Num_Periodi_Regime = 50;
    else
        Num_Periodi_Regime = 60; 
    end
   

    disp(Num_Periodi_Regime);

    T_transitorio = 1 * T_period; % Transitorio sicuro
    Tf = T_transitorio + Num_Periodi_Regime * T_period;
    t_vec = 0:Ts:Tf; 
    
    %input generation
    tau_Jd_vec = A_tau * sin(w_test * t_vec);
    q_vals = -tau_Jd_vec / (M * w_test^2);
    theta_d_vec = q_vals + tau_Jd_vec / K;
    
    theta_d_input = [t_vec(:), theta_d_vec(:)];
    tau_jd_input  = [t_vec(:), tau_Jd_vec(:)];
    
    assignin('base', 'theta_d_input', theta_d_input);
    assignin('base', 'tau_jd_input', tau_jd_input);
    
    try
        simOut = sim(mdl, 'StopTime', num2str(Tf), 'SaveOutput', 'on');
        
        if isfield(simOut, 'tau_J')
             tau_meas = simOut.tau_J.Data;
             t_rec    = simOut.tau_J.Time;
        else 
             logs = simOut.logsout;
             tau_meas = logs.get('tau_J').Values.Data;
             t_rec    = logs.get('tau_J').Values.Time;
        end
        
        tau_ref = A_tau * sin(w_test * t_rec);
        
    catch
        warning('failed simulation %.2f', w_test);
        continue;
    end
    


    y_raw = tau_meas;  
    u_raw = tau_ref;   
    
    u_raw_sync = A_tau * sin(w_test * t_rec); 


    Time_Start_Steady = T_transitorio;
    Time_End_Steady   = t_rec(end);

    t_uniform = (Time_Start_Steady : Ts : Time_End_Steady)';
    

    u_clean = interp1(t_rec, u_raw_sync, t_uniform, 'linear', 'extrap');
    y_clean = interp1(t_rec, y_raw,      t_uniform, 'linear', 'extrap');
     y_noisy   = y_clean   + STD_NOISE_TAU   * an(size(y_clean));

    Samples_Per_Period = (2*pi/w_test) / Ts;

    Num_Full_Cycles = floor( length(t_uniform) / Samples_Per_Period );
    
    if Num_Full_Cycles < 1
        warning('no data for w=%.2f', w_test);
        continue;
    end
    
   
    N_samples_exact = round(Num_Full_Cycles * Samples_Per_Period);
    

    u_final = u_clean(1:N_samples_exact);
    y_final = y_noisy(1:N_samples_exact);
    t_final = t_uniform(1:N_samples_exact);   
    

  
    FFT_u = fft(u_final);
    FFT_y = fft(y_final);
    
    [~, idx_peak] = max(abs(FFT_u));
    H_complex(k) = FFT_y(idx_peak) / FFT_u(idx_peak);
    Freqs_rad(k) = w_test;
    
  
    fprintf('Freq: %.2f | Cicles used: %d (Tagliato da %.2f)\n', ...
            w_test, Num_Full_Cycles, length(t_uniform)/Samples_Per_Period);

    
end



disp('--- Bode diagrams ---');


valid_idx = Freqs_rad > 0;
w_final   = Freqs_rad(valid_idx);
H_final   = H_complex(valid_idx);

if isempty(w_final)
    error('No data for the Bode plots.');
end


Mag_dB    = 20 * log10(abs(H_final));
Phase_deg = rad2deg(unwrap(angle(H_final))); 


figure('Name', 'Experimental Bode', 'Color', 'w');

% --- Magnitude Plot ---
subplot(2,1,1);
semilogx(w_final, Mag_dB, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'MarkerSize', 4);
grid on; grid minor;
ylabel('Magnitude [dB]', 'FontSize', 12);
title('Bode Diagram - Magnitude', 'FontSize', 14);
xlim([min(w_final)*0.9, max(w_final)*1.1]);

% --- Phase Plot ---
subplot(2,1,2);
semilogx(w_final, Phase_deg, 'r-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'MarkerSize', 4);
grid on; grid minor;
xlabel('Frequency [rad/s]', 'FontSize', 12);
ylabel('Phase [deg]', 'FontSize', 12);
title('Bode Diagram - Phase', 'FontSize', 14);
xlim([min(w_final)*0.9, max(w_final)*1.1]);


figure('Name', 'Input Amplitude Profile', 'Color', 'w');


semilogx(frequency_to_test, Amps_Vector, 'b-', 'LineWidth', 2); hold on;
semilogx(frequency_to_test, Limit_Dyn_Vector, 'r--', 'LineWidth', 1.5);
semilogx(frequency_to_test, Limit_Kin_Vector, 'g--', 'LineWidth', 1.5);


yline(0.15, 'k:', 'LineWidth', 1.5, 'Label', 'Min Torque Floor');


grid on; grid minor;
xlabel('Frequency [rad/s]', 'FontSize', 12);
ylabel('Input Torque Amplitude [Nm]', 'FontSize', 12);
title('Excitation Signal Amplitude Selection', 'FontSize', 14);

legend('Selected Amplitude (A_{tau})', ...
       'Dynamic Limit (Motor Saturation)', ...
       'Kinematic Limit (Workspace)', ...
       'Location', 'best');

ylim([0, MAX_TAU_AMP * 1.2]);


%% ========================================================================
%  FASE 3: FITTING WLS (Weighted Least Squares) - ESSENZIALE
% ========================================================================
disp('--- Fitting WLS ---');
na = 4; nb = 2; 
N = length(Freqs_rad);
num_vars = (nb + 1) + na; 

weights = ones(N, 1);
idx_dc = find(Freqs_rad < 3.0);
weights(idx_dc) = 1000; 






W_diag = zeros(2*N, 1);
for k = 1:N
    W_diag(2*k-1) = weights(k);
    W_diag(2*k)   = weights(k);
end

A_mat = zeros(2 * N, num_vars);
B_vec = zeros(2 * N, 1);

for k = 1:N
    s_val = 1j * Freqs_rad(k);
    h_val = H_complex(k);
    
    reg_num = zeros(1, nb+1);
    for i = 0:nb, reg_num(i+1) = s_val^i; end
    
    reg_den = zeros(1, na);
    for i = 1:na, reg_den(i) = -h_val * (s_val^i); end
    
    row_complex = [reg_num, reg_den];
    
    A_mat(2*k-1, :) = real(row_complex);
    B_vec(2*k-1)    = real(h_val);
    A_mat(2*k, :)   = imag(row_complex);
    B_vec(2*k)      = imag(h_val);
end


Theta = (A_mat .* W_diag) \ (B_vec .* W_diag);

b_ls = Theta(1 : nb+1);       
a_ls = [1; Theta(nb+2 : end)]; 

sys_est = tf(flip(b_ls).', flip(a_ls).');
disp('Fitting Completato.');


%% ========================================================================
%   extraction of controller gains
% ========================================================================
disp('--- controller parameter extraction ---');

% 1. retrieval of coefficients from the identified model
% let's ensure we take the normalized denominator (s^4 has coeff 1)
[~, den_est] = tfdata(sys_est, 'v'); 
den_est = den_est / den_est(1); % normalization

% coefficient mapping: den = [1, a3, a2, a1, a0]
a3 = den_est(2);
a2 = den_est(3);
a1 = den_est(4);
a0 = den_est(5);

% 2. physical parameters (already defined, recalling them for clarity)
% ensure they are present in the workspace
% M = 1.0; Mm = 2.0; K = 1000.0; 

% 3. inverse calculation of gains

% k_p_theta (from a0)
Kp_theta_est = (a0 * M * Mm) / K;

% k_d_theta (from a1)
Kd_theta_est = (a1 * M * Mm) / K;

% k_d_tau (from a3)
% note: a3 * mm = k * kd_tau + kd_theta
Kd_tau_est = (a3 * Mm - Kd_theta_est) / K;

% k_p_tau (from a2)
% note: mechanical term to subtract = k * (m + mm) / m
term_mech = K * (M + Mm) / M;
Kp_tau_est = (a2 * Mm - Kp_theta_est - term_mech) / K;

% 4. results visualization
fprintf('\n--- gain identification results ---\n');
fprintf('k_p_theta: real= %.4f | estimated = %.4f\n', K_p_theta, Kp_theta_est);
fprintf('k_d_theta: real = %.4f | estimated = %.4f\n', K_d_theta, Kd_theta_est);
fprintf('k_p_tau:   real = %.4f | estimated = %.4f\n', K_p_tau,   Kp_tau_est);
fprintf('k_d_tau:   real = %.4f | estimated = %.4f\n', K_d_tau,   Kd_tau_est);
disp('--------------------------------------------');

%% ========================================================================
%  extra: zero-pole analysis
% ========================================================================
disp('--- zero and pole analysis ---');

% 1. graphical plot
figure('Name', 'zero-pole map (pz map)', 'Color', 'w');
pzmap(sys_est);
grid on;
title('map of zeros (o) and poles (x) of the estimated model');
axis equal; % keeps real proportions to clearly see damping

% 2. numerical value extraction
Poli = pole(sys_est);
Zeri = zero(sys_est);

fprintf('\n>>> pole list (stability) <<<\n');
disp(Poli);

fprintf('\n>>> zero list (zero dynamics) <<<\n');
disp(Zeri);

% 3. physical analysis (natural frequency and damping)
fprintf('\n>>> physical analysis (natural frequencies and damping) <<<\n');
disp('--- poles (denominator) ---');
damp(Poli);

fprintf('\n--- zeros (numerator) ---');
damp(Zeri);


%% ========================================================================
%   VALIDATION WITH A REST-TO-REST TRAJECTORY
% ========================================================================
disp('--- Starting Rest-to-Rest Validation ---');

Movement_Amplitude = 0.5;   % [rad]
Movement_Duration  = 1.0;   % [s]
Pause_Duration     = 1.5;   % [s]
T_val_total        = Movement_Duration + Pause_Duration;
t_val = 0:Ts:T_val_total;
t_val = t_val(:); % Column vector

% 2. Trajectory generation
% q(t) = a0 + a1*t + a2*t^2 + a3*t^3 + a4*t^4 + a5*t^5
% For Rest-to-Rest (0 -> Amp) in time T, the normalized form (tau = t/T) is:
% s(tau) = 10*tau^3 - 15*tau^4 + 6*tau^5
q_val   = zeros(size(t_val));
v_val   = zeros(size(t_val)); % Velocity
acc_val = zeros(size(t_val)); % Acceleration

% Indices during the movement phase
idx_move = t_val <= Movement_Duration;
tau = t_val(idx_move) / Movement_Duration; % Normalized time 0..1

% Profile calculation
% Position: Distance * (10*tau^3 - 15*tau^4 + 6*tau^5)
q_val(idx_move) = Movement_Amplitude * (10*tau.^3 - 15*tau.^4 + 6*tau.^5);
q_val(~idx_move) = Movement_Amplitude; % Keep final position

% Velocity: (Distance/T) * (30*tau^2 - 60*tau^3 + 30*tau^4)
v_val(idx_move) = (Movement_Amplitude / Movement_Duration) * (30*tau.^2 - 60*tau.^3 + 30*tau.^4);

% Acceleration: (Distance/T^2) * (60*tau - 180*tau^2 + 120*tau^3)
acc_val(idx_move) = (Movement_Amplitude / Movement_Duration^2) * (60*tau - 180*tau.^2 + 120*tau.^3);

% 3. Input and Feedforward Calculation
% Torque needed to move the link inertia
tau_jd_val_vec = M * acc_val;

% Desired Theta (Motor = Link + Spring Deflection)
theta_d_val_vec = q_val + tau_jd_val_vec / K;

% --- TORQUE SAFETY CHECK ---
max_tau_required = max(abs(tau_jd_val_vec));
fprintf('Generated Trajectory:\n');
fprintf('  Displacement: %.2f rad\n', Movement_Amplitude);
fprintf('  Duration: %.2f s\n', Movement_Duration);
fprintf('  Peak Feedforward Torque: %.2f Nm\n', max_tau_required);
if max_tau_required > 10.0
    warning('WARNING');
end

% 4. Timeseries creation
ts_theta_val = timeseries(theta_d_val_vec, t_val);
ts_tau_val   = timeseries(tau_jd_val_vec, t_val);

assignin('base', 'theta_d_input', ts_theta_val);
assignin('base', 'tau_jd_input',  ts_tau_val);

% Make sure the physical parameters are correctly assigned to workspace
assignin('base', 'K', K); 
assignin('base', 'Mm', Mm);
assignin('base', 'M', M);

% 5. REAL SYSTEM SIMULATION (Simulink)
try
    % Remove initial conditions for a clean simulation
    simOut_val = sim(mdl, 'StopTime', num2str(T_val_total), 'SaveOutput', 'on');
    
    if isfield(simOut_val, 'tau_J')
         y_real = simOut_val.tau_J.Data;
         t_real = simOut_val.tau_J.Time;
    else
         logs = simOut_val.logsout;
         y_real = logs.get('tau_J').Values.Data;
         t_real = logs.get('tau_J').Values.Time;
    end
    
    % Interpolate for comparison (not strictly required with timeseries)
    y_real_interp = interp1(t_real, y_real, t_val, 'linear', 'extrap');
    
catch ME
    disp(ME.message);
    warning('Validation simulation failed');
    y_real_interp = zeros(size(t_val));
end

% 6. Simulated model response
[y_est, ~] = lsim(sys_est, tau_jd_val_vec, t_val);

% 7. FIT computation (Normalized Root Mean Square Error)
idx_eval = t_val > 0.05;
y_ref_eval = y_real_interp(idx_eval);
y_est_eval = y_est(idx_eval);

num = norm(y_ref_eval - y_est_eval);
den = norm(y_ref_eval - mean(y_ref_eval));

% Protect against near-zero denominator
if den < 1e-6
    FIT_val = 0;
else
    FIT_val = 100 * (1 - num/den);
end

fprintf('\n>>> REST-TO-REST VALIDATION <<<\n');
fprintf('Trajectory FIT: %.2f %%\n', FIT_val);

% 8. Plot results
figure('Name', 'Rest-to-Rest Validation', 'Color', 'w');

% Subplot 1: Desired Link Position
subplot(3,1,1);
plot(t_val, q_val, 'k', 'LineWidth', 1.5); grid on;
ylabel('Position [rad]'); 
title('Desired Trajectory (5th-order Polynomial)');

% Subplot 2: Feedforward Torque Input
subplot(3,1,2);
plot(t_val, tau_jd_val_vec, 'k--'); grid on;
ylabel('Torque Input [Nm]');
title('Computed Feedforward Torque');

% Subplot 3: Real vs Estimated Output Torque
subplot(3,1,3);
plot(t_val, y_real_interp, 'b', 'LineWidth', 1.5); hold on;
plot(t_val, y_est, 'r--', 'LineWidth', 1.5);
ylabel('Output Torque [Nm]'); xlabel('Time [s]');
legend('Simulink (Real System)', 'Estimated Model');
grid on;
title(sprintf('Model Comparison – FIT: %.2f %%', FIT_val));

% Parameter comparison table
fprintf('\n--- PARAMETER COMPARISON ---\n');
fprintf('Kp_theta: Real=%.1f | Estimated=%.1f\n', K_p_theta, Kp_theta_est);
fprintf('Kd_theta: Real=%.1f | Estimated=%.1f\n', K_d_theta, Kd_theta_est);
fprintf('Kp_tau:   Real=%.1f | Estimated=%.1f\n', K_p_tau,   Kp_tau_est);
fprintf('Kd_tau:   Real=%.1f | Estimated=%.1f\n', K_d_tau,   Kd_tau_est);
disp('---------------------------');

disp('Estimated Transfer Function from the experimental data:');
disp(sys_est);
