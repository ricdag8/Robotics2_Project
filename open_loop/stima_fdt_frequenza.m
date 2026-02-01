clc; clear; 
clear; clc; close all;
%% ========================================================================
%   SYSTEM PARAMETERS
% ========================================================================
Mm = 1;  
M  = 2;   
K  = 1000;  
rng(42);


Ts = 0.001;      
        


%we want to estimate the TF of the open loop system. for this reason we are
%going to use a systematic approach based on the system's frequency
%response and then we are going to compute the TF by means of the LS
%algorithm. This approach is completely black box and does not rely on any
%knowledge of the model. 



mdl = 'modello_giunto';
load_system(mdl); 

% Noise parameters
STD_NOISE_THETA = 0.001;  
STD_NOISE_TAU   = 0.01; 


%% ========================================================================
%   SIMULATION 
% ======================================================================== 

f_start_hz = 1;
f_end_hz   = 7;

w_start = f_start_hz * 2 * pi; 
w_end   = f_end_hz   * 2 * pi; 

frequncy_to_test = logspace(log10(w_start), log10(w_end), 11);

results = struct(); 

for k = 1:length(frequncy_to_test)
   
   
    w_test = frequncy_to_test(k); %we take the current frequency
    A = 10.0; 
    T_period = 2*pi / w_test;
    
    % we are going to use different periods depending on the frequency 

    if w_test < 20.0

        Num_per = 10;

    else
         Num_per = 100; 
    end
    
    Num_Periodi_Transitorio = 0;

    T_transitorio = Num_Periodi_Transitorio * T_period;
    T_regime      = Num_per * T_period;
    Tf = T_transitorio + T_regime;
    
    % input signal we feed to the model
    dt_sim = 0.001; 
    t_vec = 0:dt_sim:Tf;
    u_values = A * cos(w_test * t_vec);
    input_signal = timeseries(u_values, t_vec);
    
    % input to simulink
    assignin('base', 'u_in', input_signal);
    assignin('base', 'Tf', Tf);
    
    fprintf('simulation %d/%d: Freq=%.2f rad/s | noise on\n', k, length(frequncy_to_test), w_test);
        
    % simulation
    simOut = sim(mdl, 'StopTime', num2str(Tf), 'SaveOutput', 'on');
    try
        logs = simOut.logsout;
        try sig_th = logs.getElement('theta'); catch, sig_th = logs.get(1); end
        try sig_tau = logs.getElement('tau_J'); catch, try sig_tau = logs.getElement('tau_j'); catch, sig_tau = logs.get(2); end; end
        
        t_rec_raw = sig_th.Values.Time;
        theta_raw = sig_th.Values.Data; 
        tau_raw   = sig_tau.Values.Data;  
        
    catch
        theta_ws = evalin('base', 'theta');
        if isa(theta_ws, 'timeseries'), t_rec_raw=theta_ws.Time; theta_raw=theta_ws.Data;
        else, t_rec_raw=theta_ws(:,1); theta_raw=theta_ws(:,2); end
        
        if evalin('base', 'exist(''tau_J'',''var'')'), tau_ws = evalin('base', 'tau_J'); else, tau_ws = evalin('base', 'tau_j'); end
        if isa(tau_ws, 'timeseries'), tau_raw=tau_ws.Data; else, tau_raw=tau_ws(:,2); end
    end

    if isrow(u_values), u_values = u_values'; end
    if isrow(theta_raw), theta_raw = theta_raw'; end
    if isrow(tau_raw), tau_raw = tau_raw'; end
    if isrow(t_rec_raw), t_rec_raw = t_rec_raw'; end

    % We compute the minimum common length to avoid dimension errors.
    len_in = length(u_values);
    len_out = length(theta_raw);
    n_common = min(len_in, len_out);
    

    input_rec   = u_values(1:n_common);
    theta_clean = theta_raw(1:n_common);
    tau_clean   = tau_raw(1:n_common);
    t_rec       = t_rec_raw(1:n_common);
  
    % we add white gaussian noise to the data
    theta_noisy = theta_clean + STD_NOISE_THETA * randn(size(theta_clean));
    tau_noisy   = tau_clean   + STD_NOISE_TAU   * randn(size(tau_clean));



% --- SPECTRAL ANALYSIS (FFT), ONLY FOR PLOTS PURPOSES ---

    figure('Name', ['Analisi FFT - Freq Input: ' num2str(w_test) ' rad/s'], 'Color', 'w');
    
 
    y_signal = theta_noisy;      
    t_signal = t_rec;            % 
    y_signal = y_signal - mean(y_signal);
    
 
    Ts_eff = mean(diff(t_signal)); % Average sampling time
    Fs = 1 / Ts_eff;               % Sampling frequency [Hz]
    L = length(y_signal);          % Signal length
    
    % Compute FFT
    Y_fft = fft(y_signal);
    
    %  Compute Amplitude Spectrum (P1)
    P2 = abs(Y_fft/L);             % Normalized two-sided spectrum
    P1 = P2(1:floor(L/2)+1);       % Take only the positive half
    P1(2:end-1) = 2*P1(2:end-1);   % Multiply by 2 (energy conservation)

    
    f_hz = Fs * (0:(L/2)) / L;     % Frequency in hz
    w_rad = 2 * pi * f_hz;         % Frequency in rad/s
    
    % --- PLOT ---
    subplot(2,1,1);
    plot(t_signal, y_signal, 'b'); 
    grid on; title('Segnale nel Tempo (Theta con Rumore)');
    xlabel('Tempo [s]'); ylabel('Ampiezza [rad]');
    
    subplot(2,1,2);
    plot(w_rad, P1, 'r', 'LineWidth', 1.5); hold on;
    
  
    xline(w_test, 'g--', 'LineWidth', 1.5, 'Label', 'Input Freq');
    
    grid on; title('Spettro di Ampiezza (FFT)');
    xlabel('Frequenza [rad/s]'); ylabel('Magnitude |P1(w)|');
    
   
    xlim([0, w_test * 5]);

   
    results(k).freq = w_test;
    results(k).T_transitorio = T_transitorio; 
    results(k).time = t_rec;
    results(k).input = input_rec;
    
    results(k).output_theta = theta_noisy; 
    results(k).output_tau = tau_noisy;
    
    

    
end





%% =========================================================================
%  ESTIMATE OF THE FREQUENCY RESPONSE (Sperimental Bode)
% =========================================================================

N_test = length(results);

%
Freqs_rad = zeros(N_test, 1);  
H_complex = zeros(N_test, 1);  

disp('--- ELABORAZIONE DATI TEMPORALI ---');

for k = 1:N_test
    
    % retrieve all the data
    t_raw = results(k).time;
    u_raw = results(k).input;
    
    % we are going to use theta as the output to estimate the coefficients
    y_raw = results(k).output_theta; 
    
    w_val = results(k).freq;
    T_cut = results(k).T_transitorio;

    
    % we only take the values after the transient 
    mask_steady = t_raw > T_cut;
    
    u_steady = u_raw(mask_steady);
    y_steady = y_raw(mask_steady);
    t_steady = t_raw(mask_steady);

    % we take the number of samples 
    N_samples = length(y_steady);
    
    % we then compute the FFT on both the input and output signals
    FFT_u = fft(u_steady);
    
    FFT_y = fft(y_steady);
    
    % %for each experiment/input we have done, we basically look for the
    %maximum value in the spectrum in order to isolate the component we
    %are looking for, namely the one used in order to stimulate the system
    [~, idx_peak] = max(abs(FFT_u(1:floor(N_samples/2)+1)));
    
    % we exctract the values of the peak
    U_phasor = FFT_u(idx_peak);
    Y_phasor = FFT_y(idx_peak);
    
    % and finally we compute the estimate of our frequency response 
    H_val = Y_phasor / U_phasor;
    
    % save everythign
    Freqs_rad(k) = w_val;
    H_complex(k) = H_val;




   

end


% we then compute dB and phase for the plot
Mag_dB = 20*log10(abs(H_complex));
Phase_deg = rad2deg(unwrap(angle(H_complex)));

% =========================================================================
% PLOT OF THE EXPERIMENTAL BODE DIAGRAMS    
% =========================================================================
figure('Name', 'Exmperimental Bode', 'Color', 'w');

subplot(2,1,1);
semilogx(Freqs_rad, Mag_dB, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
grid on; 
ylabel('Magnitudo [dB]'); 
title('Estimate of the Frequency Response computed from the data');

subplot(2,1,2);
semilogx(Freqs_rad, Phase_deg, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
grid on; 
ylabel('Phase [deg]'); 
xlabel('Frequency [rad/s]');


disp('FDT Identification');

% one could either use the dynamic model of the system in order to choose
% the order of the model to be estimated, which clearly guides us to choose
% na = 4 and nb = 2. also, by looking at the Bode plots of the frequency
% response, we see both a resonance and anti-resonance component, which
% suggest to take 2 zeros (since the the resonance phenomena is due to an
% imaginary pair of zeros and also a 4 poles, 2 in the origin and 2 complex
% conjugate.


na = 4; % den
nb = 2; % num

disp(['identification with: ' num2str(na) ' Poles, ' num2str(nb) ' Zeros']);


w = Freqs_rad;
H = H_complex;
N = length(w);

% Initialize matrices for the system Ax = B
% Number of unknowns: (nb+1) 'b' coefficients + na 'a' coefficients
num_vars = (nb + 1) + na; 
A_mat = zeros(2 * N, num_vars); % 2 rows for each frequency (Real + Imag)
B_vec = zeros(2 * N, 1);

% Row-by-row construction of the system
for k = 1:N
    s_val = 1j * w(k); % Laplace variable s = jw
    h_val = H(k);      % Measured value
    
    % --- CONSTRUCTION OF THE REGRESSOR ROW ---
    % We want to write: b0 + b1*s + ... – a1*h*s – a2*h*s^2 = h
    
    % Part for the b coefficients (Numerator): [1, s, s^2 ...]
    reg_num = zeros(1, nb+1);
    for i = 0:nb
        reg_num(i+1) = s_val^i;
    end
    
    % Part for the a coefficients (Denominator): [-h*s, -h*s^2 ...]
    reg_den = zeros(1, na);
    for i = 1:na
        reg_den(i) = -h_val * (s_val^i);
    end
    
    % Full complex row
    row_complex = [reg_num, reg_den];
    rhs_complex = h_val; % Right-hand side (a0 * h = 1 * h)
    
    % Real part row
    A_mat(2*k-1, :) = real(row_complex);
    B_vec(2*k-1)    = real(rhs_complex);
    
    % Imaginary part row
    A_mat(2*k, :)   = imag(row_complex);
    B_vec(2*k)      = imag(rhs_complex);
end

%% --- ADDITION: WEIGHTED LEAST SQUARES (WLS) – MANUAL ASSIGNMENT ---
disp('Applying manual weights to selected frequency regions...');

% 1. Create base weight vector (all equal to 1)
W = ones(N, 1);

% 2. Define weight factor (how much more important these zones are)
W_factor = 100;

% Increase weight around first region of interest
W(Freqs_rad > 5.5 & Freqs_rad < 6.5) = W_factor;

% Increase weight around second region of interest
W(Freqs_rad >= 20 & Freqs_rad <= 25) = W_factor;

figure('Name', 'Weights distribution', 'Color', 'w');
stem(Freqs_rad, W, 'filled');
title('Weight assignment');
xlabel('Frequency [rad/s]'); ylabel('Weight');
grid on;

% Expand weight vector to match real/imag rows
W_expanded = zeros(2*N, 1);
for k = 1:N
    W_expanded(2*k-1) = W(k);
    W_expanded(2*k)   = W(k);
end

% 5. Solve the WEIGHTED system
A_weighted = A_mat .* W_expanded;
B_weighted = B_vec .* W_expanded;

Theta = pinv(A_weighted) * B_weighted;

disp('Manual WLS estimation completed.');








% --- COEFFICIENT EXTRACTION ---
% The first (nb+1) elements are the numerator (b)
b_raw = Theta(1 : nb+1);       

% The remaining elements are the denominator (a), recalling that a0 = 1
a_raw = [1; Theta(nb+2 : end)]; 




num_est = flip(b_raw).'; 
den_est = flip(a_raw).';

% Creation of the Transfer Function
sys_est = tf(num_est, den_est);

disp(sys_est)


%%  POLE-ZERO ANALYSIS
disp('--- MODAL ANALYSIS ---');
p = pole(sys_est);
z = zero(sys_est);

disp('System poles (damping ratio and natural frequency):');
damp(sys_est);

disp('System zeros:');
disp(z);



%% =========================================================================
%  EXTRACTION OF PHYSICAL PARAMETERS (Inverse Dynamics)
% =========================================================================
disp('--- EXTRACTION OF PHYSICAL PARAMETERS (Coefficient Method) ---');

% 1. Extract Coefficient Vectors
[num, den] = tfdata(sys_est, 'v');

% 2. Normalization
% Ensure that the denominator polynomial starts with 1 (monic)
% H(s) = (num/den(1)) / (den/den(1))
norm_factor = den(1);
num_n = num / norm_factor;
den_n = den / norm_factor;

% 3. Coefficient Mapping
% For a 4th-order system, MATLAB vectors are structured as:
% Index:     1     2     3     4     5
% Powers:   s^4   s^3   s^2   s^1   s^0

% Coefficient of s^2 in the numerator (Corresponds to 1/Mm)
n2 = num_n(3);

% Constant term in the numerator (Corresponds to K / (Mm*M))
n0 = num_n(5);

% Coefficient of s^2 in the denominator (Corresponds to K(Mm+M)/(Mm*M))
d2 = den_n(3);

% 4. Computation of Physical Parameters
% From theory: n2 = 1 / Mm
Mm_est = 1 / n2;

% From theory: omega_z^2 = n0 / n2
% Derivation: (K/(Mm*M)) / (1/Mm) = K/M
omega_z_sq = n0 / n2;

% From theory: omega_p^2 = d2
% The s^2 coefficient in the denominator is the squared resonance frequency
omega_p_sq = d2;

% Compute K
% We know that: d2 = K/M + K/Mm  (which equals omega_p^2)
% And: omega_z_sq = K/M
% Therefore: K/Mm = d2 - omega_z_sq
K_est = Mm_est * (omega_p_sq - omega_z_sq);

% Compute M (Link Inertia)
% We know that: M = K / omega_z^2
M_est = K_est / omega_z_sq;

% 5. Display Results
disp('---------------------------------------------');
fprintf('PARAMETERS EXTRACTED FROM COEFFICIENTS:\n');
fprintf('Motor Inertia (Mm): %.6f [kg m^2]\n', Mm_est);
fprintf('Link Inertia (M):   %.6f [kg m^2]\n', M_est);
fprintf('Joint Stiffness (K): %.4f [Nm/rad]\n', K_est);
disp('---------------------------------------------');


fprintf('Resonance check (sqrt(d2)):           %.2f rad/s\n', sqrt(d2));
fprintf('Anti-resonance check (sqrt(n0/n2)):   %.2f rad/s\n', sqrt(omega_z_sq));

%% =========================================================================
% VISUAL ANALYSIS OF THE TRANSIENT (Red = Discarded, Blue = Used)
% =========================================================================

disp('Identified Transfer Function (Manual):');

% --- VALIDATION (GRAPHICAL COMPARISON) ---
w_verify = logspace(log10(min(Freqs_rad)), log10(max(Freqs_rad)), 100);

% 1. Standard Bode Computation
[mag_est, phase_est] = bode(sys_est, w_verify);
mag_est = squeeze(mag_est); 
phase_est = squeeze(phase_est);

% --- FINAL SOLUTION: ADAPTIVE PHASE ALIGNMENT ---
% 1. Create a "reference curve" by interpolating the experimental data
%    onto the dense frequency grid (pchip or linear interpolation)
phase_target = interp1(Freqs_rad, Phase_deg, w_verify, 'linear', 'extrap');

% 2. Force each model phase point to the closest 360° branch of the target
phase_est_aligned = zeros(size(phase_est));

for i = 1:length(phase_est)
    % Compute the difference between desired phase (data) and current phase (model)
    diff = phase_target(i) - phase_est(i);
    
    % Find how many full 360° rotations to apply (nearest integer)
    k_shifts = round(diff / 360);
    
    % Apply the shift
    phase_est_aligned(i) = phase_est(i) + k_shifts * 360;
end

% --- OVERLAID PLOT ---
figure('Name', 'Identification Validation (Aligned)', 'Color', 'w');

% Magnitude
subplot(2,1,1);
semilogx(Freqs_rad, 20*log10(abs(H_complex)), 'bo', 'LineWidth', 1.5, ...
    'DisplayName', 'Experimental Data'); 
hold on;
semilogx(w_verify, 20*log10(mag_est), 'r-', 'LineWidth', 2, ...
    'DisplayName', 'Identified Model (LS)');
grid on; legend('Location', 'best');
ylabel('Magnitude [dB]'); 
title(['Model Fit (Order ' num2str(na) '/' num2str(nb) ')']);

% Phase (Now perfectly aligned)
subplot(2,1,2);
semilogx(Freqs_rad, Phase_deg, 'bo', 'LineWidth', 1.5); hold on;
semilogx(w_verify, phase_est_aligned, 'r-', 'LineWidth', 2);
grid on; 
ylabel('Phase [deg]'); 
xlabel('Frequency [rad/s]');

% 3. Additional Overlaid Plot
figure('Name', 'Identification Validation', 'Color', 'w');

% Magnitude
subplot(2,1,1);
semilogx(Freqs_rad, 20*log10(abs(H_complex)), 'bo', 'LineWidth', 1.5, ...
    'DisplayName', 'Experimental Data');
hold on;
semilogx(w_verify, 20*log10(mag_est), 'r-', 'LineWidth', 2, ...
    'DisplayName', 'Identified Model (LS)');
grid on; legend('Location', 'best');
ylabel('Magnitude [dB]');
title(['Model Fit (Order ' num2str(na) '/' num2str(nb) ')']);

% Phase
subplot(2,1,2);
semilogx(Freqs_rad, rad2deg(unwrap(angle(H_complex))), 'bo', 'LineWidth', 1.5);
hold on;

% Unwrap handling for model
if abs(mean(phase_est) - mean(rad2deg(unwrap(angle(H_complex))))) > 180
    phase_est = phase_est + 360;
end

semilogx(w_verify, phase_est, 'r-', 'LineWidth', 2);
grid on; 
ylabel('Phase [deg]'); 
xlabel('Frequency [rad/s]');










%% =========================================================================
% VISUALIZZAZIONE RISPOSTA NEL TEMPO PER 3 FREQUENZE
% =========================================================================

% Frequenze test
target_freqs = [22.36, 30, 38];
colors = {'b', 'g', 'r'}; 


figure('Name', 'Input and Output signals', 'Color', 'w');

for i = 1:length(target_freqs)
    w_target = target_freqs(i);
    

    [~, idx_sim] = min(abs([results.freq] - w_target));
    

    t_plot = results(idx_sim).time;
    u_plot = results(idx_sim).input;
    y_plot = results(idx_sim).output_theta;

    T_p = 2*pi / results(idx_sim).freq;
    t_start_plot = t_plot(end) - 10*T_p;
    mask_plot = t_plot > t_start_plot;

    subplot(3, 1, i);
   
    yyaxis left
    plot(t_plot(mask_plot), y_plot(mask_plot), 'Color', colors{i}, 'LineWidth', 1.2);
    ylabel('Output \theta [rad]');
    
    yyaxis right
    plot(t_plot(mask_plot), u_plot(mask_plot), 'k--', 'LineWidth', 1.0);
    ylabel('Input Torque [Nm]');
    

    grid on;
    title(sprintf('%.1f rad/s', results(idx_sim).freq));
    xlabel('Time [s]');
    legend('Output \theta (Measured)', 'Input ', 'Location', 'best');
end

% =========================================================================






    %% --- CODICE DI DEBUG: PLOT FFT ---

if k == N_test  
    
    figure('Name', ['Spectral analysis - Freq: ' num2str(w_val) ' rad/s'], 'Color', 'w');
    
  
    Fs = 1 / 0.001; 
    df = Fs / N_samples;               
    f_axis = (0:N_samples-1) * df;     
    

    f_target_hz = w_val / (2*pi);

    subplot(2,1,1);

    mag_u = abs(FFT_u) / (N_samples/2); 
    
    plot(f_axis, mag_u, 'b-', 'LineWidth', 1); hold on;
   
    plot(f_axis(idx_peak), mag_u(idx_peak), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    
    xline(f_target_hz, 'k--', 'Freq. Target');
    title(['Input Spectra (Torque) - Target: ' num2str(f_target_hz, '%.2f') ' Hz']);
    xlabel('Frequency [Hz]'); ylabel('Amplitude [Nm]');
    grid on;
    
    xlim([0, f_target_hz * 5]); 

    subplot(2,1,2);
    mag_y = abs(FFT_y) / (N_samples/2);
    
    plot(f_axis, mag_y, 'r-', 'LineWidth', 1); hold on;

    plot(f_axis(idx_peak), mag_y(idx_peak), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    
    xline(f_target_hz, 'k--', 'Freq. Target');
    title('Output Spectra (Theta position with noise)');
    xlabel('Frequency [Hz]'); ylabel('Amplitude [rad]');
    grid on;
    xlim([0, f_target_hz * 5]);
end
    