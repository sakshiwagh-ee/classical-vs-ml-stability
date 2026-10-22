function ClassicalStabilityPredictionDemo()
% ClassicalStabilityPredictionDemo.m
% Save this file exactly as "ClassicalStabilityPredictionDemo.m" (must start with a letter).
% Run by typing: ClassicalStabilityPredictionDemo
% Requires Control System Toolbox for tf, feedback, margin, nyquist, bode
%
% This version is a function (file name must match the function name) to avoid
% MATLAB filename errors on startup.

clear; close all; clc;

%%\ ---------- USER INPUT (edit these) ----------
% Plant: second-order example: G(s) = Kp * wn^2 / (s^2 + 2*zeta*wn*s + wn^2)
wn   = 5;        % natural frequency (rad/s)
zeta = 0.35;     % damping ratio
Kp   = 1;        % plant DC gain

% Controller choice: uncomment one controller definition below or define your own
% 1) Proportional controller
Kc = 2.0;
C = tf(Kc);

% % 2) PI controller example (uncomment to use)
% Kp_c = 1.5; Ti = 0.5;
% C = tf([Kp_c*Ti Kp_c],[Ti 0]);

% % 3) PID controller example (uncomment to use)
% Kp_c = 1.0; Ki = 2.0; Kd = 0.01; N = 20; % N is derivative filter
% C = tf([Kd*N Kp_c*Kd*N+Kp_c Kp_c*Ki],[N 0]);

% Frequency grid for plots
w = logspace(-1, 2.5, 1000); % rad/s

%%\ ---------- BUILD PLANT AND CLOSED-LOOP ----------
G = tf(Kp*wn^2, [1 2*zeta*wn wn^2]);   % second-order plant
L = series(C, G);                       % loop transfer function (controller * plant)
T = feedback(L, 1);                     % closed-loop transfer function (unity feedback)
S = feedback(1, L);                     % sensitivity function

%%\ ---------- ROUTH-HURWITZ (closed-loop characteristic polynomial) ----------
% Closed-loop characteristic polynomial denominator (as coefficients)
[~, denT] = tfdata(T, 'v');  % den of closed-loop transfer function
denT = denT(:).';
[routhTable, routhStable] = routh(denT);

%%\ ---------- NYQUIST & MARGINS ----------
[GM, PM, Wcg, Wcp] = margin(L); % GM (linear), PM (deg), Wcg (rad/s), Wcp (rad/s)

% Convert GM to dB (if finite)
if isfinite(GM) && GM>0
    GM_dB = 20*log10(GM);
else
    GM_dB = sign(GM)*Inf; % preserve Inf or -Inf
end

%%\ ---------- BODE STABILITY METRICS ----------
[mag, phase, wout] = bode(L, w);
mag = squeeze(mag);
phase = squeeze(phase); % in degrees

if isempty(Wcp)
    PM_calc = NaN;
else
    PM_calc = PM; % degrees
end

%%\ ---------- UNIFIED STABILITY INDEX (USI) ----------
% Term 1: Routh (binary)
t1 = double(routhStable);

% Term 2: gain margin normalized: GM_dB mapped so that >= 12 dB => 1, <= 0 dB => 0
if isfinite(GM_dB)
    t2 = (GM_dB - 0) / (12 - 0);
    t2 = max(0, min(1, t2));
else
    t2 = 1 * sign(GM_dB);
    if t2 < 0; t2 = 0; end
end

% Term 3: phase margin normalized: PM >= 45 deg => 1, <= 0 deg => 0
if ~isnan(PM_calc)
    t3 = PM_calc / 45;
    t3 = max(0, min(1, t3));
else
    t3 = 0;
end

% Weights (sum to 1)
w1 = 0.4; w2 = 0.3; w3 = 0.3;
USI = w1*t1 + w2*t2 + w3*t3;

%%\ ---------- REPORT RESULTS ----------
fprintf('--- Classical Stability Criteria Report ---\n');
fprintf('Plant: wn = %.3g rad/s, zeta = %.3g, Kp = %.3g\n', wn, zeta, Kp);
fprintf('Controller type: %s\n', class(C));
fprintf('Routh-Hurwitz: Closed-loop %s.\n', tern(routhStable,'stable','unstable'));
fprintf('Gain Margin: %.3g (%.2f dB) at w = %.3g rad/s\n', GM, GM_dB, Wcg);
fprintf('Phase Margin: %.3g deg at w = %.3g rad/s\n', PM, Wcp);
fprintf('Unified Stability Index (USI): %.3f (0..1, higher = more robust)\n', USI);
fprintf('  Components: Routh=%g, GM_norm=%.3f, PM_norm=%.3f\n', t1, t2, t3);
fprintf('--------------------------------------------\n');

%%\ ---------- PLOTS ----------
figure('Name','Bode Plot (Loop Transfer L)','NumberTitle','off');
margin(L); grid on;
title('Bode & Margins of Loop Transfer Function L(s)');

figure('Name','Nyquist Plot (Loop Transfer L)','NumberTitle','off');
nyquist(L, {0.01, 200}); hold on;
plot(-1, 0, 'r+', 'MarkerSize', 10, 'LineWidth', 1.5); % -1 point
title('Nyquist Plot of L(s) (red + marks -1 point)');
grid on; hold off;

figure('Name','Closed-loop Step Response T(s)','NumberTitle','off');
tstep = 0:0.01:5;
step(T, tstep); grid on; title('Closed-loop Step Response T(s)');

figure('Name','Sensitivity Functions','NumberTitle','off');
subplot(2,1,1);
bode(S, w); grid on; title('Sensitivity S(s)');
subplot(2,1,2);
bode(T, w); grid on; title('Complementary Sensitivity T(s)');

%%\ ---------- DISPLAY Routh TABLE ----------
disp('Routh-Hurwitz table:');
disp(routhTable);

end % function ClassicalStabilityPredictionDemo

%%\ ---------------- Local helper functions ----------------
function s = tern(cond, a, b)
    if cond, s = a; else s = b; end
end

function [R, stable] = routh(coeffs)
% routh - Build Routh-Hurwitz table and return table and stability flag
% coeffs - polynomial coefficients in descending powers [a_n ... a0]

    a = coeffs(:).';
    n = length(a)-1; % degree
    if all(a==0)
        error('Polynomial coefficients cannot be all zero.');
    end

    rows = n+1;
    cols = ceil((n+1)/2);

    R = zeros(rows, cols);

    % Fill first two rows
    R(1,1:cols) = a(1:2:end);
    if length(a) > 1
        R(2,1:cols) = a(2:2:end);
    else
        R(2,1:cols) = 0;
    end

    % Fill remaining rows
    for i = 3:rows
        for j = 1:cols-1
            top = R(i-1,1);
            if abs(top) < eps
                top = eps; % avoid division by zero
            end
            R(i,j) = ((R(i-1,1)*R(i-2,j+1)) - (R(i-2,1)*R(i-1,j+1))) / top;
        end
    end

    % Handle special case: an entire row of zeros
    for i = 1:rows
        if all(abs(R(i,:)) < 1e-12)
            order = rows - i;
            aux = R(i-1,1:cols);
            p = zeros(1, order+1);
            idx = 1;
            for k = order:-1:0
                p(idx) = aux(1 + (order - k));
                idx = idx + 1;
            end
            dp = polyder(p);
            fill = dp;
            fill = [fill zeros(1, cols - length(fill))];
            R(i,1:length(fill)) = fill(1:cols);
        end
    end

    firstCol = R(:,1);
    firstCol(abs(firstCol) < 1e-12) = 0;

    signs = sign(firstCol);
    signs(signs==0) = 1; % treat zeros as positive
    signChanges = sum(abs(diff(signs))>0);
    stable = (signChanges == 0) && (firstCol(1) > 0);
end

%[appendix]{"version":"1.0"}
%---
