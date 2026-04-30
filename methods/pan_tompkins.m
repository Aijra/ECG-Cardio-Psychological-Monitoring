clc; clear; close all;

%% ========= 1. DEMOGRAPHICS (as provided in the original PDF) =========
demo = containers.Map('KeyType','double','ValueType','any');

demo(100)=struct('age',69,'gender','M'); demo(101)=struct('age',75,'gender','F');
demo(102)=struct('age',84,'gender','F'); demo(103)=struct('age',NaN,'gender','M');
demo(104)=struct('age',66,'gender','F'); demo(105)=struct('age',73,'gender','F');
demo(106)=struct('age',24,'gender','F'); demo(107)=struct('age',63,'gender','M');
demo(108)=struct('age',76,'gender','F'); demo(109)=struct('age',77,'gender','M');
demo(111)=struct('age',50,'gender','F'); demo(112)=struct('age',39,'gender','M');
demo(113)=struct('age',74,'gender','F'); demo(114)=struct('age',69,'gender','F');
demo(115)=struct('age',72,'gender','F'); demo(116)=struct('age',66,'gender','M');
demo(117)=struct('age',66,'gender','M'); demo(118)=struct('age',66,'gender','M');
demo(119)=struct('age',63,'gender','M'); demo(121)=struct('age',51,'gender','M');
demo(122)=struct('age',63,'gender','M'); demo(123)=struct('age',59,'gender','M');
demo(124)=struct('age',65,'gender','M');

demo(200)=struct('age',64,'gender','M'); demo(201)=struct('age',64,'gender','M');
demo(202)=struct('age',64,'gender','M'); demo(203)=struct('age',74,'gender','F');
demo(205)=struct('age',66,'gender','M'); demo(207)=struct('age',61,'gender','M');
demo(208)=struct('age',69,'gender','M'); demo(209)=struct('age',65,'gender','F');
demo(210)=struct('age',64,'gender','M'); demo(212)=struct('age',67,'gender','M');
demo(213)=struct('age',66,'gender','F'); demo(214)=struct('age',53,'gender','M');
demo(215)=struct('age',55,'gender','M'); demo(217)=struct('age',66,'gender','M');
demo(219)=struct('age',53,'gender','F'); demo(220)=struct('age',64,'gender','M');
demo(221)=struct('age',70,'gender','M'); demo(222)=struct('age',67,'gender','M');
demo(223)=struct('age',59,'gender','F'); demo(228)=struct('age',22,'gender','F');
demo(230)=struct('age',71,'gender','F'); demo(231)=struct('age',66,'gender','F');
demo(232)=struct('age',76,'gender','F'); demo(233)=struct('age',70,'gender','F');
demo(234)=struct('age',55,'gender','M');

%% ========= 2. SETTINGS =========
% List of records to process (skip 102 if missing, as in original)
records = [100:101,103:124,200:234];
Fs = 360;   % sampling frequency (Hz)

% Folder where .dat and .csv files are located
data_folder = pwd;   % change to your folder if needed, e.g. 'C:/MyData'

% Storage for overall evaluation
total_TP = 0;
total_FP = 0;
total_FN = 0;
per_record_F1 = [];

% Storage for HRV classification
psych_labels = {};
cardiac_labels = {};

%% ========= 3. MAIN LOOP OVER RECORDS =========
for i = 1:length(records)
    rec = records(i);
    fprintf('\n========================================\n');
    fprintf('Processing Record %d...\n', rec);
    
    % ----- Load ECG from .dat file (MIT format) -----
    dat_file = fullfile(data_folder, sprintf('%03d.dat', rec));
    fid = fopen(dat_file, 'r');
    if fid == -1
        warning('Missing .dat file: %s', dat_file);
        continue;
    end
    raw = fread(fid, [3 inf], 'uint8');
    fclose(fid);
    
    % Reconstruct first channel (usually MLII) – 12-bit signed
    ecg = double(bitshift(raw(2,:),4) + bitshift(raw(1,:),-4));
    ecg(ecg > 2047) = ecg(ecg > 2047) - 4096;
    
    % ----- Load annotations from .csv file -----
    csv_file = fullfile(data_folder, sprintf('%03d_annotations_1.csv', rec));
    if ~exist(csv_file, 'file')
        warning('Missing .csv file: %s', csv_file);
        continue;
    end
    ann = readtable(csv_file);
    
    % Keep only beat types that are considered "normal" or "true" QRS
    % According to MIT-BIH: N (normal), L (LBBB), R (RBBB), A (APC), V (PVC)
    valid_symbols = {'N','L','R','A','V'};
    mask = ismember(ann.annotation_symbol, valid_symbols);
    true_peaks = double(ann.index(mask));
    
    % Sanity check
    true_peaks(true_peaks > length(ecg)) = [];
    if isempty(true_peaks)
        warning('No valid annotations for record %d', rec);
        continue;
    end
    
    % ----- Improved QRS detection -----
    detected_peaks = qrs_detection_pan_tompkins(ecg, Fs);
    
    % ----- Evaluation (match detected to true) -----
    tol = round(0.1 * Fs);   % 100 ms tolerance
    matched = false(size(true_peaks));
    TP = 0; FP = 0;
    
    for d = 1:length(detected_peaks)
        [min_dist, idx] = min(abs(true_peaks - detected_peaks(d)));
        if min_dist <= tol && ~matched(idx)
            TP = TP + 1;
            matched(idx) = true;
        else
            FP = FP + 1;
        end
    end
    FN = sum(~matched);
    
    % Update totals
    total_TP = total_TP + TP;
    total_FP = total_FP + FP;
    total_FN = total_FN + FN;
    
    % Compute metrics for this record
    Se = TP / (TP + FN + eps);
    P  = TP / (TP + FP + eps);
    F1 = 2 * Se * P / (Se + P + eps);
    per_record_F1(end+1) = F1;
    
    fprintf('  TP=%d, FP=%d, FN=%d\n', TP, FP, FN);
    fprintf('  Sensitivity (Se) = %.3f, Positive Predictivity (P) = %.3f, F1 = %.3f\n', Se, P, F1);
    
    % % ----- Plot a 10-second segment if F1 > 0.95 (optional) -----
    % if F1 > 0.95 && ~isempty(detected_peaks)
    %     plot_start = max(1, detected_peaks(1) - round(1*Fs));
    %     plot_end   = min(length(ecg), plot_start + round(10*Fs) - 1);
    %     t = (plot_start:plot_end) / Fs;
    %     figure('Name', sprintf('Record %d (F1=%.2f%%)', rec, F1*100));
    %     plot(t, ecg(plot_start:plot_end), 'b'); hold on;
    % 
    %     % Detected peaks in this window
    %     det_win = detected_peaks(detected_peaks >= plot_start & detected_peaks <= plot_end);
    %     plot(det_win/Fs, ecg(det_win), 'ro', 'MarkerFaceColor','r');
    % 
    %     % True peaks in this window
    %     true_win = true_peaks(true_peaks >= plot_start & true_peaks <= plot_end);
    %     plot(true_win/Fs, ecg(true_win), 'g*');
    % 
    %     xlabel('Time (s)'); ylabel('Amplitude');
    %     title(sprintf('Record %d: QRS detection (F1 = %.2f%%)', rec, F1*100));
    %     legend('ECG','Detected','True','Location','best');
    %     grid on;
    % end
    % 
    % ----- HRV analysis and classification -----
    if length(detected_peaks) > 5
        RR = diff(detected_peaks) / Fs;
        % Remove physiologically impossible RR intervals (0.3 to 1.5 s)
        RR = RR(RR > 0.3 & RR < 1.5);
        if length(RR) > 5
            % Mild median filtering to remove outliers
            RR = medfilt1(RR, 5);
            mean_HR = 60 / mean(RR);
            SDNN = std(RR);
            RMSSD = sqrt(mean(diff(RR).^2));
            
            % Age-based thresholds (derived from literature)
            age = demo(rec).age;
            if isnan(age), age = 50; end
            
            if age > 70
                RMSSD_low = 0.012; RMSSD_high = 0.035; SDNN_low = 0.020;
            elseif age > 50
                RMSSD_low = 0.015; RMSSD_high = 0.040; SDNN_low = 0.025;
            else
                RMSSD_low = 0.020; RMSSD_high = 0.050; SDNN_low = 0.030;
            end
            
            % Psychological state classification
            if mean_HR > 100 && RMSSD < RMSSD_low
                psych = 'ACUTE STRESS';
            elseif mean_HR > 90 && RMSSD < RMSSD_low
                psych = 'MILD STRESS';
            elseif RMSSD >= RMSSD_low && RMSSD <= RMSSD_high && mean_HR < 80
                psych = 'RELAXED';
            elseif RMSSD > RMSSD_high
                psych = 'RELAXED';
            else
                psych = 'NORMAL';
            end
            
            % Cardiac state classification
            if mean_HR > 100
                cardiac = 'TACHYCARDIA';
            elseif mean_HR < 55
                cardiac = 'BRADYCARDIA';
            elseif SDNN < SDNN_low
                cardiac = 'LOW HRV';
            elseif SDNN > 0.1
                cardiac = 'ARRHYTHMIA RISK';
            else
                cardiac = 'NORMAL';
            end
            
            psych_labels{end+1} = psych;
            cardiac_labels{end+1} = cardiac;
        end
    end
end

%% ========= 4. FINAL EVALUATION (overall) =========
overall_Se = total_TP / (total_TP + total_FN);
overall_P  = total_TP / (total_TP + total_FP);
overall_F1 = 2 * overall_Se * overall_P / (overall_Se + overall_P);

fprintf('\n================ FINAL RESULTS ================\n');
fprintf('Total TP: %d\n', total_TP);
fprintf('Total FP: %d\n', total_FP);
fprintf('Total FN: %d\n', total_FN);
fprintf('Overall Sensitivity (Se) = %.4f\n', overall_Se);
fprintf('Overall Positive Predictivity (P) = %.4f\n', overall_P);
fprintf('Overall F1-score = %.4f\n', overall_F1);

% Histogram of per-record F1 scores
figure;
histogram(per_record_F1, 10);
xlabel('F1-score'); ylabel('Number of records');
title('Distribution of QRS detection F1-score across MIT-BIH records');

%% ========= PSYCHOLOGICAL DISTRIBUTION =========
[unique_psych, ~, idx] = unique(psych_labels);
counts_psych = accumarray(idx, 1);

percent_psych = (counts_psych / sum(counts_psych)) * 100;

fprintf('\n===== PSYCHOLOGICAL STATE DISTRIBUTION =====\n');
for i = 1:length(unique_psych)
    fprintf('%s: %d (%.2f%%)\n', ...
        unique_psych{i}, counts_psych(i), percent_psych(i));
end

% ========= CARDIAC DISTRIBUTION =========
[unique_card, ~, idx2] = unique(cardiac_labels);
counts_card = accumarray(idx2, 1);

percent_card = (counts_card / sum(counts_card)) * 100;

fprintf('\n===== CARDIAC STATE DISTRIBUTION =====\n');
for i = 1:length(unique_card)
    fprintf('%s: %d (%.2f%%)\n', ...
        unique_card{i}, counts_card(i), percent_card(i));
end

%Psychological pie chart
figure;
p = pie(counts_psych);

title('Psychological State Distribution');

% Remove default text
delete(findobj(p,'Type','text'));

% Create legend labels
total = sum(counts_psych);
legend_labels = cell(length(unique_psych),1);

for i = 1:length(unique_psych)
    percent = (counts_psych(i)/total)*100;
    legend_labels{i} = sprintf('%s - %.1f%% (%d)', ...
        unique_psych{i}, percent, counts_psych(i));
end

legend(legend_labels, 'Location', 'eastoutside');

%Cardiac Pie Chart
figure;
p = pie(counts_card);

title('Cardiac Condition Distribution');

% Remove default text
delete(findobj(p,'Type','text'));

% Create legend labels
total = sum(counts_card);
legend_labels = cell(length(unique_card),1);

for i = 1:length(unique_card)
    percent = (counts_card(i)/total)*100;
    legend_labels{i} = sprintf('%s - %.1f%% (%d)', ...
        unique_card{i}, percent, counts_card(i));
end

legend(legend_labels, 'Location', 'eastoutside');

%% ========================================================================
%  FUNCTION: qrs_detection_pan_tompkins
%  A robust implementation of the Pan-Tompkins algorithm.
%  ========================================================================
function peaks = qrs_detection_pan_tompkins(ecg, Fs)
    % Remove baseline wander (0.5 Hz high-pass)
    [b_hp, a_hp] = butter(2, 0.5/(Fs/2), 'high');
    ecg_filt = filtfilt(b_hp, a_hp, double(ecg));
    
    % Bandpass filter 5-15 Hz
    [b_bp, a_bp] = butter(2, [5 15]/(Fs/2), 'bandpass');
    ecg_bp = filtfilt(b_bp, a_bp, ecg_filt);
    
    % Derivative
    deriv = diff(ecg_bp);
    deriv(end+1) = deriv(end);
    
    % Squaring
    squared = deriv .^ 2;
    
    % Moving window integration (150 ms)
    win_len = round(0.15 * Fs);
    mwi = movmean(squared, win_len);
    
    % Initial thresholds (first 2 seconds)
    init_len = min(2*Fs, length(mwi));
    SPKI = 0.25 * max(mwi(1:init_len));
    NPKI = 0.5 * mean(mwi(1:init_len));
    THR_SIG = NPKI + 0.25 * (SPKI - NPKI);
    THR_NOISE = 0.5 * THR_SIG;
    
    % Find candidate peaks (minimum distance 200 ms)
    [pks, locs] = findpeaks(mwi, 'MinPeakDistance', round(0.2*Fs));
    
    detected = [];
    last_peak = -inf;
    RR_list = [];
    RR_avg = 0;
    
    for i = 1:length(pks)
        pk = pks(i);
        loc = locs(i);
        
        % Refractory period
        if (loc - last_peak) < round(0.2*Fs)
            continue;
        end
        
        % Apply threshold
        if pk > THR_SIG
            % T‑wave discrimination (only if RR is very short)
            if ~isempty(detected)
                prev = detected(end);
                RR = (loc - prev) / Fs;
                if RR < 0.36
                    % Compare slopes in a 40 ms window around each peak
                    win_slope = round(0.04*Fs);
                    l1 = max(1, prev-win_slope);
                    r1 = min(length(ecg_bp), prev+win_slope);
                    l2 = max(1, loc-win_slope);
                    r2 = min(length(ecg_bp), loc+win_slope);
                    slope1 = max(abs(diff(ecg_bp(l1:r1))));
                    slope2 = max(abs(diff(ecg_bp(l2:r2))));
                    mean_RR = mean(RR_list + eps);
                    if slope2 < 0.5*slope1 || RR < 0.5*mean_RR
                        continue;   % likely a T‑wave
                    end
                end
            end
            % Accept QRS
            detected(end+1) = loc;
            last_peak = loc;
            SPKI = 0.125 * pk + 0.875 * SPKI;
            
            % Update RR list
            if length(detected) > 1
                RR = (detected(end) - detected(end-1)) / Fs;
                RR_list(end+1) = RR;
                if length(RR_list) > 8
                    RR_list(1) = [];
                end
                RR_avg = mean(RR_list);
            end
        else
            % Update noise estimate
            NPKI = 0.125 * pk + 0.875 * NPKI;
        end
        
        % Update thresholds
        THR_SIG = NPKI + 0.25 * (SPKI - NPKI);
        THR_NOISE = 0.5 * THR_SIG;
        
        % Searchback for missed beats (if interval > 1.66 * average RR)
        if ~isempty(detected) && RR_avg > 0
            if (loc - last_peak) > round(1.66 * RR_avg * Fs)
                search_start = last_peak;
                search_end = loc;
                if search_end > search_start
                    [missed_pk, idx] = max(mwi(search_start:search_end));
                    if missed_pk > THR_NOISE
                        true_loc = search_start + idx - 1;
                        detected(end+1) = true_loc;
                        last_peak = true_loc;
                        SPKI = 0.25 * missed_pk + 0.75 * SPKI;
                        % Recompute thresholds
                        THR_SIG = NPKI + 0.25 * (SPKI - NPKI);
                        THR_NOISE = 0.5 * THR_SIG;
                    end
                end
            end
        end
    end
    
    % Refine peak positions: find the highest point within 50 ms of each detection
    refined = [];
    for i = 1:length(detected)
        win = round(0.05 * Fs);
        l = max(1, detected(i) - win);
        r = min(length(ecg), detected(i) + win);
        [~, idx] = max(ecg(l:r));
        refined(end+1) = l + idx - 1;
    end
    peaks = unique(refined);
end

