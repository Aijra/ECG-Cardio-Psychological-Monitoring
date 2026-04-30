%% 
echo off

%% ================== LOAD ECG ==================
filename = '/MATLAB Drive/mit-bih-arrhythmia-database-1.0.0/101.dat';

fid = fopen(filename, 'r');
data = fread(fid, [2, inf], 'bit12'); %MIT-BIH format
fclose(fid);

ecg = double(data(1, :)); % channel 1
Fs = 360;
t = (0:length(ecg)-1)/Fs;

sec = 10;
N = sec * Fs;

figure;
plot(t(1:N), ecg(1:N));
title('Original ECG Signal(Record101) - First 10 Seconds');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% ================== PREPROCESSING ==================
% Remove DC
ecg = ecg - mean(ecg);

% Bandpass filter (BEST RANGE)
[b, a] = butter(3, [0.5 40]/(Fs/2), 'bandpass');
ecg_filtered = filtfilt(b, a, ecg);

% Optional smoothing (improves detection)
ecg_filtered = medfilt1(ecg_filtered, 5);

sec = 10;
N = sec * Fs;

figure;
plot(t(1:N), ecg_filtered(1:N));
title('Preprocessed ECG Signal(Record101) - First 10 Seconds');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% ================== LOAD ANNOTATIONS ==================
ann = readtable('101_annotations_1.csv');
true_peaks = double(ann.index);
true_peaks(true_peaks > length(ecg_filtered)) = [];

%% ================== PAN-TOMPKINS ==================

% Derivative
derivative = diff(ecg_filtered);
derivative(end+1) = derivative(end);

% Squaring
squared = derivative .^ 2;

% Moving Window Integration
window_size = round(0.12 * Fs); % 120 ms
mwi = conv(squared, ones(1,window_size)/window_size, 'same');

% Normalize
mwi = mwi / max(mwi);

%% ================== PEAK DETECTION ==================

threshold = mean(mwi) + 0.6 * std(mwi);
min_distance = round(0.25 * Fs); % 250 ms

[~, locs_mwi] = findpeaks(mwi, ...
    'MinPeakHeight', threshold, ...
    'MinPeakDistance', min_distance);

%% ================== PEAK LOCALIZATION ==================

search_window = round(0.08 * Fs);
locs = zeros(size(locs_mwi));

for i = 1:length(locs_mwi)
    start_idx = max(locs_mwi(i) - search_window, 1);
    end_idx   = min(locs_mwi(i) + search_window, length(ecg_filtered));
    
    segment = ecg_filtered(start_idx:end_idx);
    
    [~, idx] = max(abs(segment)); % more robust
    locs(i) = start_idx + idx - 1;
end

detected_peaks = unique(locs);

%% ================== PEAK REFINEMENT ==================

RR = diff(detected_peaks) / Fs;
RR = RR(:); % force column

% Remove unrealistic beats
valid = [true; (RR > 0.35 & RR < 1.5)];
detected_peaks = detected_peaks(valid);

disp(['Detected Peaks: ', num2str(length(detected_peaks))]);

%% ================== ACCURACY ==================

tolerance = round(0.1 * Fs); % 100 ms (strict)

TP = 0;
FP = 0;
matched = zeros(length(true_peaks),1);

for i = 1:length(detected_peaks)
    
    diff_array = abs(true_peaks - detected_peaks(i));
    [min_diff, idx] = min(diff_array);
    
    if min_diff <= tolerance && matched(idx) == 0
        TP = TP + 1;
        matched(idx) = 1;
    else
        FP = FP + 1;
    end
end

FN = sum(matched == 0);

Accuracy = TP / (TP + FP + FN);
Sensitivity = TP / (TP + FN);
Precision = TP / (TP + FP);


fprintf('\n===== PERFORMANCE =====\n');
fprintf("Accuracy: %.2f%%\n", Accuracy*100);
fprintf("Sensitivity: %.2f%%\n", Sensitivity*100);
fprintf("Precision: %.2f%%\n", Precision*100);

%% ================== DEBUG PLOT ==================

sec = 10;
N = sec * Fs;

figure;
plot(t(1:N), ecg_filtered(1:N)); hold on;

% filter detected peaks within 15 sec
det_peaks_10 = detected_peaks(detected_peaks <= N);
plot(t(det_peaks_10), ecg_filtered(det_peaks_10), 'ro');

% filter true peaks within 15 sec
true_peaks_10 = true_peaks(true_peaks <= N);
plot(t(true_peaks_10), ecg_filtered(true_peaks_10), 'g*');

legend('ECG','Detected R-peaks','Annotated R-peaks');
title('R-Peak Alignment(Record101) (First 10 Seconds)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% ================== HRV ==================

RR = diff(detected_peaks) / Fs;

% Filter RR
RR = RR(RR > 0.4 & RR < 1.2);

median_RR = median(RR);
RR = RR(abs(RR - median_RR) < 0.15);

RR = movmedian(RR, 3);

% Compute instantaneous HR
HR = 60 ./ RR;

% --- Define 10-second windows ---
window_size = 10;  % seconds
t_min = floor(min(t_rr));
t_max = ceil(max(t_rr));

edges = t_min:window_size:t_max;

% Preallocate
HR_10s = [];
t_10s = [];

% --- Compute mean HR in each window ---
for i = 1:length(edges)-1
    idx = t_rr >= edges(i) & t_rr < edges(i+1);
    
    if any(idx)
        HR_10s(end+1) = mean(HR(idx)); %#ok<SAGROW>
        t_10s(end+1) = edges(i) + window_size/2; %#ok<SAGROW>
    end
end

% --- Classification on 10s averaged HR ---
tachy_idx = HR_10s > 100;
brady_idx = HR_10s < 60;
normal_idx = HR_10s >= 60 & HR_10s <= 100;

% --- Plot ---
figure;
plot(t_10s, HR_10s, 'm', 'LineWidth', 2); hold on;

scatter(t_10s(tachy_idx), HR_10s(tachy_idx), 40, 'r', 'filled');
scatter(t_10s(brady_idx), HR_10s(brady_idx), 40, 'b', 'filled');
scatter(t_10s(normal_idx), HR_10s(normal_idx), 30, 'g', 'filled');

yline(100, '--r', 'Tachycardia');
yline(60, '--b', 'Bradycardia');

xlabel('Time (s)');
ylabel('Heart Rate (BPM)');
title('Heart Rate (10-Second Averaged)');
legend('HR (10s avg)','Tachycardia','Bradycardia','Normal');

grid on;

mean_HR = mean(HR);
SDNN = std(RR);
RMSSD = sqrt(mean(diff(RR).^2));

%% ================== RESULTS ==================

disp(['Mean HR: ', num2str(mean_HR)]);
disp(['SDNN: ', num2str(SDNN)]);
disp(['RMSSD: ', num2str(RMSSD)]);
