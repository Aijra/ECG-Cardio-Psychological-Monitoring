
clc; clear; close all;

%% ========= DEMOGRAPHICS =========
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

%% ========= AGE GROUP FUNCTION =========
get_age_group = @(age) ...
    (isnan(age))*"Unknown" + ...
    (age<30 & ~isnan(age))*"Young" + ...
    (age>=30 & age<50)*"Adult" + ...
    (age>=50 & age<70)*"Senior" + ...
    (age>=70)*"Elderly";

%% ========= SETTINGS =========
records = [100:101 103:124 200:234];
Fs = 360;

results = [];
results_struct = [];
k = 0;


total_TP = 0;
total_FP = 0;
total_FN = 0;
total_TN = 0;

% Psychological stats
psych_labels = {};
cardiac_labels = {};

%% ========= MAIN LOOP =========
for i = 1:length(records)

    rec = records(i);
    fprintf('\nProcessing Record %d...\n', rec);

    %% LOAD ECG
    filename = ['/MATLAB Drive/mit-bih-arrhythmia-database-1.0.0/' num2str(rec) '.dat'];
    fid = fopen(filename,'r');

    if fid == -1
        warning('Missing ECG file'); continue;
    end

    raw = fread(fid,[3 inf],'uint8'); fclose(fid);

    ecg = bitshift(raw(2,:),4) + bitshift(raw(1,:),-4);
    ecg(ecg > 2047) = ecg(ecg > 2047) - 4096;

    %% LOAD CSV
    csv_file = [num2str(rec) '_annotations_1.csv'];
    if ~isfile(csv_file)
        warning('Missing CSV'); continue;
    end

    ann = readtable(csv_file);
    valid_symbols = {'N','L','R','A','V'};
    mask = ismember(ann.annotation_symbol, valid_symbols);
    true_peaks = double(ann.index(mask));
    true_peaks(true_peaks > length(ecg)) = [];

   
    %% DETECT QRS
    detected_peaks = qrs_detect(ecg, Fs);

    %% EVALUATION
    tol = round(0.1*Fs);
    TP=0; FP=0; matched=zeros(length(true_peaks),1);

    for j=1:length(detected_peaks)
        [d, idx] = min(abs(true_peaks - detected_peaks(j)));

        if d<=tol && matched(idx)==0
            TP=TP+1; matched(idx)=1;
        else
            FP=FP+1;
        end
    end

    FN = sum(matched==0);

    % TN approximation (non-peak samples)
    %TN = length(ecg) - (TP + FP + FN); -- gives fake accuracy
    %Ignoring TN completely -we're gonna use F1

    %F1 = 2*TP / (2*TP + FP + FN);
    
    total_TP = total_TP + TP;
    total_FP = total_FP + FP;
    total_FN = total_FN + FN;
    %total_TN = total_TN + TN;

    Se = TP/(TP+FN+eps);
    P  = TP/(TP+FP+eps);
    Acc = TP/(TP+FP+FN+eps);

    fprintf('Rec %d → Se: %.3f | P: %.3f | Acc: %.3f\n', rec, Se, P, Acc);

    results = [results; rec Se P Acc];
    %% 
% %% %% ========= STEP 7: PLOT IF ACCURACY >= 95% =========
% if Acc >= 0.95
%     % Define zoom window (10 seconds around first detected peak)
%     zoom_duration = 10; % seconds
%     if ~isempty(detected_peaks)
%         start_idx = max(1, detected_peaks(1) - round(Fs*1)); % 1 sec before first peak
%         end_idx   = min(length(ecg), start_idx + round(Fs*zoom_duration) - 1);
%     else
%         start_idx = 1;
%         end_idx = min(length(ecg), round(Fs*zoom_duration));
%     end
% 
%     time = (start_idx:end_idx)/Fs;  % time vector
%     segment = ecg(start_idx:end_idx); 
% 
%     figure('Name',['Record ' num2str(rec)], 'NumberTitle','off');
%     plot(time, segment, 'b'); hold on;
% 
%     % Plot detected peaks within segment
%     seg_detected = detected_peaks(detected_peaks >= start_idx & detected_peaks <= end_idx);
%     plot(seg_detected/Fs, ecg(seg_detected), 'ro', 'MarkerFaceColor','r');
% 
%     % Plot true peaks within segment
%     seg_true = true_peaks(true_peaks >= start_idx & true_peaks <= end_idx);
%     plot(seg_true/Fs, ecg(seg_true), 'g*');
% 
%     xlabel('Time (s)'); ylabel('Amplitude');
%     title(['ECG Record ' num2str(rec) ' → Acc: ' num2str(Acc*100,'%.2f') '%']);
%     legend('ECG','Detected Peaks','True Peaks');
%     grid on;
% end
    %% HRV CALCULATION (INSIDE LOOP)
    if length(detected_peaks) > 5

        RR = diff(detected_peaks)/Fs;
        RR = RR(RR>0.4 & RR<1.2);

        if ~isempty(RR)

            RR = movmedian(RR,3);

            HR = 60./RR;
            mean_HR = mean(HR);
            SDNN = std(RR);
            RMSSD = sqrt(mean(diff(RR).^2));

            info = demo(rec);
            age = info.age;

            %% Adaptive thresholds
            if age > 70
                RMSSD_low=0.010; RMSSD_high=0.035; SDNN_low=0.020;
            elseif age > 50
                RMSSD_low=0.015; RMSSD_high=0.040; SDNN_low=0.025;
            else
                RMSSD_low=0.020; RMSSD_high=0.050; SDNN_low=0.030;
            end

            %% Psychological state
            % Better HRV-based classification
            if mean_HR > 100 && RMSSD < RMSSD_low
                psych = 'ACUTE STRESS';
            elseif RMSSD < RMSSD_low
                psych = 'MILD STRESS';
            elseif RMSSD >= RMSSD_low && RMSSD <= RMSSD_high
                psych = 'NORMAL';
            else
                % Add stricter condition for relaxed
                if mean_HR < 75
                    psych = 'RELAXED';
                else
                    psych = 'NORMAL';
                end
            end

            %% Cardiac state
            if mean_HR>100
                cardiac='TACHYCARDIA';
            elseif mean_HR<55
                cardiac='BRADYCARDIA';
            elseif SDNN<SDNN_low
                cardiac='LOW HRV';
            elseif SDNN>0.1
                cardiac='ARRHYTHMIA RISK';
            else
                cardiac='NORMAL';
            end

            psych_labels{end+1} = psych;
            cardiac_labels{end+1} = cardiac;

            %% STORE RESULTS
            k = k + 1;
            results_struct(k).record = rec;
            results_struct(k).HR = mean_HR;
            results_struct(k).RMSSD = RMSSD;
            results_struct(k).SDNN = SDNN;
            results_struct(k).psych = psych;
            results_struct(k).cardiac = cardiac;
        end
    end
end

%% ========= FINAL OUTPUT =========
fprintf('\nFINAL RESULTS:\n');
disp(results);

fprintf('Mean Sensitivity: %.3f\n', mean(results(:,2)));
fprintf('Mean Precision: %.3f\n', mean(results(:,3)));
fprintf('Mean Accuracy: %.3f\n', mean(results(:,4)));

results_table = struct2table(results_struct);
disp(results_table);

%% ========= CONFUSION MATRIX =========
fprintf('\n===== OVERALL CONFUSION MATRIX =====\n');
fprintf('TP: %d\n', total_TP);
fprintf('FP: %d\n', total_FP);
fprintf('FN: %d\n', total_FN);
fprintf('TN: %d\n', total_TN);

overall_accuracy = (total_TP + total_TN) / ...
                   (total_TP + total_FP + total_FN + total_TN);

fprintf('Overall Accuracy: %.4f\n', overall_accuracy);

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
%% ========= IMPROVED QRS FUNCTION =========
function detected_peaks = qrs_detect(ecg, Fs)

    % --- Step 1: Preprocessing ---
    ecg = ecg - mean(ecg);
    ecg = detrend(ecg);

    [b,a] = butter(3,[5 20]/(Fs/2),'bandpass');
    ecg_f = filtfilt(b,a,ecg);

    % --- Step 2: Envelope detection ---
    env = abs(hilbert(ecg_f));

    % Smooth envelope
    win = round(0.12*Fs);
    env_smooth = movmean(env, win);

    % Normalize
    env_smooth = env_smooth / max(env_smooth + eps);

    % --- Step 3: Adaptive threshold ---
    thresh = mean(env_smooth) + 0.15*std(env_smooth);

    % --- Step 4: Peak detection ---
    [~, locs] = findpeaks(env_smooth, ...
        'MinPeakHeight', thresh, ...
        'MinPeakDistance', round(0.25*Fs));

    % --- Step 5: Map to actual R-peaks ---
    detected_peaks = [];

    for i = 1:length(locs)

        left = max(1, locs(i) - round(0.1*Fs));
        right = min(length(ecg), locs(i) + round(0.1*Fs));

        [~, idx] = max(ecg(left:right));
        detected_peaks(end+1) = left + idx - 1;
    end

    detected_peaks = unique(detected_peaks);

    % --- Step 6: Remove false peaks using RR ---
    if length(detected_peaks) > 1
        RR = diff(detected_peaks)/Fs;

        % keep physiologically valid intervals
        valid = [true RR > 0.35 & RR < 1.5];
        detected_peaks = detected_peaks(valid);
    end
end

fprintf('\nF1 Score: %.4f\n', ...
    2*(total_TP)/(2*total_TP + total_FP + total_FN));


F1_score = 2*total_TP / (2*total_TP + total_FP + total_FN);

fprintf('\n===== FINAL MODEL PERFORMANCE =====\n');
fprintf('F1 Score: %.4f\n', F1_score);
fprintf('Precision: %.4f\n', total_TP/(total_TP+total_FP));
fprintf('Sensitivity: %.4f\n', total_TP/(total_TP+total_FN));


