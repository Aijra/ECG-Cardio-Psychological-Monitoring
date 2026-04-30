%% 
clc; clear; close all;

%% ================= DEMOGRAPHICS =================
demo = containers.Map('KeyType','double','ValueType','any');
demo(100)=struct('age',69); demo(101)=struct('age',75); demo(102)=struct('age',84);
demo(103)=struct('age',NaN); demo(104)=struct('age',66); demo(105)=struct('age',73);
demo(106)=struct('age',24); demo(107)=struct('age',63); demo(108)=struct('age',76);
demo(109)=struct('age',77); demo(111)=struct('age',50); demo(112)=struct('age',39);
demo(113)=struct('age',74); demo(114)=struct('age',69); demo(115)=struct('age',72);
demo(116)=struct('age',66); demo(117)=struct('age',66); demo(118)=struct('age',66);
demo(119)=struct('age',63); demo(121)=struct('age',51); demo(122)=struct('age',63);
demo(123)=struct('age',59); demo(124)=struct('age',65);

demo(200)=struct('age',64); demo(201)=struct('age',64); demo(202)=struct('age',64);
demo(203)=struct('age',74); demo(205)=struct('age',66); demo(207)=struct('age',61);
demo(208)=struct('age',69); demo(209)=struct('age',65); demo(210)=struct('age',64);
demo(212)=struct('age',67); demo(213)=struct('age',66); demo(214)=struct('age',53);
demo(215)=struct('age',55); demo(217)=struct('age',66); demo(219)=struct('age',53);
demo(220)=struct('age',64); demo(221)=struct('age',70); demo(222)=struct('age',67);
demo(223)=struct('age',59); demo(228)=struct('age',22); demo(230)=struct('age',71);
demo(231)=struct('age',66); demo(232)=struct('age',76); demo(233)=struct('age',70);
demo(234)=struct('age',55);

%% ================= SETTINGS =================
records = [100:124 200:234];
Fs = 360;

%% ================= STORAGE =================
Se_PT=[]; P_PT=[]; F1_PT=[]; Err_PT=[];
Se_WT=[]; P_WT=[]; F1_WT=[]; Err_WT=[];
Se_HT=[]; P_HT=[]; F1_HT=[]; Err_HT=[];

Se_HY=[]; P_HY=[]; F1_HY=[]; Err_HY=[];

%% ================= MAIN LOOP =================
for i = 1:length(records)

    rec = records(i);
    fprintf('Processing Record %d\n', rec);

    %% -------- LOAD ECG --------
    file = ['/MATLAB Drive/mit-bih-arrhythmia-database-1.0.0/' num2str(rec) '.dat'];
    fid = fopen(file,'r');
    if fid==-1, continue; end

    raw = fread(fid,[3 inf],'uint8'); fclose(fid);

    ecg = bitshift(raw(2,:),4) + bitshift(raw(1,:),-4);
    ecg(ecg>2047)=ecg(ecg>2047)-4096;

    %% -------- LOAD ANNOTATIONS --------
    csv = [num2str(rec) '_annotations_1.csv'];
    if ~isfile(csv), continue; end

    ann = readtable(csv);
    mask = ismember(ann.annotation_symbol,{'N','L','R','A','V'});
    true_peaks = double(ann.index(mask));
    true_peaks(true_peaks>length(ecg))=[];

    %% ================= METHOD 1: PAN-TOMPKINS =================
    det1 = pan_tompkins(ecg,Fs);
    [Se,P,F1,Err] = evaluate(det1,true_peaks,Fs);
    Se_PT=[Se_PT Se]; P_PT=[P_PT P]; F1_PT=[F1_PT F1]; Err_PT=[Err_PT Err];

    %% ================= METHOD 2: WAVELET =================
    det2 = wavelet_qrs(ecg,Fs);
    [Se,P,F1,Err] = evaluate(det2,true_peaks,Fs);
    Se_WT=[Se_WT Se]; P_WT=[P_WT P]; F1_WT=[F1_WT F1]; Err_WT=[Err_WT Err];

    %% ================= METHOD 3: HILBERT =================
    det3 = hilbert_qrs(ecg,Fs);
    [Se,P,F1,Err] = evaluate(det3,true_peaks,Fs);
    Se_HT=[Se_HT Se]; P_HT=[P_HT P]; F1_HT=[F1_HT F1]; Err_HT=[Err_HT Err];

end

%% ================= FINAL RESULTS =================
fprintf('\n=========== FINAL COMPARISON ===========\n');

fprintf('\nPan-Tompkins:\nSe=%.3f | P=%.3f | F1=%.3f | Err=%.3f\n', ...
    mean(Se_PT),mean(P_PT),mean(F1_PT),mean(Err_PT));

fprintf('\nWavelet:\nSe=%.3f | P=%.3f | F1=%.3f | Err=%.3f\n', ...
    mean(Se_WT),mean(P_WT),mean(F1_WT),mean(Err_WT));

fprintf('\nHilbert:\nSe=%.3f | P=%.3f | F1=%.3f | Err=%.3f\n', ...
    mean(Se_HT),mean(P_HT),mean(F1_HT),mean(Err_HT));

%% ================= F1 COMPARISON =================
figure;
bar([mean(F1_PT) mean(F1_WT) mean(F1_HT)]);
set(gca,'XTickLabel',{'Pan-Tompkins','Wavelet','Hilbert'});
ylabel('F1 Score'); title('QRS Detection Comparison');
grid on;

%% ================= METRICS =================
figure;

subplot(1,3,1);
bar([mean(Se_PT) mean(Se_WT) mean(Se_HT)]);
title('Sensitivity'); grid on;

subplot(1,3,2);
bar([mean(P_PT) mean(P_WT) mean(P_HT)]);
title('Precision'); grid on;

subplot(1,3,3);
bar([mean(F1_PT) mean(F1_WT) mean(F1_HT)]);
set(gca,'XTickLabel',{'Pan-Tompkins','Wavelet','Hilbert'});
title('F1 Score'); grid on;

%% ================= ERROR PLOT =================
figure;
bar([mean(Err_PT) mean(Err_WT) mean(Err_HT)]);
set(gca,'XTickLabel',{'Pan-Tompkins','Wavelet','Hilbert'});
ylabel('Error Rate');
title('Error Rate Comparison');
grid on;

%% ================= EVALUATION =================
function [Se,P,F1,Err] = evaluate(det, true_peaks, Fs)

tol = round(0.1*Fs);
TP=0; FP=0;
matched=zeros(length(true_peaks),1);

for i=1:length(det)
    [d,idx]=min(abs(true_peaks-det(i)));

    if d<=tol && matched(idx)==0
        TP=TP+1; matched(idx)=1;
    else
        FP=FP+1;
    end
end

FN=sum(matched==0);

Se = TP/(TP+FN+eps);
P  = TP/(TP+FP+eps);
F1 = 2*TP/(2*TP+FP+FN+eps);

% ✅ ERROR RATE
Err = (FP + FN)/(TP + FP + FN + eps);

end

%% ================= PAN-TOMPKINS =================
function r_locs = pan_tompkins(ecg, Fs)

[b,a] = butter(2, [5 15]/(Fs/2), 'bandpass');
ecg_f = filtfilt(b,a,ecg);

d = diff(ecg_f); d(end+1)=d(end);
sq = d.^2;
mwi = movmean(sq, round(0.15*Fs));

thr = mean(mwi) + 0.5*std(mwi);

[~,locs] = findpeaks(mwi,'MinPeakHeight',thr,'MinPeakDistance',round(0.25*Fs));

r_locs = refine_peaks(ecg, locs, Fs);

end

%% ================= WAVELET =================
function r_locs = wavelet_qrs(ecg, Fs)

[c,l] = wavedec(ecg,5,'db4');
qrs = wrcoef('d',c,l,'db4',3) + wrcoef('d',c,l,'db4',4);

energy = movmean(qrs.^2, round(0.12*Fs));
thr = mean(energy) + 0.6*std(energy);

[~,locs] = findpeaks(energy,'MinPeakHeight',thr,'MinPeakDistance',round(0.25*Fs));

r_locs = refine_peaks(ecg, locs, Fs);

end

%% ================= HILBERT =================
function r_locs = hilbert_qrs(ecg, Fs)

[b,a] = butter(2,[5 20]/(Fs/2),'bandpass');
ecg_f = filtfilt(b,a,ecg);

env = abs(hilbert(ecg_f));
env = movmean(env, round(0.1*Fs));

thr = mean(env) + 0.7*std(env);

[~,locs] = findpeaks(env,'MinPeakHeight',thr,'MinPeakDistance',round(0.25*Fs));

r_locs = refine_peaks(ecg, locs, Fs);

end

%% ================= PEAK REFINEMENT =================
function r_locs = refine_peaks(ecg, locs, Fs)

r_locs = zeros(size(locs));
win = round(0.08*Fs);

for i=1:length(locs)
    left=max(locs(i)-win,1);
    right=min(locs(i)+win,length(ecg));
    [~,idx]=max(ecg(left:right));
    r_locs(i)=left+idx-1;
end

end

