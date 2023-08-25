SDKPATH = 'C:\Users\jakek\Downloads\TDTMatlabSDK\TDTSDK';
addpath(genpath(SDKPATH));

str_root_dir = 'C:\';
str_user = 'Users';
str_user_name = 'jakek';
str_downloads = 'Downloads';
str_context = 'behavior_scs_experiments';
str_exp = '01_13';
str_channel = 'channel1';
str_plot= 'sem';
d_input = fullfile(str_root_dir, str_user, str_user_name, str_downloads, str_context, str_exp, str_channel);
muscle_type = {'bicep', 'fcr', 'ecr'};

data = TDTbin2mat('G:\My Drive\cl_data\01_12\2023-01-12_grasphold002_sub-cl10\2023-01-12_sub-cl10_s02_r01_aie');
T = load('G:\My Drive\cl_data\01_12\2023-01-12_grasphold002_sub-cl10\2023-01-12_sub-cl10_s02_r01_aie\2023-01-12_sub-cl10_s02_r01_aie.mat');


str_target_stream = 'EMGw';
% channel 1 - biceps; channel 2 - fcr; channel 3 - ecr
fs_target = data.streams.(str_target_stream).fs;
t_target = 0:1/fs_target:(length(data.streams.(str_target_stream).data)-1)/fs_target;
y_target = data.streams.(str_target_stream).data;
system_latency = data.epocs.TSyn.onset;

% m = 1e-3;
% y_target(y_target>+m)= +m;
% y_target(y_target<-m)= -m;

temp_table = T.data.trial;

% sanity check for total number of led trials to compare with actual # of
% registered trials
bb_trial = 0;
for ix_bb_trial = 1:height(temp_table)
    if strcmp(temp_table.stim_type(ix_bb_trial), 'led000') == 1
        bb_trial = bb_trial + 1;
    elseif strcmp(temp_table.stim_type(ix_bb_trial), 'led010') == 1
        bb_trial = bb_trial + 1;
    end
end

% initializing table to store Synapse-stored onset times that match actual
% succesful trials stored in post-execution exported file
temp_table.lever_on_marker = zeros(height(temp_table), 1);
for ix_good_reach = 1:length(data.epocs.TLeN.onset)
    if temp_table.time_trigger_lever_hold(ix_good_reach) > 0 
       if data.epocs.TLeN.onset(ix_good_reach) - temp_table.time_trigger_lever_on(ix_good_reach) < 2
          temp_table.lever_on_marker(ix_good_reach) = data.epocs.TLeN.onset(ix_good_reach); 
       end
    else
        continue
    end
end

% now iterating through recorded beam break times and case matching with
% lever hold times
temp_table.beam = zeros(height(temp_table), 1);
for ix_beam = 1:length(data.epocs.bbe_.onset)
    for ix_good_reach = 1:length(data.epocs.TLeN.onset)
        if abs(data.epocs.bbe_.onset(ix_beam) - data.epocs.TLeN.onset(ix_good_reach)) < 1
           temp_table.beam(ix_good_reach) = data.epocs.bbe_.onset(ix_beam);
        else
            continue
        end
    end
end


figure(1);
naming_counter = 1;
y_filt_bs = NaN(3, length(y_target));
for ix_channel = 1:3
    % Yao 2021 movement onset detection for EMG

    % 5th order bandpass filter (1); passband ripple = 0.2hz due to low
    % frequency artifact observed in the signal


    % Notch filter at 60 hz to remove line noise
    % Boxcar filter for smoothing
    
    dee = designfilt('highpassiir', 'FilterOrder', 5, 'PassbandFrequency', 200, 'PassbandRipple', 0.2, ...
    'SampleRate', fs_target);

    y_filt = filtfilt(dee, double(y_target(ix_channel, :)));

    dee2 = designfilt('highpassiir', 'StopbandFrequency', 59, 'PassbandFrequency', 61, ...
    'DesignMethod', 'butter', 'SampleRate', fs_target);

    y_filt_bs(ix_channel, :) = filtfilt(dee2, double(y_filt));
    
    %   Sanity check
%     plot(t_target, y_filt_bs*1000000);
end

% Initializing storage array for calculating across bb trials
bb_count = 0;
for ix_trials = 1:height(temp_table.lever_on_marker)
    if temp_table.lever_on_marker(ix_trials) > 0 && temp_table.beam(ix_trials) > 0
       t_start = temp_table.beam(ix_trials);
       t_lever = temp_table.lever_on_marker(ix_trials);
       [~, ix_start] = min(abs(t_target - t_start));
       [~, ix_lever] = min(abs(t_target - t_lever));
       t_noise = t_start - 2;
       [~, ix_plot] = min(abs(t_target - t_noise));
       vec_plot = ix_plot:(ix_start + fs_target * 2);

       epoch_length = length(vec_plot);
       bb_count = bb_count + 1;
    else
        continue
    end
end

%Create array to store epochs for each bb event
mean_bicep_response = NaN(epoch_length, bb_count);
mean_fcr_response = NaN(epoch_length, bb_count);
mean_ecr_response = NaN(epoch_length, bb_count);

sem_bicep_response = NaN(bb_count, epoch_length);
sem_fcr_response = NaN(bb_count, epoch_length);
sem_ecr_response = NaN(bb_count, epoch_length);


%% visualizing trials on case by case basis and storing each event 
figure(1);clf;
counter = 1;
for ix_trials = 1:height(temp_table.lever_on_marker)
    if temp_table.lever_on_marker(ix_trials) > 0 && temp_table.beam(ix_trials) > 0
       t_start = temp_table.beam(ix_trials);
       t_lever = temp_table.lever_on_marker(ix_trials);
       [~, ix_start] = min(abs(t_target - t_start));
       [~, ix_lever] = min(abs(t_target - t_lever));
       t_noise = t_start - 2;
       [~, ix_plot] = min(abs(t_target - t_noise));
       vec_plot = ix_plot:(ix_start + fs_target * 2);

       subplot(3, 1, 1);
       sem_bicep_response(counter, :) = y_filt_bs(1, vec_plot);
       smoothed_bicep = smooth(abs(y_filt_bs(1, vec_plot)), fs_target * 0.020);
       mean_bicep_response(:, counter) = smoothed_bicep;
       plot(t_target(1, vec_plot), smoothed_bicep);
       xline(t_target(:, ix_start), 'k--', {'BB onset'});
       xline(t_target(:, ix_lever), 'k--', {'Lever onset'});
       xlabel('Time (s)');
       ylabel('Amplitude (\muV)');
       str_title = append('EMG response for', ' ', muscle_type(1, 1));
       title(str_title);

       subplot(3, 1, 2);
       sem_fcr_response(counter, :) = y_filt_bs(2, vec_plot);
       smoothed_fcr = smooth(abs(y_filt_bs(2, vec_plot)), fs_target * 0.020);
       mean_fcr_response(:, counter) = smoothed_fcr;
       plot(t_target(1, vec_plot), smoothed_fcr);
       xline(t_target(:, ix_start), 'k--', {'BB onset'});
       xline(t_target(:, ix_lever), 'k--', {'Lever onset'});
       xlabel('Time (s)');
       ylabel('Amplitude (\muV)');
       str_title = append('EMG response for', ' ', muscle_type(1, 2));
       title(str_title);

       subplot(3, 1, 3)
       sem_ecr_response(counter, :) = y_filt_bs(3, vec_plot);
       smoothed_ecr = smooth(abs(y_filt_bs(3, vec_plot)), fs_target * 0.020);
       mean_ecr_response(:, counter) = smoothed_ecr;
       plot(t_target(1, vec_plot), abs(smoothed_ecr));
       xline(t_target(:, ix_start), 'k--', {'BB onset'});
       xline(t_target(:, ix_lever), 'k--', {'Lever onset'});   
       xlabel('Time (s)');
       ylabel('Amplitude (\muV)');
       str_title = append('EMG response for', ' ', muscle_type(1, 3));
       title(str_title);
       counter = counter + 1;

       str_title_store = sprintf('Trial #%0.1d.png', ix_trials);
       fullFileName = fullfile(d_input, str_title_store);
       print(fullFileName, '-dpng');

    else
        continue
    end
end


%Plotting all trials and mean
figure(2);clf;hold on
tc_emg = [-2, 2];
for ix = 1:width(mean_bicep_response)
    idx_window = round(tc_emg(1) * fs_target):round(tc_emg(2) * fs_target);
    t_plot = idx_window / fs_target;
    
    subplot(3, 1, 1);hold on;
    plot(t_plot, mean_bicep_response(:, ix));
    subplot(3, 1, 2);hold on;
    plot(t_plot, mean_fcr_response(:, ix));
    subplot(3, 1, 3);hold on;
    plot(t_plot, mean_ecr_response(:, ix));

end

%must calculate std/sem before averaging matrix (mean behavior)
y_bicep_std = std(mean_bicep_response, 0, 2);
y_bicep_sem = y_bicep_std./sqrt(width(mean_bicep_response));

y_fcr_std = std(mean_fcr_response, 0, 2);
y_fcr_sem = y_fcr_std./sqrt(width(mean_bicep_response));

y_ecr_std = std(mean_ecr_response, 0, 2);
y_ecr_sem = y_ecr_std./sqrt(width(mean_bicep_response));

mean_bicep_response = mean(mean_bicep_response, 2);
mean_fcr_response = mean(mean_fcr_response, 2);
mean_ecr_response = mean(mean_ecr_response, 2);

subplot(3, 1, 1);hold on
% mean_bicep_response = smooth(abs(mean_bicep_response), fs_target * 0.020);
plot(t_plot, mean_bicep_response, 'k', 'Linewidth', 1.5);
title('Averaged bicep activity relative to beam break');
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
xline(0, 'k--', {'BB onset'});


subplot(3, 1, 2);hold on
% mean_fcr_response = smooth(abs(mean_fcr_response), fs_target * 0.020);
plot(t_plot, mean_fcr_response, 'k', 'Linewidth', 1.5);
title('Averaged FCR activity relative to beam break');
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
xline(0, 'k--', {'BB onset'});

subplot(3, 1, 3);hold on
% mean_ecr_response = smooth(abs(mean_ecr_response), fs_target * 0.020);
plot(t_plot, mean_ecr_response, 'k', 'Linewidth', 1.5);
title('Averaged ECR activity relative to beam break');
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
xline(0, 'k--', {'BB onset'});


str_title_store = sprintf('Trial #%0.1d.png', ix_trials);
fullFileName = fullfile(d_input, str_title_store);
print(fullFileName, '-dpng');

%% plotting sem w/ciplot
figure(4);hold on;
plot(t_plot, mean_bicep_response, 'r');
plot(t_plot, mean_bicep_response - y_bicep_sem, 'k--');
plot(t_plot, mean_bicep_response + y_bicep_sem, 'k--');
xline(0, 'k--', {'BB onset'});
title('SEM of averaged bicep activity during beam break');
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
ciplot((mean_bicep_response - y_bicep_sem), (mean_bicep_response + y_bicep_sem), t_plot, [0, 0, 1], 0.5);
muscle = muscle_type(1);
fullFileName = fullfile(str_root_dir, str_user, str_user_name, str_downloads, str_context, str_exp, str_plot, muscle);
print(fullFileName, '-dpng');


clf;hold on
plot(t_plot, mean_fcr_response, 'r');
plot(t_plot, mean_fcr_response - y_fcr_sem, 'k--');
plot(t_plot, mean_fcr_response + y_fcr_sem, 'k--');
xline(0, 'k--', {'BB onset'});
title('SEM of averaged FCR activity during beam break');
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
ciplot((mean_fcr_response - y_fcr_sem), (mean_fcr_response + y_fcr_sem), t_plot, [0, 0, 1], 0.5);
muscle = muscle_type(2);
fullFileName = fullfile(str_root_dir, str_user, str_user_name, str_downloads, str_context, str_exp, str_plot, muscle);
print(fullFileName, '-dpng');



clf;hold on
plot(t_plot, mean_ecr_response, 'r');
plot(t_plot, mean_ecr_response - y_ecr_sem, 'k--');
plot(t_plot, mean_ecr_response + y_ecr_sem, 'k--');
xline(0, 'k--', {'BB onset'});
title('SEM of averaged ECR activity during beam break');
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
ciplot((mean_ecr_response - y_ecr_sem), (mean_ecr_response + y_ecr_sem), t_plot, [0, 0, 1], 0.5);
muscle = muscle_type(3);
fullFileName = fullfile(str_root_dir, str_user, str_user_name, str_downloads, str_context, str_exp, str_plot, muscle);
print(fullFileName, '-dpng');


%% EMG thresholding to determine start of movement

% 1-D median filtering

for ix = 1:3
    y_filt_bs(ix, :) = medfilt1(abs(y_filt_bs(ix, :)), round(fs_target * 0.05));
    y_filt_bs(ix, :) = (y_filt_bs(ix, :) - median(y_filt_bs(ix, :)))./mad(y_filt_bs(ix, :));
end


% Medium absolute deviation - use to discriminate physiological signal from
% noise artifact spikes in EMG

figure(3);clf;hold on;
for ix_trials = 1:height(temp_table.lever_on_marker)
    if temp_table.lever_on_marker(ix_trials) > 0 && temp_table.beam(ix_trials) > 0
       t_start = temp_table.beam(ix_trials);
       t_lever = temp_table.lever_on_marker(ix_trials);
       [~, ix_start] = min(abs(t_target - t_start));
       [~, ix_lever] = min(abs(t_target - t_lever));
       t_noise = t_start - 1;
       [~, ix_plot] = min(abs(t_target - t_noise));
       vec_plot = ix_plot:(ix_start + fs_target * 2);
    
       subplot(3, 1, 1);
       plot(t_target(1, vec_plot), y_filt_bs(1, vec_plot)./mad(y_filt_bs(1, vec_plot)));
       plot(t_target(1, vec_plot), y_filt_bs(1, vec_plot));
       xline(t_target(:, ix_start), 'k--', {'BB onset'});
       xline(t_target(:, ix_lever), 'k--', {'Lever onset'});
        
       subplot(3, 1, 2);
       plot(t_target(1, vec_plot), y_filt_bs(2, vec_plot) ./mad(y_filt_bs(2, vec_plot)));
       plot(t_target(1, vec_plot), y_filt_bs(2, vec_plot));
       xline(t_target(:, ix_start), 'k--', {'BB onset'});
       xline(t_target(:, ix_lever), 'k--', {'Lever onset'});
        
       subplot(3, 1, 3);
       plot(t_target(1, vec_plot), y_filt_bs(3, vec_plot) ./mad(y_filt_bs(3, vec_plot)));
       plot(t_target(1, vec_plot), y_filt_bs(3, vec_plot));
       xline(t_target(:, ix_start), 'k--', {'BB onset'});
       xline(t_target(:, ix_lever), 'k--', {'Lever onset'});
 
       str_title_store = sprintf('Trial #%0.1d.png', ix_trials);
       fullFileName = fullfile(d_input, str_title_store);
       print(fullFileName, '-dpng');
   
    else
        continue
    end
end

% Plot all trials on a single plot

figure(4);clf;hold on
tc_emg = [-2, 2];
for ix_trials = 1:height(temp_table.lever_on_marker)
    if temp_table.beam(ix_trials) > 0
       idx_window = round(tc_emg(1) * fs_target):round(tc_emg(2) * fs_target);
       t_start = temp_table.beam(ix_trials);
       [~, ix_start] = min(abs(t_target - t_start));
       ix_plot = ix_start + idx_window;
       t_plot = idx_window / fs_target;
        
       plot(t_plot, y_filt_bs(1, ix_plot)./mad(y_filt_bs(1, ix_plot)));

    else
       continue
    end
end

plot(t_plot, mean_bicep_response(1, :)./mad(mean_bicep_response(1, :)), 'k', 'LineWidth', 1.5);

xline(0, 'k--', {'BB onset'});
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
title('Bicep EMG across recorded beam break trials');

%%
clf;hold on;
tc_emg = [-2, 2];
for ix_trials = 1:height(temp_table.lever_on_marker)
    if temp_table.beam(ix_trials) > 0
       idx_window = round(tc_emg(1) * fs_target):round(tc_emg(2) * fs_target);
       t_start = temp_table.beam(ix_trials);
       [~, ix_start] = min(abs(t_target - t_start));
       ix_plot = ix_start + idx_window;
       t_plot = idx_window / fs_target;
        
       plot(t_plot, y_filt_bs(2, ix_plot)./mad(y_filt_bs(2, ix_plot)));

       % Calculation of mean EMG behavior for beam break marked trials
       

    else
       continue
    end
end

plot(t_plot, mean_fcr_response(1, :)./mad(mean_fcr_response(1, :)), 'k', 'LineWidth', 1.5);

xline(0, 'k--', {'BB onset'});
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
title('FCR EMG across recorded beam break trials');

%%
clf;hold on;
tc_emg = [-2, 2];
for ix_trials = 1:height(temp_table.lever_on_marker)
    if temp_table.beam(ix_trials) > 0
       idx_window = round(tc_emg(1) * fs_target):round(tc_emg(2) * fs_target);
       t_start = temp_table.beam(ix_trials);
       [~, ix_start] = min(abs(t_target - t_start));
       ix_plot = ix_start + idx_window;
       t_plot = idx_window / fs_target;
        
       plot(t_plot, y_filt_bs(3, ix_plot)./mad(y_filt_bs(3, ix_plot)));

       % Calculation of mean EMG behavior for beam break marked trials
       

    else
       continue
    end
end

plot(t_plot, mean_ecr_response(1, :)./mad(mean_ecr_response(1, :)), 'k', 'LineWidth', 1.5);

xline(0, 'k--', {'BB onset'});
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
title('ECR EMG across recorded beam break trials');

    %% Note: ciplot must be before graph parameters or else graph specifications get wiped out
ciplot(y_mea - y_sem, y_mea + y_sem, t_slice * 1000, [0, 0, 0], 0.5);
