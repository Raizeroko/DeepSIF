function process_raw_nmm(filename, varargin)
% Scan through the raw NMM data to find the spike data

% %%%%%%%%%%%%%% SETUP PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p是一个参数解析器对象，由于matlab中不像c++那样在函数传参要指明参数类型，
% 所以就用一个参数解析器判断对应名称的参数是不是想要的类型
p = inputParser;
%判断名称叫'filename'的参数，即函数传的第一个参数是否是个字符串类型
addRequired(p,'filename',@ischar);
% 在参数解析器中添加一个名为 'leadfield_name' 的参数，设置其默认值为 'leadfield_75_20k.mat'
% leadfield_75_20k.mat:leadfield matrix for fsaverage5, 75 channels
addParameter(p,'leadfield_name','leadfield_75_20k.mat', @ischar);
parse(p, filename, varargin{:})
filename = p.Results.filename;
headmodel = load(['../anatomy/' p.Results.leadfield_name]);
fwd = headmodel.fwd;
savefile_path = '../source/';

% ------------------------------------------------------------------------%
%------------------python里只生成了0-2-------------------------------------%
%iter_list = 1:5;   
iter_list = 0:2;   % the iter during NMM generation.
previous_iter_spike_num = zeros(1, 994);
for i_iter = 1:length(iter_list)
    iter = iter_list(i_iter);
    if isempty(dir([savefile_path 'nmm_' filename '/clip_info/iter' int2str(iter)]))
       mkdir([savefile_path 'nmm_' filename '/clip_info/iter' int2str(iter)])
    end

    % ------- Resume running if the process was interupted ----------------
    done = dir([savefile_path 'nmm_' filename '/clip_info/iter' int2str(iter) '/iter_' int2str(iter) '_i_*']);
    finished_regions = zeros(1, length(done));
    %-----------这个循环什么作用-------------------------------------------%
    for i = 1:length(done)
        finished_regions(i) = str2num(done(i).name(10:end-3));
    end
    %994个区域还有哪些区域没处理
    remaining_regions = setdiff(1:994, finished_regions+1);
    if isempty(remaining_regions)
        continue;
    end

    % -------- start the main progress -----------------------------------%
    %----------不能理解remaining_regions的作用，无法修改--------------------%
    %----------不能理解183:183,明显有错-------------------------------------%
    %----------待确认ii的范围----------------------------------------------%
    % for ii = 183:183%length(remaining_regions)    
    for ii = 1:10%length(remaining_regions)

        i = remaining_regions(ii);
        % creat folders to save nmm files
        if isempty(dir([savefile_path 'nmm_' filename '/a' int2str(i-1)]))
            mkdir([savefile_path 'nmm_' filename '/a' int2str(i-1)])
        end
        
        %fn = [savefile_path 'raw_nmm_200s/a' int2str(i-1) '/mean_iter_' int2str(iter) '_a_iter_' int2str(i-1)];
        fn = [savefile_path 'raw_nmm/a' int2str(i-1) '/mean_iter_' int2str(iter) '_a_iter_' int2str(i-1)];
        raw_data = load([fn '_ds.mat']);
        nmm = raw_data.data;
        % nmm = downsample(nmm, 4);
        [spike_time, spike_chan] = find_spike_time(nmm);                   % Process raw tvb output to find the spike peak time
        
        % ----------- select the spikes we want to extract ---------------%
        %----------这里的i用的不合理，待确认i的用途--------------------------%
        rule1 =  (spike_chan == i);                                        % there is spike in the source region
        start_time = floor(spike_time(rule1)/500) * 500 + 1;               % there is no source in other region in the clip
        rep = repmat(start_time', [900, 1]);
        clear_ind = rep + (-200:699)';                                     % 900 * num_spike
        rule2 = (sum(ismember(clear_ind, spike_time(~rule1)), 1) == 0);    % there are no other spikes in the clip
        spike_time = spike_time(rule1);
        %---------------------待确认用途-----------------------------------%
        %spike_time = spike_time(rule2);
               
        % ----------- Optional :  Scale the NMM here----------------------%
        alpha_value = find_alpha(nmm, fwd, i, spike_time, 15);
        nmm = rescale_nmm_channel(nmm, i, spike_time, alpha_value);       
        % ------------Save Spike NMM Data --------------------------------%
        start_time = floor(spike_time/500) * 500 + 1;
        rep = repmat(start_time', [500, 1]);
        spike_ind = rep + (0:499)';       
%         start_time = floor((spike_time+200)/500) * 500 + 1 - 200;        % start time can be changed
%         start_time = max(start_time, 101);
%         spike_ind = repmat(start_time, [500, 1]) + (0:499)';

        nmm_data = reshape(nmm(spike_ind,:), 500, [], size(nmm,2));        % size: time * num_spike * channel
        save_spikes_(nmm_data, [savefile_path 'nmm_' filename '/a' int2str(i-1) '/nmm_'], previous_iter_spike_num(i));
        previous_iter_spike_num(i) = previous_iter_spike_num(i) + length(spike_time);
        % Save something in clip info, so that we can make sure we finish this process
        save_struct = struct();
        save_struct.num_spike = previous_iter_spike_num(i);
        save_struct.spike_time = spike_time;
        parsave([savefile_path 'nmm_' filename '/clip_info/iter' int2str(iter) '/iter_' int2str(iter) '_i_' int2str(i-1) '.mat'], save_struct)
        sprintf(['iter_' int2str(iter) '_i_%d is done\n'], i-1)
    end % END REGION
end % END ITER
end % ENG FUNCTION



%% --------------------- Helper functions ------------------------------ %%
function parsave(fname, mapObj)
    save(fname,'-struct', 'mapObj');
end


function save_spikes_(spike_data, savefile_path, previous_iter_spike_num)
% Save the spike data into seperate files
% INPUTS: spike_data: time * num_spikes * channel; extracted spike data
%         savefile_path: string
    for iii = 1:size(spike_data,2)
        % The raw data
        data = squeeze(spike_data(:,iii,:));
        save([savefile_path int2str(iii+previous_iter_spike_num) '.mat'], 'data', '-v7')
    end
end


function [spike_time, spike_chan] = find_spike_time(nmm)
% Process raw tvb output to find the spike peak time.
%
% INPUTS:
%     - nmm        : (Downsampled) raw tvb output, time * channel
% OUTPUTS:
%     - spike_time : the spike peak time in the (downsampled) NMM data
%     - spike_chan : the spike channel for each spike

    spikes_nmm = nmm;
    spikes_nmm(nmm < 8) = 0;                                               % find the spiking activity stronger than the background
    local_max = islocalmax(spikes_nmm);                                    % find the peak
    [spike_time, spike_chan] = find(local_max);

    % 先排除峰值
    use_ind = (spike_time-249 > 0) & (spike_time+250 < size(nmm, 1)) & ([1 diff(spike_time)'>100]');
    spike_time = spike_time(use_ind);
    spike_chan = spike_chan(use_ind);

    % 然后再进行排序
    [spike_time, sort_ind] = sort(spike_time);
    spike_chan = spike_chan(sort_ind);


    %先排序再排除峰值明显有错%
    %[spike_time, sort_ind] = sort(spike_time);
    %spike_chan = spike_chan(sort_ind);                                     % sort the activity based on time
    %spike_time-249>0:排除前250s的峰值
    %spike_time+250<size(nmm,1):排除后250s的峰值
    %[1 diff(spike_time)'>100']:排除靠的太近的峰值
    % ignore the spikes at the beginning or end of the signal 
    % ignore peaks close together for now (will have signals with close peaks in multi-source condition)
    % 找到错误：全部变成0了
    %use_ind = (spike_time-249 > 0) & (spike_time+250 < size(nmm, 1)) & ([1 diff(spike_time)'>100]');                                   
    %spike_time = spike_time(use_ind)';
    %spike_chan = spike_chan(use_ind)';
    
end



function [alpha] = find_alpha(nmm, fwd, region_id, time_spike, target_SNR)
% Find the scaling factor for the NMM channels.
%
% INPUTS:
%     - nmm        : NMM data with single activation, time * channel
%     - fwd        : leadfield matrix, num_electrode * num_regions
%     - region_id  : source regions, start at 1, region_id(1) is the center
%     - time_spike : spike peak time
%     - target_SNR : set snr between signal and the background activity.
% OUTPUTS:
%     - spike_time : the spike peak time in the (downsampled) NMM data
%     - spike_chan : the spike channel for each spike
%     - alpha      : the scaling factor for one patch source
    rep = repmat(time_spike', [200, 1]);
    spike_ind = rep + (-99:100)';
    spike_ind = min(max(spike_ind(:),0), size(nmm,1));                     % make sure the index is not out of range
%     spike_ind = max(0, time_spike-100): max(time_spike+100,size(nmm,1));   % make sure the index is not out of range                     
    spike_shape = nmm(:,region_id(1)); %/max(nmm(:,region_id(1)));
    nmm(:, region_id) = repmat(spike_shape,1,length(region_id));
    % calculate the scaling factor
    [Ps, Pn, ~] = calcualate_SNR(nmm, fwd, region_id, spike_ind);
    alpha = sqrt(10^(target_SNR/10)*Pn/Ps);
end


function scaled_nmm = rescale_nmm_channel(nmm, region_id, spike_time, alpha_value)
% Re-scaling NMM channels in source channels
%
% INPUTS:
%     - nmm        : NMM data with single activation, time * channel
%     - spike_time : spike peak time
%     - region_id  : source regions, start at 1
%     - alpha_value: scaling factor
% OUTPUTS:
%     - scaled_nmm : scaled NMM in the source region; time * channel

    nmm_rm = nmm - mean(nmm, 1);
    for i=1:length(spike_time)
        sig = nmm_rm(spike_time(i)-249:spike_time(i)+250, region_id);      % one second data around the peak

        thre = 0.1;
        small_ind = find(abs(sig)<thre);                                   % the index that the values are close to 0
        small_ind((small_ind>450) | (small_ind < 50)) = [];
        start_ind = find((small_ind-250)<0);                               % spike start time
        % test 1
        while isempty(start_ind)
            thre = thre+0.05;
            small_ind = find(abs(sig)<thre);
            small_ind((small_ind>450) | (small_ind < 50)) = [];
            start_ind = find((small_ind-250)<0);
        end
        start_sig = small_ind(start_ind(end));

        % test 2
        [~, min_ind] = min(sig(301:400));
        min_ind = min_ind + 301;
        end_ind = find((small_ind-min_ind)>0);
        while isempty(end_ind)
            thre = thre+0.05;
            small_ind = find(abs(sig)<thre);
            small_ind((small_ind>450) | (small_ind < 50)) = [];
            end_ind = find((small_ind-min_ind)>0);
        end
        end_sig = small_ind(end_ind(1));                                   % spike end time

        sig(start_sig:end_sig) = sig(start_sig:end_sig) * alpha_value;     % scale the signal
        nmm_rm(spike_time(i)-249:spike_time(i)+250, region_id) = sig;
    end
    scaled_nmm = nmm_rm + mean(nmm, 1);


end


function [Ps, Pn, cur_snr] = calcualate_SNR(nmm, fwd, region_id, spike_ind)
% Caculate SNR at sensor space.
%
% INPUTS:
%     - nmm        : NMM data with single activation, time * channel
%     - fwd        : leadfield matrix, num_electrode * num_regions
%     - region_id  : source regions, start at 1, region_id(1) is the center
%     - spike_ind  : index to calculate the spike snr
% OUTPUTS:
%     - Ps         : signal power
%     - Pn         : noise power
%     - cur_snr    : current SNR in dB

    sig_eeg = (fwd(:, region_id)*nmm(:, region_id)')';   % time * channel
    sig_eeg_rm = sig_eeg - mean(sig_eeg, 1);
    dd = 1:size(nmm,2);
    dd(region_id) = [];
    noise_eeg = (fwd(:,dd)*nmm(:,dd)')';
    noise_eeg_rm = noise_eeg - mean(noise_eeg, 1);

    Ps = norm(sig_eeg_rm(spike_ind,:),'fro')^2/length(spike_ind);
    Pn = norm(noise_eeg_rm(spike_ind,:),'fro')^2/length(spike_ind);
    cur_snr = 10*log10(Ps/Pn);
end