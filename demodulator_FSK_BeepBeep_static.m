function demodulator_FSK_BeepBeep_static(MODE, pause_time_s, L_or_R)
global lastx
global x
global f0
global f1
global mode

global start_demodulator
global start_demodulator_stamp
global allcodes

global last_impulse_f0
global last_impulse_f1

global string_send_back

global sample_num_stamp
global send_TOF_tag_stamp

global send_id

global recording_data
global cur_recording_index

mode = MODE;
cur_recording_index = 1;

string_send_back = '';

allcodes = [];
last_impulse_f0 = [];
last_impulse_f1 = [];

start_demodulator = 0;
start_demodulator_stamp = 0;

sampleRate = 48000;
windows_size = 1024;
f0 = 4000;
f1 = 6000;
sample_num_stamp = 0;
send_id = 0;


recObj = audiorecorder(sampleRate, 24, 2);

disp('record start.');
record(recObj, 5);
pause(pause_time_s);
send_str(MODE, L_or_R);
pause(5.1 - pause_time_s);
recording_data = getaudiodata(recObj);
audiowrite('Cget.wav',recording_data,sampleRate);

% 
[recording_data, fs] = audioread('Aget.wav');

recording_data = recording_data(:,1); % 只取一个声道的麦克风声音。
hd = design(fdesign.bandpass('N,F3dB1,F3dB2',6,5500,12500,fs),'butter');
hd0 = design(fdesign.bandpass('N,F3dB1,F3dB2',6,5500,6500,fs),'butter');
hd1 = design(fdesign.bandpass('N,F3dB1,F3dB2',6,3500,4500,fs),'butter');
hold on;
plot(filter(hd1,recording_data));
plot(filter(hd0,recording_data));
recording_data = filter(hd0,recording_data) + filter(hd1,recording_data);
recording_data = filter(hd,recording_data);

% calc_relativity();

disp(fs);
disp('record finished. start analysis.');
% plot(recording_data);

lastx = zeros(1,4096);

while (1)
    x = my_record();
    
    if isempty(x)
        break;
    end
        
    lengthx = length(x);
    
    x = x';

    update_decode_fast(windows_size, sampleRate);
    
    sample_num_stamp = sample_num_stamp + length(x);

    lastx = x;
end
end

% function calc_relativity()
% global recording_data
% global relativity
% relativity = conv2(recording_data, modulator_FSK([0,1,0,1,0,1,0,1],-1));
% disp(max(relativity));
% end

function data=my_record()
global recording_data
global cur_recording_index

right = cur_recording_index+1024-1;
if right > length(recording_data)
    data = [];
else
    data = recording_data(cur_recording_index:right);
    cur_recording_index = cur_recording_index + 1024;
end
end

function impulse = Get_y_impulse(f, y, fs)
index_f = round(f/fs*size(y,2));
impulse = max(y(1:end,index_f-10:index_f+10)')';
end

function impulse = Get_single_y_impulse(f, y, fs)
index_f = round(f/fs*size(y,2));
impulse = max(y(1:end,index_f-2:index_f+2));
end


% function impulse = Get_impulse(f, samples, fs)
% index_f = round(f/fs*length(samples));
% y = abs(fft(samples));
% impulse = max(y(index_f-2:index_f+2, 1:end)')';
% end

function update_decode_fast(windows_size, sampleRate)
global lastx
global x
global f0
global f1
global start_demodulator
global start_demodulator_stamp
    
if start_demodulator
    demodulator_after_preamble(x, sampleRate, windows_size);
    return;
end


%% just for test
% lastx = modulator_FSK([0,1,0,1]);
% x = modulator_FSK([0,1,0,1]);
%% 根据录到的样本，得到所有窗口的样本们。存在矩阵里。

samples_matrix = zeros(length(lastx), windows_size);
idx = 0;

for i=length(lastx)-windows_size+2:length(lastx)
    idx = idx+1;
    samples_matrix(idx,:)=[lastx(i:end),x(1:windows_size-(length(lastx)-i+1))];
end

for i=1:1:length(x)-windows_size+1
    idx = idx+1;
    samples_matrix(idx,:)=x(i:i+windows_size-1);
end

%% 根据样本做fft，得到每个窗口的频谱图。
ys = abs(fft(samples_matrix')');

global last_impulse_f0
global last_impulse_f1

%% 获取各个窗口在频率 f0, f1 处的振幅
impulse_f0 = Get_y_impulse(f0, ys, sampleRate);
impulse_f1 = Get_y_impulse(f1, ys, sampleRate);

last_impulse_f0 = [last_impulse_f0; impulse_f0];
last_impulse_f1 = [last_impulse_f1; impulse_f1];

if length(last_impulse_f0) > 16 * windows_size
    last_impulse_f0 = last_impulse_f0(end-16*windows_size+1:end);
    last_impulse_f1 = last_impulse_f1(end-16*windows_size+1:end);
end

%% 设置根据 impulse 判断 0/1 的规则，无法识别的置为 -1

l = 1.6;
ri = 100;
while ri - l > 0.005
    mid = (l + ri) / 2;
    codes0 = last_impulse_f0 - mid * last_impulse_f1;
    codes1 = last_impulse_f1 - mid * last_impulse_f0;               
    codes = 0 * (codes0 > 0) + 1 * (codes1 > 0) + (-1) * ((codes0 < 0) .* (codes1 < 0));

    %% 将codes reshape 成一个矩阵，使得每一列为一个 id 所对应的所有 codes 
    codes = reshape(codes, windows_size, [])';          % 分成各个 id 的 f0, f1 振幅
    find_preamble_flag = find_preamble(codes, windows_size, sampleRate, 0);

    v = max(max(find_preamble_flag));

    if v == 8
        l = mid;
    else
        ri = mid;
    end
end

codes0 = last_impulse_f0 - l * last_impulse_f1;
codes1 = last_impulse_f1 - l * last_impulse_f0;               
codes = 0 * (codes0 > 0) + 1 * (codes1 > 0) + (-1) * ((codes0 < 0) .* (codes1 < 0));

%% 将codes reshape 成一个矩阵，使得每一列为一个 id 所对应的所有 codes 
codes = reshape(codes, windows_size, [])';          % 分成各个 id 的 f0, f1 振幅
find_preamble_flag = find_preamble(codes, windows_size, sampleRate, 1);

if start_demodulator == 1
    disp("the midsearch l = " + l);
last_impulse_f0 = [];
last_impulse_f1 = [];
end

end

%% 寻找前导码
function find_preamble_flag = find_preamble(allcodes, windows_size, sampleRate, auto_start_demodulator)
global x
global start_demodulator_stamp
global sample_num_stamp
global recording_data

codes = allcodes;

find_preamble_flag = zeros(size(allcodes));
for i=1:8
    correct_value = mod(i+1,2);                    % 前导码当前bit应该是多少
    find_preamble_flag = find_preamble_flag + (codes == correct_value);
    codes(1,:) = [];
    codes(end+1,:) = zeros(1, size(codes,2))-1;
end

%% 查找有哪些位置是符合前导码条件的
[pr, pc] = find(find_preamble_flag == 8);

% 如果存在位置能找到前导码
if (length(pr) > 0 && auto_start_demodulator)
    % disp("length of p = " + length(pr));
    p = 1;
    r = pr(p);
    c = pc(p);
    offset_from_begin = (r - 1) * windows_size + c + 7 * windows_size;
    offset_to_end = numel(allcodes) - offset_from_begin + 1;
    start_demodulator_stamp = sample_num_stamp + length(x) - offset_to_end + 1;
    start_preamble();
    demodulator_after_preamble(x(end-offset_to_end+1:end), sampleRate, windows_size);
end
end

%% 退出接收消息的状态
function clear_preamble()
global start_demodulator
global code_after_preamble
global preamble
preamble = [];
code_after_preamble = [];
start_demodulator = 0;
end

%% 进入接收消息的状态
function start_preamble()
global start_demodulator
start_demodulator = 1;
global code_after_preamble
global preamble
preamble = [];
code_after_preamble = [];
end

function res = array_to_int(data)
res = 0;
for i=1:length(data)
    res = res * 2 + data(i);
end
end

%% 如果当前处于接收消息的状态才会执行。用于将录制到的样本解调成数据包
function demodulator_after_preamble(datas, sampleRate, windows_size)
global preamble
global code_after_preamble
global f0
global f1
global start_demodulator
global length_of_LengthCode
global sample_num_stamp
global start_demodulator_stamp
global recording_data

length_of_LengthCode = 10;


if start_demodulator == 0
    return;
end

preamble = [preamble, datas];
endi = 1-windows_size;

%% 将录制到的样本解调成 01 码
for i=1:windows_size:length(preamble)-windows_size+1
    x = preamble(i:i+windows_size-1);
    b = start_demodulator_stamp + length(code_after_preamble)*windows_size;
    e = b + windows_size - 1;
    xr = recording_data(b:e)';
    delta = b - start_demodulator_stamp;
    
    if max(abs(xr - x)) ~= 0
        error("x and recording data not match.");
    end
    
    x = recording_data(b:e)';
    
    y = abs(fft(x')');
    impulse_f0 = Get_single_y_impulse(f0, y, sampleRate);
    impulse_f1 = Get_single_y_impulse(f1, y, sampleRate);
    % disp(impulse_f0 + " vs " + impulse_f1);
    

    plot(b:e,x/2,'g');
    
    if (impulse_f0 > impulse_f1)
        code_after_preamble = [code_after_preamble, 0];
    else
        code_after_preamble = [code_after_preamble, 1];
    end
    endi = i;
end
preamble = preamble(endi+windows_size:end);

%% 根据01码得到长度码、编码信息
if (length(code_after_preamble) > length_of_LengthCode + 8)
    length_to_recv = array_to_int(code_after_preamble(9:length_of_LengthCode + 8));
    if (length(code_after_preamble) > length_to_recv + length_of_LengthCode + 8)
        code_after_preamble = code_after_preamble(1:length_to_recv + length_of_LengthCode + 8);
        disp("recv finished, length = " + length(code_after_preamble));
        for i=-8:0
            plot([start_demodulator_stamp+i*windows_size,start_demodulator_stamp+i*windows_size],[0,0.02],'m','linewidth',2);
        end
        for i=1:length(code_after_preamble)
            b = start_demodulator_stamp + (i-1) * 1024;
            e = start_demodulator_stamp + (i-0) * 1024;
            if (code_after_preamble(i) == 1)
            plot([b,b], [0, 0.02], 'm', 'linewidth', 2);
            end
        end
        recv_code(code_after_preamble(19:end), sample_num_stamp + length(x) - start_demodulator_stamp);
        clear_preamble();
    end
end
end

%% 二进制串传字符串，用于处理接收到的01码
function [ str ] = bin2string( binary )
%UNTITLED2 此处显示有关此函数的摘要
%   把二进制串转化为字符串
binary = binary(1:end-mod(length(binary),8));
L = length(binary);
str = [];
binary = reshape(binary',[8,L/8]);
binary = binary';
for i=1:L/8
    s= 0;
    for j = 1:8
        s = s+2^(8-j)*binary(i,j);
    end
    str = [str,char(s)];
end
end

%% 接口函数。设置在收到码后做哪些事情
function recv_code(code, num_of_samples_during_recv)
global mode
global string_send_back
global start_demodulator_stamp
global send_TOF_tag_stamp
global sample_num_stamp
    str_recv = bin2string(code);
    disp("start_demodulator_stamp = " + start_demodulator_stamp);
    disp("str = "+ str_recv);
    if strcmp(mode, 'recv')
        disp("[RECV mode] num of samples_during_recv" + num_of_samples_during_recv);
        if strcmp(str_recv, 'ToF')
            string_send_back = num2str(num_of_samples_during_recv);
        end
    end
    
    if strcmp(mode, 'send') && ~isempty(str_recv)
        time_during_recv = str2num(str_recv);
        if ~isempty(time_during_recv)
            delta = start_demodulator_stamp - send_TOF_tag_stamp - time_during_recv;
            disp("start_demodulator_stamp = " + start_demodulator_stamp);
            disp("send_TOF_tag_stamp = " + send_TOF_tag_stamp);
            disp("time_during_recv = " + time_during_recv);
            disp("delta = " + delta);
        end
    end
end
