function [positions_of_premble] = demodulator_new_version(filename, fs, windows_size, f0, f1, premble_array, length_of_length_code, i_channel)
    [data, ~] = audioread(filename);
    
    data = data';
    
    if ~isempty(i_channel) && i_channel == 2
        data = data(2,:);
    else
        data = data(1,:);
    end
    
    plot(data);
    hold on;
    
%     hd0 = design(fdesign.bandpass('N,F3dB1,F3dB2',6,min(f0,f1)-100,max(f0,f1)+100,fs),'butter');
%     hd1 = design(fdesign.bandpass('N,F3dB1,F3dB2',6,f1-100,f1+100,fs),'butter');
%     data = filter(hd0,data);
    
%     hd0 = design(fdesign.bandpass('N,F3dB1,F3dB2',6,f0-100,f0+100,fs),'butter');
%     hd1 = design(fdesign.bandpass('N,F3dB1,F3dB2',6,f1-100,f1+100,fs),'butter');
%     data = filter(hd0,data)+filter(hd1,data);
%     plot(data);
    positions_of_premble=demodulator_data(data, fs, windows_size, f0, f1, premble_array, length_of_length_code);
end

function positions_of_premble=demodulator_data(data, fs, windows_size, f0, f1, premble_array, length_of_length_code)
    data_conved_0 = conv(data, flip(modulator_FSK_new_version([0], fs, windows_size, f0, f1)), 'valid');
    data_conved_1 = conv(data, flip(modulator_FSK_new_version([1], fs, windows_size, f0, f1)), 'valid');
    data_conved_self_pre = data .* data;
    
    for i=2:length(data)
        data_conved_self_pre(i) = data_conved_self_pre(i) + data_conved_self_pre(i-1);
    end

    data_conved_self = data_conved_self_pre(windows_size * length(premble_array) + 1:end) - data_conved_self_pre(1:end - windows_size * length(premble_array));
    data_conved_premble = zeros(1, length(data));
    
    premble_sample = modulator_FSK_new_version(premble_array, fs, windows_size, f0, f1);
    premble_sample_conved_sum = sum(premble_sample .* premble_sample);
    
    for i=1:length(data) - windows_size * length(premble_array)
        for j=1:length(premble_array)
            index = i+(j-1)*windows_size;
            if index > length(data_conved_0)
                break;
            end
            if premble_array(j) == 0
                delta = data_conved_0(index);
            else
                delta = data_conved_1(index);
            end
            data_conved_premble(i) = data_conved_premble(i) + delta;
        end
        data_conved_premble(i) = data_conved_premble(i) / sqrt(data_conved_self(i)) / sqrt(premble_sample_conved_sum);
    end
    
    [relative_values,positions] = sort(data_conved_premble);
    relative_values = flip(relative_values);
    positions = flip(positions);
    
    disable_premble = zeros(1,length(data));
        
    positions_of_premble = [];
    
    for p_index=1:length(positions)
        p = positions(p_index);
        if ~disable_premble(p)
            %%%%
            
            codes = demodulator_after_preamble(data(p:end), fs, windows_size, f0, f1, length(premble_array), length_of_length_code);
            disable_premble(max(1,p-windows_size*length(premble_array)+1):p+length(codes)*windows_size-1) = 1;
            positions_of_premble = [positions_of_premble, p];
            % disp(p);
            plot([p,p],[0,1],'m','linewidth',2);
            
            %%%% 
           

%             codes = demodulator_after_preamble(data(p:end), fs, windows_size, f0, f1, length(premble_array), length_of_length_code);
%             str = bin2string(codes(length(premble_array) + length_of_length_code+1:end));
%             if strcmp(str, 'send')
%                 disp(str);
%                 disable_premble(max(1,p-windows_size*length(premble_array)+1):p+length(codes)*windows_size-1) = 1;
%                 positions_of_premble = [positions_of_premble, p];
%                 
%             disp(p);
%             end
%             plot([p,p],[0,0.2],'m','linewidth',2);
%             pause(0.01);

            %%%%
        end
        
        if (relative_values(p_index) < 0.3)
            break;
        end
    end
    
    disp(length(positions_of_premble));
    
    [positions_of_premble, ~] = sort(positions_of_premble);
    disp(positions_of_premble');
end


function str=demodulator_after_premble_to_str(data, fs, windows_size, f0, f1, length_of_premble, length_of_length_code)
codes = demodulator_after_preamble(datas, sampleRate, windows_size, f0, f1, length_of_premble, length_of_length_code);
str=bin2string(codes(length_of_premble + length_of_length_code+1:end));
end


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

function codes=demodulator_after_preamble(datas, sampleRate, windows_size, f0, f1, length_of_premble, length_of_length_code)
endi = 1-windows_size;
codes = [];

%% 将录制到的样本解调成 01 码
for i=1:windows_size:length(datas)-windows_size+1
    x = datas(i:i+windows_size-1);
    y = abs(fft(x')');
    impulse_f0 = Get_single_y_impulse(f0, y, sampleRate);
    impulse_f1 = Get_single_y_impulse(f1, y, sampleRate);
    if (impulse_f0 > impulse_f1)
        codes = [codes, 0];
    else
        codes = [codes, 1];
    end
    
    if (length(codes) >= length_of_premble + length_of_length_code)
        length_to_recv = array_to_int(codes(length_of_premble+1:length_of_premble + length_of_length_code));
        if (length(codes) >= length_of_premble + length_of_length_code + length_to_recv)
            break;
        end
    end
end
end

% function p=find_next_premble(data, fs, windows_size, f0, f1, premble_array, length_of_length_code, max_find_length)
% p = find_next_premble_by_relative(data, fs, windows_size, f0, f1, premble_array, length_of_length_code, max_find_length);
% end


% %% 根据01码得到长度码、编码信息
% if (length(code_after_preamble) > length_of_LengthCode + 8)
%     length_to_recv = array_to_int(code_after_preamble(9:length_of_LengthCode + 8));
%     if (length(code_after_preamble) > length_to_recv + length_of_LengthCode + 8)
%         code_after_preamble = code_after_preamble(1:length_to_recv + length_of_LengthCode + 8);
%         disp("recv finished, length = " + length(code_after_preamble));
%         for i=-8:0
%             plot([start_demodulator_stamp+i*windows_size,start_demodulator_stamp+i*windows_size],[0,0.02],'m','linewidth',2);
%         end
%         for i=1:length(code_after_preamble)
%             b = start_demodulator_stamp + (i-1) * 1024;
%             e = start_demodulator_stamp + (i-0) * 1024;
%             if (code_after_preamble(i) == 1)
%             plot([b,b], [0, 0.02], 'm', 'linewidth', 2);
%             end
%         end
%         recv_code(code_after_preamble(19:end), sample_num_stamp + length(x) - start_demodulator_stamp);
%         clear_preamble();
%     end
% end


function res = array_to_int(data)
res = 0;
for i=1:length(data)
    res = res * 2 + data(i);
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


% function p=find_next_premble_by_FFT(data, fs, windows_size, f0, f1, premble_array, ~, max_find_length)
% samples_matrix = [];
% idx = 0;
% 
% for i=1:1:min(length(data)-windows_size+1,max_find_length)
%     idx = idx+1;
%     samples_matrix(idx,:)=data(i:i+windows_size-1);
% end
% 
% ys = abs(fft(samples_matrix')');
% 
% impulse_f0 = Get_y_impulse(f0, ys, fs);
% impulse_f1 = Get_y_impulse(f1, ys, fs);
% 
% l = 1.6;
% ri = 100;
% while ri - l > 0.005
%     mid = (l + ri) / 2;
%     codes0 = impulse_f0 - mid * impulse_f1;
%     codes1 = impulse_f1 - mid * impulse_f0;               
%     codes = 0 * (codes0 > 0) + 1 * (codes1 > 0) + (-1) * ((codes0 < 0) .* (codes1 < 0));
% 
%     codes = reshape(codes, windows_size, [])';
%     [find_preamble_flag, ~] = find_preamble(codes, windows_size, fs, premble_array);
% 
%     v = max(max(find_preamble_flag));
%     if v > 1
%         disp(v);
%     end
%     if v == length(premble_array)
%         l = mid;
%     else
%         ri = mid;
%     end
% end
% 
% codes0 = impulse_f0 - l * impulse_f1;
% codes1 = impulse_f1 - l * impulse_f0;               
% codes = 0 * (codes0 > 0) + 1 * (codes1 > 0) + (-1) * ((codes0 < 0) .* (codes1 < 0));
% 
% codes = reshape(codes, windows_size, [])';
% [~, p] = find_preamble(codes, windows_size, fs, premble_array);
% 
% end
% 
% 
% function [find_preamble_flag, p] = find_preamble(allcodes, windows_size, ~, premble_array, ~)
% codes = allcodes;
% 
% find_preamble_flag = zeros(size(allcodes));
% for i=1:length(premble_array)
%     find_preamble_flag = find_preamble_flag + (codes == premble_array(i));
%     codes(1,:) = [];
%     codes(end+1,:) = zeros(1, size(codes,2))-1;
% end
% 
% [pr, pc] = find(find_preamble_flag == length(premble_array));
% 
% if (~isempty(pr))
%     r = pr(1);
%     c = pc(1);
%     p = (r - 1) * windows_size + c;
% else
%     p = -1;
% end
% 
% end
% 

