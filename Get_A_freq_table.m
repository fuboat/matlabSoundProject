% [data, fs] 为音频数据（只支持单声道）。
% 输出的 [A, f] 表示在各个频率的幅值。如 freq = f(1) 处的幅值为 A(1).

function [A, f] = Get_A_freq_table(data, fs)
NFFT = size(data,1);
A = fft(data,NFFT);
f = ((0:NFFT-1)*fs/NFFT);

A = A(1:floor(NFFT/2))*2/NFFT;  % 幅值归一化
f = f(1:floor(NFFT/2));         % 因为FFT的特性，后半部分数据没有更多信息，舍弃。
end
