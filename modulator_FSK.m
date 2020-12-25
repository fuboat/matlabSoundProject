function data = modulator_FSK(codes, f0, f1, L_or_R)
sampleRate = 48000;
windows_size = 256;

t = 0:1/sampleRate:1;

data = zeros(1, windows_size * length(codes));

t0 = 0;
t1 = 0;

for i=1:length(codes)
    if codes(i) == 0
        f = f0;
        next_t_begin = t0;
    else
        f = f1;
        next_t_begin = t1;
    end
    % next_t_begin = 0;
    t = next_t_begin:1/sampleRate:next_t_begin+1/sampleRate*(windows_size+1);
    t = t(1:windows_size);
    data((i-1)*windows_size+1:i*windows_size) = sin(2*pi*f*t);
    if f == f0
        t0 = t0 + windows_size / sampleRate;
        t1 = t1 + windows_size / sampleRate * f0 / f1;
    else
        t1 = t1 + windows_size / sampleRate;
        t0 = t0 + windows_size / sampleRate *  f1 / f0;
    end
end

data = [zeros(1,1000),data,zeros(1,1000)];

if L_or_R == 1
 data = [data', zeros(length(data),1)];
end
if L_or_R == 2
 data = [zeros(length(data),1), data'];
end

player = audioplayer(data, sampleRate);
playblocking(player);
audiowrite('sound.wav',data,sampleRate);
end