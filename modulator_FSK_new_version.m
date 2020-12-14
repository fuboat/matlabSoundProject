function data = modulator_FSK_new_version(codes, sampleRate, windows_size, f0, f1)
t = 0:1/sampleRate:1;
t = t(1:windows_size);
data = zeros(1, windows_size * length(codes));

t0 = 0;
t1 = 0;

for i=1:length(codes)
    if codes(i) == 0
        f = f0;
    else
        f = f1;
    end
    data((i-1)*windows_size+1:i*windows_size) = sin(2*pi*f*t);
end
end
