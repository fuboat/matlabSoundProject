function data = modulator_FSK(codes)
sampleRate = 48000;
windows_size = 256;
f0 = 20000;
f1 = 21000;

t = 0:1/sampleRate:1;
t = t(1:windows_size);

data = zeros(1, windows_size * length(codes));
for i=1:length(codes)
    if codes(i) == 0
        data((i-1)*windows_size+1:i*windows_size) = sin(2*pi*f0*t);
    elseif codes(i) == 1
        data((i-1)*windows_size+1:i*windows_size) = sin(2*pi*f1*t);
    end
end

sound(data, sampleRate);

end