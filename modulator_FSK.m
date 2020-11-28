function data = modulator_FSK(codes)
sampleRate = 48000;
windows_size = 512;
f0 = 20000;
f1 = 18000;

t = 0:1/sampleRate:1;

data = zeros(1, windows_size * length(codes));

lastA = 0;

for i=1:length(codes)
    if codes(i) == 0
        f = f0;
    else
        f = f1;
    end
    t = asin(lastA)/(2*pi*f):1/sampleRate:asin(lastA)/(2*pi*f)+1/sampleRate*(windows_size+1);
    t = t(1:windows_size);
    data((i-1)*windows_size+1:i*windows_size) = sin(2*pi*f*t);
    lastA = sin(2*pi*f*t(windows_size));
end

sound(data, sampleRate);
plot(data);

end