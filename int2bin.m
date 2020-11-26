function [binary] = int2bin(x, len)
binary = [];
while x
    binary = [mod(x,2),binary];
    x = (x - mod(x,2)) / 2;
end
binary = [zeros(1,len - length(binary)), binary];
end
