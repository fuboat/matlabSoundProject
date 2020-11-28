function [allcode] = send_str(str)
code = string2bin(str)';
lengthCode = int2bin(length(code), 10);
preamble = [0,1,0,1,0,1,0,1,1,1,1,1,1,1,1,1];
allcode = [zeros(1,100), preamble, lengthCode, code, zeros(1,100)];
modulator_FSK(allcode);
end

function [binary] = int2bin(x, len)
binary = [];
while x
    binary = [mod(x,2),binary];
    x = (x - mod(x,2)) / 2;
end
binary = [zeros(1,len - length(binary)), binary];
end

function [ binary ] = string2bin( str )
%   把字符串转换成二进制串
ascii = abs(str);
L = length(ascii);
binary = zeros(L,8);
for i=1:L
    binary_str = dec2bin(ascii(i));
    binary_str_index = length(binary_str);
    for j = 8:-1:1
        if binary_str_index >0
            binary(i,j) = str2num(binary_str(binary_str_index));
        else
            binary(i,j) = 0;
        end
        binary_str_index = binary_str_index-1;
    end
end
binary = reshape(binary',[L*8,1]);
end