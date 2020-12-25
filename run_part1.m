function run_part1(csvfilename, wavfilename)
outputfilename = 'output.csv';
fs = csvread(csvfilename,0,1,[0,1,0,1]);
symbol_duration = csvread(csvfilename,1,1,[1,1,1,1]);
windows_size = fs * symbol_duration;
f0 = csvread(csvfilename,2,1,[2,1,2,1]);
f1 = csvread(csvfilename,3,1,[3,1,3,1]);
tic;
[positions,strs,codess] = demodulator_new_version(wavfilename,fs,windows_size,f0,f1,repmat([0,1],[1,10]),8,1);
toc;
write_to_csv(outputfilename, codess);

packet_total = 0;
packet_get = 0;
bit_total = 0;
bit_right = 0;

id = 1;

datas = csvread(csvfilename,5,0);
max_id = size(datas,1);

for id=1:1:max_id
    packet_total = packet_total + 1;
    len = csvread(csvfilename,4+id,2,[4+id,2,4+id,2]);
    real_codes = csvread(csvfilename,4+id,4,[4+id,4,4+id,4+len-1]);
    offset = csvread(csvfilename,4+id,3,[4+id,3,4+id,3]);
    p = find(positions >= offset-100 & positions <= offset + 100,1);
    % p = id;
    if ~isempty(p)
        % disp(offset);
        % disp(positions(p));
        packet_get = packet_get + 1;
        codes = codess{p};
        bit_total = bit_total + len;
        de_len = codes(20+1:20+8);
        de_len = array_to_int(de_len);
        check_len = min([de_len,len]);
        bit_right = bit_right + sum(codes(20+8+1:20+8+check_len) == real_codes(1:check_len));
    end
    id = id + 1;
end
disp("总bit数为 " + bit_total + "，正确bit数为 " + bit_right + "，误码率=" + (bit_total-bit_right) / bit_right);
disp("总数据包个数为 " + max_id + "，检测到的数据个数为 " + packet_get + "，丢包率=" + (max_id-packet_get) / max_id);
disp("解码结果已保存至 " + outputfilename);
end

function res = array_to_int(data)
res = 0;
for i=1:length(data)
    res = res * 2 + data(i);
end
end

function write_to_csv(csvfilename, codess)
fid = fopen(csvfilename,'w');
for i=1:1:length(codess)
   codes = codess{i};
   fprintf(fid,'%d',length(codes)-28);
   for j=29:1:length(codes)
       fprintf(fid,',%d',codes(j));
   end
   fprintf(fid,'\n');
end
fclose(fid);
end