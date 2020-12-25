function [positions_of_premble, strs] = recv_str_new_version(filename, fs, windows_size, f0, f1, premble_array, length_of_length_code, i_channel)
recObj = audiorecorder(fs, 24, 2);
recordblocking(recObj,15);
recording_data = getaudiodata(recObj);
audiowrite(filename, recording_data, fs);
[positions_of_premble, strs] = demodulator_new_version_for_recv(filename, fs, windows_size, f0, f1, premble_array, length_of_length_code, i_channel);
end
