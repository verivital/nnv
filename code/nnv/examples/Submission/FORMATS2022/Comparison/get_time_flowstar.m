function  time = get_time_flowstar(fname)
ttt = fileread(fname);
idx = strfind(ttt,'Total time');
if idx
    time = str2num(ttt(idx+21:idx+27));
else
    time = 'Timeout';
end

