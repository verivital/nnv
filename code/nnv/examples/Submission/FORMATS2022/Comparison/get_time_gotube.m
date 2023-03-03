function time = get_time_gotube(fname)
str = fileread(fname); % dedicated for reading files as text 
data = jsondecode(str);
time = data.notes.total_time;
end