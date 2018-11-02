function out=mycellstring2double(in)
for k=1:length(in)
    out{k}=str2num(in{k});
end