function pf=myreadjson(file)
fid = fopen(file); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
pf = jsondecode(str);
end