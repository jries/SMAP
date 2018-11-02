function output=field2cell(input)
ind=strfind(input,'.');
ind=[0 ind length(input)+1];
for k=1:length(ind)-1
    output{k}=input(ind(k)+1:ind(k+1)-1);
end