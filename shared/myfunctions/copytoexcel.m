function str=copytoexcel(in,flip)
if isstruct(in)
    in=struct2table(in);
elseif iscell(in)
    in=cell2table(in);
elseif ~istable(in)
    in=array2table(in);
end

if nargin>1 && strcmp(flip,'flip')
    in2=array2table(table2array(in)');
    in2.Properties.RowNames=in.Properties.VariableNames;
else
    in2=in;
end

rownames=in2.Properties.RowNames;
varnames=in2.Properties.VariableNames;
var=table2cell(in2);%in2.Variables;

str=sprintf('\t');
str = [str sprintf('%s\t', varnames{:})];
newline = sprintf('\n');
str=[str newline];
for i = 1:size(var,1) 
  row = [rownames{i} sprintf('\t') sprintf('%f\t', var{i,:})];
  row(end) = newline;
  str = [str row]; 
end
str(end)=[];
clipboard('copy',str)