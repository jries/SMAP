function p=makerelativetopwr(p)
smapdir=pwd;
if strcmp(p(1:length(smapdir)),smapdir)
    p=p(length(smapdir)+1:end);
end
if strcmp(p(1),filesep)
    p=p(2:end);
end