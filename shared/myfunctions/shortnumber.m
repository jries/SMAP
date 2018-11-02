function outstr=shortnumber(in)
if abs(in)<1000
    outstr=num2str(in);
elseif abs(in)<1e6
    outstr=[num2str(in/1000,3) 'k'];
else
    outstr=[num2str(in/1e6,3) 'M'];
end