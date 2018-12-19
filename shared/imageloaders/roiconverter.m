function roi=roiconverter(roiin,roimode)
roi=roiin;
if isempty(roimode)
    return
end
switch roimode
    case 'Evolve-Normal'
        
    case 'Evolve-EM'
        roi(1)=512-roi(1)-roi(3)+2; %somewhere this 2 comes from... 
    case 'Andor'
        roi(1)=512-roi(1)-roi(3); %em mode is reference for roi, gets transformed
    otherwise
end
end