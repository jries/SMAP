function par=addFields(par1,par2)
    par=par1;
    allFields=fieldnames(par2);
    for k=1:length(allFields)
        par.(allFields{k})=par2.(allFields{k});
    end
end