function writehelpfile(helpfile,guidef)
if isempty(helpfile)
    return
end
if ~contains(helpfile, filesep) %not full path
    outdir='/Users/ries/Documents/MATLAB/SMAPhelp/';
    helpfile=[outdir helpfile];
end
fid=fopen(helpfile,'w');
if isfield(guidef.plugininfo,'description')
    description=guidef.plugininfo.description;
else
    description = '';
end
fprintf(fid,description);
fprintf(fid,'\n');
fprintf(fid,'\n');
fprintf(fid,'gui:Parameters:');
fprintf(fid,'\n');

fn=fieldnames(guidef);
fn=setdiff(fn,'plugininfo');
for k=1:length(fn)
    if ~isfield(guidef.(fn{k}),'object')
        continue     
    end
    fprintf(fid,'gui:');
    fprintf(fid,fn{k});
    fprintf(fid,' ');
    
    fh=guidef.(fn{k}).object;
    if isfield(fh,'Tooltip')
        fprintf(fid,fh.Tooltip);
    elseif isfield(fh,'TooltipString')
        fprintf(fid,fh.TooltipString);
    elseif isfield(guidef.(fn{k}),'Tooltip')
        fprintf(fid,guidef.(fn{k}).Tooltip);
    elseif isfield(guidef.(fn{k}),'TooltipString')
        fprintf(fid,guidef.(fn{k}).TooltipString);
    end
    fprintf(fid,'\n');
end
fclose(fid)
end