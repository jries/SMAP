function  [file,dato]=saveLocalizationsCSV(locData,file,saveroi,numberOfLayers,sr_layerson,indg)
if nargin<3
    saveroi=false;
end


    
if saveroi
    %find if any layer is grouped
    grouped=false;
    nongrouped=false;
    for k=1:numberOfLayers
        if sr_layerson(k)
            grouped=grouped|locData.isgrouped(k);
            nongrouped=nongrouped|(~locData.isgrouped(k));
        end
    end
    if grouped && nongrouped
         warning('save visible works only for either grouped or ungrouped data, not for mixed')
    end
    if grouped
        locext='grouploc';
    else 
        locext='loc';
    end
    
    [~,indg]=locData.getloc('xnm','position','roi','grouping',grouped);  
    for k=1:numberOfLayers
        if sr_layerson(k)
           [~,indgh]=locData.getloc('xnm','layer',k); 
%            if length(indgh)~=length(indg)
%                disp('save visible works only for ungrouped data')
%            end
           indg=indg&indgh;
       
        end
    end
    
else
    grouped=false;
%     indg=[];
end
loc=locData.savelocs([],indg,[],grouped).loc; 
numlocs=length(loc.frame);
dato=zeros(numlocs,6);
dato(:,1)=1:numlocs;
dato(:,2)=loc.frame;
dato(:,3)=loc.xnm;
dato(:,4)=loc.ynm;
if isfield(loc,'znm')
    dato(:,5)=loc.znm;
end
dato(:,6)=loc.phot;
dato(:,7)=loc.bg;
dato(:,8)=loc.locprecnm;

fid=fopen(file,'w');
fprintf(fid, 'ID,frame,xnm,ynm,znm,phot,bg,locprecnm \n');
fclose(fid);
dlmwrite(file,dato, '-append');

if 0
 %save description file
 path=fileparts(file);
fdesc=[path filesep 'file-description.xml'];
fid=fopen(fdesc,'w');
fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?> \n <description> \n');
fprintf(fid,' <firstrow>0</firstrow>\n');
fprintf(fid,' <frame>1</frame>\n');
fprintf(fid,' <xnano>2</xnano>\n');
fprintf(fid,' <ynano>3</ynano>\n');
fprintf(fid,' <znano>4</znano>\n');
fprintf(fid,' <intensity>5</intensity>\n');
fprintf(fid,' <separator>COMMA</separator>\n');
fprintf(fid,' <xshift>0</xshift>\n');
fprintf(fid,' <yshift>0</yshift>\n');
fprintf(fid,' <zshift>0</zshift>\n');
fprintf(fid,'</description>');
fclose(fid);
end
end