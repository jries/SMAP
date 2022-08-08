function  savesml(locData,file,p,excludesavefields)
if p.saveroi
    [~,indgroi]=locData.getloc('xnm','position','roi');  
    indgl=false(size(indgroi));
    for k=1:p.numberOfLayers
        if p.sr_layerson(k)
           [locs,indgh]=locData.getloc('inlayeru','layer',k); 
%            if length(indgh)~=length(indgl)
%                disp('save visible works only for ungrouped data')
%            end
           indgl=indgl|locs.inlayeru{1};
       
        end
    end
    indg=indgl&indgroi;
    
else
%     indg=true(size(locData.loc.xnm));
    indg=[];
end

if isfield(p,'savefile') && p.savefile
    filenumber=p.dataselect.Value;
%     indg=indg&locData.loc.filenumber==filenumber;
elseif isfield(p,'saveSepFile') && p.saveSepFile
    lastfile=locData.files.filenumberEnd;
else
    filenumber=[];
end

% if all(indg)%save all
%     indg=[];
% end



% if nargin>3
% locData.loc=rmfield(locData.loc,excludesavefields);
% end
% 
% if isfield(p,'savefile') && p.savefile
%     saveloc.file=saveloc.file(filenumber);
%     saveloc.history=saveloc.history(filenumber);
%     saveloc.loc.filenumber=ones(size(locData.loc.filenumber));
% end

if isfield(p,'saveSepFile') && p.saveSepFile % Yu-Le added
    for k = 1:lastfile
        filenumber = k;
        [path,fileName] = fileparts(file);
        lSameFolder__ = startsWith(fileName,'__');
        lOwnFolder__ = startsWith(fileName,'own__');
        fileName = regexprep(fileName, '^(|own)__', '_');
        if endsWith(fileName,'_se')
            fileName = regexprep(fileName, '(_{1,2}se)$', '');
            char_fromTheLast = 6;
        else
            fileName = regexprep(fileName, '(_{1,2}sml)$', '');
            char_fromTheLast = 7;
        end
        if lSameFolder__||lOwnFolder__
            % different names, same suffix
            [ownPath,oriName] = fileparts(locData.files.file(filenumber).name);
            newName = [oriName(1:end-4) fileName file(end-char_fromTheLast:end)]; % take out '_sml'
            if lOwnFolder__
                oneFile = [ownPath filesep newName];
            else
                oneFile = [path filesep newName];
            end
        else
            % same name, different suffix
            oneFile = [file(1:end-8) '_' num2str(filenumber) file(end-7:end)];
        end
        switch p.pluginpath{end}
            case 'SMLMsaver'
                saveloc=locData.savelocs(oneFile,indg,[],[],excludesavefields,filenumber); % BETA , maybe problematic with more than 1 file: this will save only displayed loicalizations
            case 'generalSeSaver'
                saveloc=locData.saveSE(oneFile,indg,[],[],excludesavefields,filenumber);
        end
        
    end
else
    switch p.pluginpath{end}
        case 'SMLMsaver'
            saveloc=locData.savelocs(file,indg,[],[],excludesavefields,filenumber); % BETA , maybe problematic with more than 1 file: this will save only displayed loicalizations
        case 'generalSeSaver'
            saveloc=locData.saveSE(file,indg,[],[],excludesavefields,filenumber);
    end
end


% rg=p.mainGui; 
% parameters=rg.saveParameters;
% fileformat.name='sml';
% out=struct('saveloc',saveloc,'fileformat',fileformat,'parameters',parameters);
% if isfield(locData.files.file(1),'transformation')
%     for k=length(locData.files.file):-1:1
%         if ~isempty(locData.files.file(k).transformation)
%             out.transformation=locData.files.file(k).transformation;
%             break
%         end
%     end
%     
% end
% 
% if locData.getGlobalSetting('saveas73') %now I use this to save with the old saver: large files are saved as v7.3, not as small parts.
%     v=saverightversion(file,out,'-v7');
%     disp(['saved as version ' v])
% %     save(file,'saveloc','fileformat','parameters','-v7.3');
% else
%     savematparts(file,out,{'saveloc','loc'});
% 
% end

% save(file,'saveloc','fileformat','parameters','-v7.3');
end