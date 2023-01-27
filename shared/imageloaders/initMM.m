function initMM(obj)
global SMAP_globalsettings
% dirs={'ij.jar'
% 'plugins/Micro-Manager/MMAcqEngine.jar'
% 'plugins/Micro-Manager/MMCoreJ.jar'
% 'plugins/Micro-Manager/MMJ_.jar'
% 'plugins/Micro-Manager/clojure.jar'
% 'plugins/Micro-Manager/bsh-2.0b4.jar'
% 'plugins/Micro-Manager/swingx-0.9.5.jar'
% 'plugins/Micro-Manager/swing-layout-1.0.4.jar'
% 'plugins/Micro-Manager/commons-math-2.0.jar'
%  'plugins/Micro-Manager/ome-xml.jar'
%  'plugins/Micro-Manager/scifio.jar'
%  'plugins/Micro-Manager/guava-17.0.jar'
%  'plugins/Micro-Manager/loci-common.jar'
%  'plugins/Micro-Manager/slf4j-api-1.7.1.jar'};

%     if ispc
%         MMpath='C:/Program Files/Fiji/scripts';
%     else
%         MMpath='/Applications/Fiji.app/scripts';
%     end

if isempty(SMAP_globalsettings)
    disp('Micro-manager java path not added, as imageloader was not called from SMAP. add manually to javaclasspath');
end
try    

MMpath=obj.getGlobalSetting('MMpath'); 
catch
    MMpath=SMAP_globalsettings.MMpath.object.String;
end

if ~exist(MMpath,'dir')       
    errordlg('cannot find Micro-Manager, please select Micro-Manager directory in menu SMAP/Preferences/Directotries2...')
    return
end


% for k=1:length(dirs)
%     dirs{k}=[MMpath filesep strrep(dirs{k},'/',filesep)];
% end
plugindir=[MMpath filesep 'plugins' filesep 'Micro-Manager' filesep];
allf=dir([plugindir '*.jar']);
dirs={allf(:).name};

for k=1:length(dirs)
    dirs{k}=[MMpath filesep 'plugins' filesep 'Micro-Manager' filesep strrep(dirs{k},'/',filesep)];
end

mmjar=[MMpath filesep 'ij.jar']; 
if ~exist(mmjar,'file')
     disp('Micro-manager V2.0 has a bug. Please install V1.4 instead.')
    mmjar=[MMpath filesep 'ImageJ.app' filesep 'Contents' filesep 'Java' filesep 'ij.jar']; 
end
if ~exist(mmjar,'file')
    disp('ij.jar not found. imageloaderMM line 159')
%     mmjar=[MMpath filesep 'ImageJ.app' filesep 'Contents' filesep 'Java' filesep 'ij.jar']; 
end

dirs{end+1}=  mmjar; 
jp=javaclasspath;
diradd=dirs(~ismember(dirs,jp));
if ~isempty(diradd)
javaaddpath(diradd);
end

end