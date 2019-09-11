%addpath for compiling (no addpath allowed if deployed)
SMAP
%mainGui
addpath('shared');
addpath(pwd);
% %fit3Dcspline:
% fit3ddir=strrep(pwd,'SMAP','fit3D');
% if exist(fit3ddir,'file')
%     addpath(fit3ddir);
% end
% 
% fit3ddir=strrep(pwd,'SMAP','fit3Dcspline');
% if exist(fit3ddir,'file')
%     addpath(fit3ddir);
% end      


%MM
global SMAP_globalsettings
try    

MMpath=g.getGlobalSetting('MMpath'); 
catch
    MMpath=SMAP_globalsettings.MMpath.object.String;
end
plugindir=[MMpath filesep 'plugins' filesep 'Micro-Manager' filesep];
allf=dir([plugindir '*.jar']);
dirs={allf(:).name};

for k=1:length(dirs)
    dirs{k}=[MMpath filesep 'plugins' filesep 'Micro-Manager' filesep strrep(dirs{k},'/',filesep)];
end

dirs{end+1}=  [MMpath filesep 'ij.jar']; 
jp=javaclasspath;
diradd=dirs(~ismember(dirs,jp));
if ~isempty(diradd)
javaaddpath(diradd);
end

%fitter 4 Pi
fitterpath=[fileparts(g.getPar('maindirectory')) filesep 'ries-private' filesep 'PSF4Pi'];
addpath(fitterpath)

%modelfit
addpath(['..' filesep 'ries-private' filesep 'SMLMModelFitter'], ['..' filesep 'ries-private' filesep 'SMLMModelFitter' filesep 'external'])
%utrack
addpath(genpath('External/u-track/software'))

%fit with vectorial PSF
newpath=strrep(pwd,'SMAP','ries-private');
addpath([newpath filesep 'VectorPSF_Fit']);

% openfiji:

fijipath=g.getGlobalSetting('fijipath');  
addpath(fijipath)
disp('close ImageJ / Fiji')
ImageJ


