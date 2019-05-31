%This is a script to start SMAP.

dirlist=genpath('shared');
addpath(dirlist)
if exist('g','var')
    
    delete(g)
end

% try
    g=gui.GuiMainSMAP;g.makeGui;       
% catch err
%     disp('ESMAPrror making the GUI. Try deleting plugins/plugin.m and the settings/temp directory.')
%     err.rethrow
% end

% display git status
if ~isdeployed
    [status,message]=system('git status');
    if status==0
        ind=find(message==10);
        disp(['git: ' message(1:ind(min(4,length(ind))))]);
    end
end
