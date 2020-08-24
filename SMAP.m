%  DESCRIPTION:   SMAP: Superresolution Microscopy Analysis Platform
%  COPYRIGHT:     Jonas Ries, 2020
%  LICENSE:       GPLv3
%  AUTHOR:        Jonas Ries, EMBL Heidelberg, ries@embl.de 27.03.2020
%                 www.rieslab.de, www.github.com/jries/SMAP
%  Please cite as: 
%  Ries, J. SMAP: a modular super-resolution microscopy analysis platform for SMLM data. Nat Methods (2020). https://doi.org/10.1038/s41592-020-0938-1
%  Please also cite the references for the plugins you use (as mentioned in
%  the plugin info).

%This is a script to start SMAP.

dirlist=genpath('shared');
addpath(dirlist)
if exist('g','var')
    warning('off','MATLAB:class:DestructorError')
    delete(g)
    warning('on','MATLAB:class:DestructorError')
end

    g=gui.GuiMainSMAP;g.makeGui;       

% display git status
if ~isdeployed
    [status,message]=system('git status');
    if status==0
        ind=find(message==10);
        disp(['git: ' message(1:ind(min(4,length(ind))))]);
    end
end