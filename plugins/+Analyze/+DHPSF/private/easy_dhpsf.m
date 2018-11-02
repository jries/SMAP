% Copyright (c)2013, The Board of Trustees of The Leland Stanford Junior
% University. All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
% Redistributions of source code must retain the above copyright notice, 
% this list of conditions and the following disclaimer.
% Redistributions in binary form must reproduce the above copyright notice, 
% this list of conditions and the following disclaimer in the documentation 
% and/or other materials provided with the distribution.
% Neither the name of the Leland Stanford Junior University nor the names 
% of its contributors may be used to endorse or promote products derived 
% from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function easy_dhpsf()
% easy_dhpsf allows scientific users to extract single-molecule
% localizations in 3D when using the double-helix point spread function
% widefield microscope. Tiff stacks of images are analyzed using template
% matching followed by double-Gaussian fitting to extract estimates of the
% molecule positions.

%% global file locations
% file where data is saved
projFile = '';

% all saveable data is originally organized in a structure 's'
% if a channel is selected, then the 's' data is duplicated to that
% channel.
% list of TIFs with raw SM data
s.smacmRawFile = {};
s.smacmRawPath = {};
% list of corresponding files with dark counts
s.smacmDarkFile = {};
% (optional) sequence file with shutter state, z position for each frame
s.smacmSifFile = {};
s.smacmSifPath = {};
% processed SMACM localizations (corresponding to each raw TIF)
s.smacmLocFile = {};
% concatenated SMACM localizations (fid corrected if available)
s.smacmFullLocFile = '';
% DHPSF calibration data, also points to where templates are saved
s.calFilePrefix = '';
% DHPSF fiducial tracking data location
s.fidFilePrefix = {};
% DHPSF SM fit data location
s.fitFilePrefix = {};

%% global fitting parameters
% status of fitting project: 
% cal, fid, thresh, template match, saved output, saved project
s.projStatus = false(1,5);
% calibration bead images
s.templateImgs = [];
% number of calibration beads to choose from
s.numCalBeads = 1;
% selected calibration bead
s.calBeadIdx = 1;
% use fiducials?
s.useFids = false;
% list of template index numbers, corresponds to first dimension of
% templateImgs[]
s.templateIdxs = [];
% list of template thresholds (columns) for each raw smacm file (rows)
s.templateThreshs = [];
% selected file for specifying template matching thresholds
s.threshFileSelIdx = 1;
% selected template
s.templateSelIdx = 1;
% locations of peaks within templates
s.templateLocs = zeros(6,5);
% region of interest inside raw TIF to process
s.smacmRawROI = [0 0 0 0];
% EM gain setting used when acquiring SMACM data
s.smacmEMGain = 300;
% photons/count, camera setting, global to all modules
s.conversionGain = 26.93; %8A % 24.7; % 8B
% imaging system property, global to all modules
s.nmPerPixel = 125.78; % 8A% 160; %8B
% channel identifier
s.channel = '0';
% [minWidth maxWidth] of the two spots of the DHPSF, units of pixels
% relative to the original value of 160 nm / pix
s.sigmaBounds = [1.0 1.5]*160/s.nmPerPixel;
% [minSpacing maxSpacing] between the two spots of the DHPSF, units of
% pixels relative to the original value of 160 nm / pix
s.lobeDistBounds = [3.5 7]*160/s.nmPerPixel;   %[3.5 10]*160/s.nmPerPixel;
% half-width of box to extract when fitting DHPSF images, units of integer pixels
s.boxRadius = round(7*160/s.nmPerPixel);
% smoothing filter width for identifying DHPSF SMs, units of pixels
s.gaussianFilterSigma = 1.5*160/s.nmPerPixel;
% minimum lateral distance between identified SMs, units of pixels
s.minDistBetweenSMs = 7.5*160/s.nmPerPixel;
%
r = s;
g = s;
%% GUI parameters
figSize = [440 610];
figMargin = 20;
panelMargin = 10;
buttonWidth = 90;

%%  Initialize and hide the GUI as it is being constructed.
hfig = figure('Visible','off','Position',[1,1,figSize],...
    'MenuBar','None','ToolBar','none',...
    'DefaultUIPanelUnits','pixels','DefaultAxesUnits','pixels',...
    'DefaultUIPanelFontName','SegoeUI','DefaultUIPanelFontSize',10,...
    'DefaultUIControlFontName','SegoeUI','DefaultUIControlFontSize',9);

%% Construct the components.
htb = uitoolbar(hfig);
% upper master controls
hbuttonNew = uipushtool(htb,'CData',f_iconRead(fullfile(matlabroot,...
    'toolbox','matlab','icons','file_new.png')),...
    'TooltipString','New project',...
    'ClickedCallback',{@buttonNewProj_Callback});
hbuttonLoad = uipushtool(htb,'CData',f_iconRead(fullfile(matlabroot,...
    'toolbox','matlab','icons','file_open.png')),...
    'TooltipString','Load project',...
    'ClickedCallback',{@buttonLoadProj_Callback});
hbuttonSave = uipushtool(htb,'CData',f_iconRead(fullfile(matlabroot,...
    'toolbox','matlab','icons','file_save.png')),...
    'TooltipString','Save project',...
    'ClickedCallback',{@buttonSaveProj_Callback});
htextProjStatus = uicontrol('Style','text',...
    'Position',[figMargin,575,figSize(1)-2*figMargin,15]);
% individual panels for each module
hpanelSetup = uipanel('Title','Setup',...
    'Position',[figMargin 490 figSize(1)-2*figMargin 65]);
hpanelCal = uipanel('Title','Calibrate DHPSF',...
    'Position',[figMargin 405 figSize(1)-2*figMargin 65]);
hpanelThresh = uipanel('Title','Calibrate SM identification',...
    'Position',[figMargin 275 figSize(1)-2*figMargin 110]);
hpanelFid = uipanel('Title','Track fiduciaries',...
    'Position',[figMargin 190 figSize(1)-2*figMargin 65]);
hpanelFit = uipanel('Title','Localize DHPSF SMs',...
    'Position',[figMargin 105 figSize(1)-2*figMargin 65]);
hpanelOut = uipanel('Title','Output DHPSF SM localizations',...
    'Position',[figMargin figMargin figSize(1)-2*figMargin 65]);
% text boxes for statuses
htextCalStatus = uicontrol('Parent',hpanelCal,'Style','text',...
    'Position',[figSize(1)-120,panelMargin+5,70,15]);
htextFidStatus = uicontrol('Parent',hpanelFid,'Style','text',...
    'Position',[figSize(1)-120,panelMargin+5,70,15]);
htextThreshStatus = uicontrol('Parent',hpanelThresh,'Style','text',...
    'Position',[figSize(1)-120,panelMargin+5,70,15]);
htextFitStatus = uicontrol('Parent',hpanelFit,'Style','text',...
    'Position',[figSize(1)-120,panelMargin+5,70,15]);
% setup controls
% channel selection - see also the popupSetupChannel_Callback function
% htextSetupChannel = uicontrol('Parent',hpanelSetup,'Style','text',...
%     'String','Channel:',...
%     'Position',[panelMargin,panelMargin,60,30],...
%     'HorizontalAlignment','left');
% hpopupSetupChannel = uicontrol('Parent',hpanelSetup,'Style','popupmenu',...
%     'Enable','off','Position',[panelMargin+60,panelMargin,40,30],...
%     'String','0|R|G',...
%     'Callback',{@popupSetupChannel_Callback});
htextSetupConv = uicontrol('Parent',hpanelSetup,'Style','text',...
    'String','Conversion gain:',...
    'Position',[panelMargin*2+100,panelMargin,70,30],...
    'HorizontalAlignment','left');
heditSetupConv = uicontrol('Parent',hpanelSetup,'Style','edit',...
    'String','26.93',...
    'Position',[panelMargin*2+170,panelMargin,50,20],...
    'HorizontalAlignment','left',...
    'Callback',{@editSetupConv_Callback});
htextSetupPixSize = uicontrol('Parent',hpanelSetup,'Style','text',...
    'String','Pixel size (nm):',...
    'Position',[panelMargin*3+220,panelMargin,60,30],...
    'HorizontalAlignment','left');
heditSetupPixSize = uicontrol('Parent',hpanelSetup,'Style','edit',...
    'String','125.78',...
    'Position',[panelMargin*3+280,panelMargin,50,20],...
    'HorizontalAlignment','left',...
    'Callback',{@editSetupPixSize_Callback});
% calibration controls
hbuttonCalRun = uicontrol('Parent',hpanelCal,'Style','pushbutton',...
    'String','Run','Position',[panelMargin,panelMargin,buttonWidth,25],...
    'Callback',{@buttonCalRun_Callback});
htextCalSel = uicontrol('Parent',hpanelCal,'Style','text',...
    'String','Calibration bead:',...
    'Position',[panelMargin*2+buttonWidth,panelMargin,70,30],...
    'HorizontalAlignment','left');
hpopupCalSel = uicontrol('Parent',hpanelCal,'Style','popupmenu',...
    'Position',[panelMargin*2+buttonWidth+70,panelMargin,50,30],...
    'Callback',{@popupCalSel_Callback});
haxesCalImg = axes('parent',hpanelCal,'position',[panelMargin*2+buttonWidth+130,panelMargin,45,45]);
% fiducial fitting controls
hbuttonFidRun = uicontrol('Parent',hpanelFid,'Style','pushbutton',...
    'String','Run','Position',[panelMargin,panelMargin,buttonWidth,25],...
    'Callback',{@buttonFidRun_Callback});
hcheckFidUse = uicontrol('Parent',hpanelFid,'Style','checkbox',...
    'String','Use fiduciaries','Min',0,'Max',1,...
    'Position',[panelMargin*2+buttonWidth,panelMargin,110,25],...
    'Callback',{@checkFidUse_Callback});
% thresholding controls
hbuttonThreshRun = uicontrol('Parent',hpanelThresh,'Style','pushbutton',...
    'String','Run','Position',[panelMargin,panelMargin+50,buttonWidth,25],...
    'Callback',{@buttonThreshRun_Callback});
htextThreshFileSel = uicontrol('Parent',hpanelThresh,'Style','text',...
    'String','Filename:',...
    'Position',[panelMargin*2+buttonWidth,panelMargin+50,70,30],...
    'HorizontalAlignment','left');
hpopupThreshFileSel = uicontrol('Parent',hpanelThresh,'Style','popupmenu',...
    'Position',[panelMargin*2+buttonWidth+70,panelMargin+50,210,30],...
    'Callback',{@popupThreshFileSel_Callback});
htextThreshSel = uicontrol('Parent',hpanelThresh,'Style','text',...
    'String','Template:',...
    'Position',[panelMargin*2+buttonWidth,panelMargin+20,70,30],...
    'HorizontalAlignment','left');
hpopupThreshSel = uicontrol('Parent',hpanelThresh,'Style','popupmenu',...
    'Position',[panelMargin*2+buttonWidth+70,panelMargin+20,50,30],...
    'Callback',{@popupThreshSel_Callback});
htextThreshVal = uicontrol('Parent',hpanelThresh,'Style','text',...
    'String','Threshold:',...
    'Position',[panelMargin*2+buttonWidth,panelMargin,70,20],...
    'HorizontalAlignment','left');
heditThreshVal = uicontrol('Parent',hpanelThresh,'Style','edit',...
    'String','0.000',...
    'Position',[panelMargin*2+buttonWidth+70,panelMargin,50,20],...
    'HorizontalAlignment','left',...
    'Callback',{@editThreshVal_Callback});
haxesThreshImg = axes('parent',hpanelThresh,'position',[panelMargin*2+buttonWidth+130,panelMargin,45,45]);
% SM fitting controls
hbuttonFitRun = uicontrol('Parent',hpanelFit,'Style','pushbutton',...
    'String','Run','Position',[panelMargin,panelMargin,buttonWidth,25],...
    'Callback',{@buttonFitRun_Callback});
hbuttonFitDebug = uicontrol('Parent',hpanelFit,'Style','pushbutton',...
    'String','Debug','Position',[panelMargin*2+buttonWidth,panelMargin,buttonWidth,25],...
    'Callback',{@buttonFitDebug_Callback});
% data output controls
hbuttonOutExport = uicontrol('Parent',hpanelOut,'Style','pushbutton',...
    'String','Export to csv','Position',[panelMargin,panelMargin,buttonWidth,25],...
    'Callback',{@buttonOutExport_Callback});
hbuttonOutScatter = uicontrol('Parent',hpanelOut,'Style','pushbutton',...
    'String','3D scatterplot',...
    'Position',[panelMargin*2+buttonWidth,panelMargin,buttonWidth,25],...
    'Callback',{@buttonOutScatter_Callback});
hbuttonOutHist = uicontrol('Parent',hpanelOut,'Style','pushbutton',...
    'String','2D histogram','Position',...
    [panelMargin*3+buttonWidth*2,panelMargin,buttonWidth,25],...
    'Callback',{@buttonOutHist_Callback});
% align([hpanelCal,hpanelFid,hpanelThresh,hpanelFit,hpanelOut],'Center','None');
newProj;

%% Initialize the GUI.
% Change units to normalized so components resize automatically.
set([hfig,htextProjStatus,...
    hpanelSetup,hpanelCal,hpanelFid,hpanelThresh,hpanelFit,hpanelOut,...
    htextCalStatus,htextFidStatus,htextThreshStatus,htextFitStatus,...
    hbuttonCalRun,hbuttonFidRun,hbuttonThreshRun,hbuttonFitRun,hbuttonFitDebug...
    hbuttonOutExport,hbuttonOutScatter,hbuttonOutHist,...
    htextSetupConv,heditSetupConv,... %     htextSetupChannel,hpopupSetupChannel,...
    htextSetupPixSize, heditSetupPixSize,htextCalSel,hpopupCalSel,...
    haxesCalImg,hcheckFidUse,...
    htextThreshFileSel,hpopupThreshFileSel,...
    htextThreshSel,hpopupThreshSel,htextThreshVal,heditThreshVal,haxesThreshImg...
    ],'Units','normalized');

% Assign the GUI a name to appear in the window title.
set(hfig,'Name','Easy-DHPSF')
% initialize axes
axis(haxesCalImg,'off');
axis(haxesThreshImg,'off');
% Move the GUI to the center of the screen.
movegui(hfig,'center')
% Make the GUI visible.
set(hfig,'Visible','on');

    %% helper functions for setting the state of the program
    function newProj
        % confirm new project if current project is unsaved
        if s.projStatus(5) == 0 && any(s.projStatus(1:5))
            prompt = {'The current project is not saved. Continue?'};
            questiondialog = questdlg(prompt,'Confirm','Yes','No','Yes');
            % Handle response
            switch questiondialog
                case 'Yes'
                case 'No'
                    return
            end
        end
        projFile = '';
        s.smacmRawFile = {};
        s.smacmDarkFile = {};
        s.smacmSifFile = {};
        s.smacmLocFile = {};
        s.smacmFullLocFile = '';
        s.calFilePrefix = '';
        s.fidFilePrefix = {};
        s.fitFilePrefix = {};
        s.projStatus = false(1,5);
        s.templateImgs = [];
        s.numCalBeads = 1;
        s.calBeadIdx = 1;
        s.useFids = false;
        s.templateIdxs = [];
        s.templateThreshs = [];
        s.templateLocs = zeros(6,5);
        s.threshFileSelIdx = 1;
        s.templateSelIdx = 1;
        s.smacmRawROI = [0 0 0 0];
        s.channel = '0';
        s.smacmEMGain = 300;
        r = s;
        g = s;
        updateGUI;
    end
    function loadProj
        % confirm new project if current project is unsaved
        if s.projStatus(5) == 0 && any(s.projStatus(1:5))
            prompt = {'The current project is not saved. Continue?'};
            questiondialog = questdlg(prompt,'Confirm','Yes','No','Yes');
            % Handle response
            switch questiondialog
                case 'Yes'
                case 'No'
                    return
            end
        end
        [projFile, projPath] = uigetfile({'*.mat';'*.*'},'Open project MAT file');
        if isequal(projFile,0)
            projFile = '';
            return;
        end
        projFile = [projPath projFile];
        temp = load(projFile); % load into structure due to restriction on static workspaces
        s = temp.s;
        clear temp

        updateGUI;
    end
    function saveProj
        [projFile, projPath] = uiputfile({'*.mat';'*.*'},'Save project MAT file');
        if isequal(projFile,0)
            return;
        end
        projFile = [projPath projFile];
        s.projStatus(5) = true;
        save(projFile,'s');
        
        updateGUI;
    end
    function updateGUI
        set(heditSetupConv,'String',num2str(s.conversionGain));
        set(heditSetupPixSize,'String',num2str(s.nmPerPixel));
        if s.projStatus(1)
            set(htextCalStatus,'String','Complete','BackgroundColor','g');
            set(hbuttonThreshRun,'Enable','on');
            popupChoices = num2str(1:s.numCalBeads,'%g|');
            set(hpopupCalSel,'Enable','on',...
                'String',popupChoices(1:length(popupChoices)-1));
            imagesc(squeeze(s.templateImgs(round(size(s.templateImgs,1)/2),:,:)),'parent',haxesCalImg);
            axis(haxesCalImg,'image');
            axis(haxesCalImg,'off');
            colormap(haxesCalImg,'hot');
        else
            set(htextCalStatus,'String','Incomplete','BackgroundColor','y');
            set(hbuttonFidRun,'Enable','off');
            set(hbuttonThreshRun,'Enable','off');
            set(hpopupCalSel,'Enable','off','String','1','Value',s.calBeadIdx);
            cla(haxesCalImg);
        end
        if s.projStatus(2)
            set(htextThreshStatus,'String','Complete','BackgroundColor','g');
            set(hbuttonFidRun,'Enable','on');
            set(hbuttonFitRun,'Enable','on');
            popupChoices = sprintf('%s|',s.smacmRawFile{:});
            set(hpopupThreshFileSel,'Enable','on',...
                'String',popupChoices(1:length(popupChoices)-1));
            popupChoices = num2str(1:length(s.templateIdxs),'%g|');
            set(hpopupThreshSel,'Enable','on',...
                'String',popupChoices(1:length(popupChoices)-1));
            set(heditThreshVal,'Enable','on','String',num2str(s.templateThreshs(...
                s.threshFileSelIdx,s.templateSelIdx)));
            imagesc(squeeze(s.templateImgs(s.templateIdxs(s.templateSelIdx),:,:)),...
                'parent',haxesThreshImg);
            axis(haxesThreshImg,'image');
            axis(haxesThreshImg,'off');
            colormap(haxesThreshImg,'hot');
        else
            set(htextThreshStatus,'String','Incomplete','BackgroundColor','y');
            set(hbuttonFitRun,'Enable','off');
            set(hpopupThreshFileSel,'Enable','off','String',' ','Value',1);
            set(hpopupThreshSel,'Enable','off','String','1','Value',s.templateSelIdx);
            set(heditThreshVal,'Enable','off','String','0.000','Value',0.000);
            cla(haxesThreshImg);
        end
        if s.projStatus(3)
            set(htextFidStatus,'String','Complete','BackgroundColor','g');
            set(hcheckFidUse,'Enable','on','Value',s.useFids);
        else
            set(htextFidStatus,'String','Incomplete','BackgroundColor','y');
            set(hcheckFidUse,'Enable','off');
        end
        if s.projStatus(4)
            set(htextFitStatus,'String','Complete','BackgroundColor','g');
            set(hbuttonOutExport,'Enable','on');
            set(hbuttonOutScatter,'Enable','on');
            set(hbuttonOutHist,'Enable','on');
            set(hbuttonFitDebug,'Enable','on');
        else
            set(htextFitStatus,'String','Incomplete','BackgroundColor','y');
            set(hbuttonOutExport,'Enable','off');
            set(hbuttonOutScatter,'Enable','off');
            set(hbuttonOutHist,'Enable','off');
            set(hbuttonFitDebug,'Enable','off');
        end
        if s.projStatus(5)
            set(htextProjStatus,'String',projFile,...
                'TooltipString',projFile,...
                'BackgroundColor',get(hfig,'Color'));
        else
            if isempty(projFile)
                set(htextProjStatus,'String','Unsaved project',...
                    'TooltipString','','BackgroundColor','y');
            else
                set(htextProjStatus,'String',['*' projFile],...
                    'TooltipString',projFile,'BackgroundColor','y');
            end
        end
    end

%% GUI function callbacks

    % master controls
    function buttonNewProj_Callback(~,~)
        newProj;
    end
    function buttonLoadProj_Callback(~,~)
        loadProj;
    end
    function buttonSaveProj_Callback(~,~)
        saveProj;
    end

    % setup controls
%     function popupSetupChannel_Callback(source,~) 
%         s.channel = get(source,'Value');
%         updateGUI;
%     end
    function editSetupConv_Callback(source,~)
        s.conversionGain = str2double(get(source,'String'));
        s.projStatus(5) = false;
        updateGUI;
    end
    function editSetupPixSize_Callback(source,~)
        s.nmPerPixel = str2double(get(source,'String'));
        % update all dependent parameters
        s.sigmaBounds = [1.0 1.5]*160/s.nmPerPixel;
        s.lobeDistBounds = [3.5 10]*160/s.nmPerPixel;
        s.boxRadius = round(7*160/s.nmPerPixel);
        s.gaussianFilterSigma = 1.5*160/s.nmPerPixel;
        s.minDistBetweenSMs = 7.5*160/s.nmPerPixel;
        s.projStatus(5) = false;
        updateGUI;
    end

    % calibration controls
    function buttonCalRun_Callback(~,~) 
        [s.calFilePrefix, s.numCalBeads] = ...
            f_calDHPSF(s.conversionGain,s.nmPerPixel,s.boxRadius,'0',s.sigmaBounds,s.lobeDistBounds);
        s.projStatus(1) = true;
        s.projStatus(5) = false;
        temp=load([s.calFilePrefix 'bead ' num2str(s.calBeadIdx) ' templates.mat'],'template');
        s.templateImgs = temp.template;
        clear temp;
        
        s.projStatus(5) = false;
        updateGUI;
    end
    function popupCalSel_Callback(source,~) 
        s.calBeadIdx = get(source,'Value');
        temp=load([s.calFilePrefix 'bead ' num2str(s.calBeadIdx) ' templates.mat'],'template');
        s.templateImgs = temp.template;
        clear temp;
        
        s.projStatus(5) = false;
        updateGUI;
    end

    % threshold controls
    function buttonThreshRun_Callback(~,~) 
        [s.templateIdxs,s.smacmRawROI,s.smacmRawFile, s.smacmRawPath, ...
             s.smacmDarkFile, s.smacmSifFile, s.smacmSifPath, s.smacmEMGain,...
             s.templateLocs] = f_calSMidentification(...
             [s.calFilePrefix 'calibration.mat'],...
             s.calBeadIdx, [s.calFilePrefix 'bead ' num2str(s.calBeadIdx) ' templates.mat'],...
            s.boxRadius,s.channel,s.sigmaBounds,s.gaussianFilterSigma,s.minDistBetweenSMs);
        s.templateThreshs = zeros(length(s.smacmRawFile),length(s.templateIdxs));
        
        s.projStatus(2) = true;
        s.projStatus(5) = false;
        updateGUI;     
    end
    function popupThreshFileSel_Callback(source,~) 
        s.threshFileSelIdx = get(source,'Value');
        updateGUI;
    end
    function popupThreshSel_Callback(source,~) 
        s.templateSelIdx = get(source,'Value');
        updateGUI;
    end
    function editThreshVal_Callback(source,~)
        s.templateThreshs(s.threshFileSelIdx,s.templateSelIdx)...
            = str2double(get(source,'String'));
        
        s.projStatus(5) = false;
        updateGUI;
    end

    % fiducial controls: gives the location of the raw fid fit files as a cell array
    function buttonFidRun_Callback(~,~) 
        [s.fidFilePrefix] = f_trackFiducials(...
            s.smacmRawFile, s.smacmRawPath, [s.calFilePrefix 'calibration.mat'],s.calBeadIdx,...
            [s.calFilePrefix 'bead ' num2str(s.calBeadIdx) ' templates.mat'],...
            s.templateIdxs, s.smacmDarkFile, s.smacmSifFile, s.smacmSifPath, s.boxRadius, ...
            s.channel, s.gaussianFilterSigma,s.minDistBetweenSMs,...
            s.lobeDistBounds,s.conversionGain,s.nmPerPixel,s.smacmEMGain,s.templateLocs,s.sigmaBounds);
        
        s.projStatus(3) = true;
        s.projStatus(5) = false;
        updateGUI;
    end
    function checkFidUse_Callback(source,~) 
        s.useFids = logical(get(source,'Value'));
        updateGUI;
    end
    
    % SM fitting controls
    function buttonFitRun_Callback(~,~)
        if  any(s.templateThreshs==0)
            msgbox(['One or more of the chosen template thresholds equals 0. '...
                    'Please define all thresholds before template matching.'],...
                    'Define thresholds', 'warn');
            return
        end
        [s.fitFilePrefix] = f_fitSMs(s.smacmRawFile, s.smacmRawPath, ...
            [s.calFilePrefix 'calibration.mat'],s.calBeadIdx,...
            [s.calFilePrefix 'bead ' num2str(s.calBeadIdx) ' templates.mat'],...
            s.templateIdxs,s.templateThreshs/10000, s.smacmDarkFile, s.smacmSifFile, s.smacmSifPath, s.boxRadius, ...
            s.channel,s.sigmaBounds, s.gaussianFilterSigma,s.minDistBetweenSMs,...
            s.lobeDistBounds,s.conversionGain,s.nmPerPixel,s.smacmEMGain,s.templateLocs,s.smacmRawROI); 
        
        s.projStatus(4) = true;
        s.projStatus(5) = false;
        updateGUI;
    end
    function buttonFitDebug_Callback(~,~) 
        totalPSFfits = f_concatSMfits(s.fitFilePrefix,s.useFids,s.fidFilePrefix);
        f_debugMoleculeFits(totalPSFfits)
        updateGUI;
    end
    % output controls
    function buttonOutExport_Callback(~,~) 
        [csvFile, csvPath] = uiputfile({'*.csv';'*.txt';'*.*'},'Save localizations as comma-separated file');
        if isequal(csvFile,0)
            return;
        end
        % makes sure that output has an extension
        if csvFile(end-3) ~= '.'
            csvFile = [csvFile '.csv'];
        end
        % open a file for writing
        [fid,message] = fopen([csvPath csvFile], 'w');
        if ~isempty(message)
            error([csvPath csvFile ': ' message]);
            %return;
        end
        % print a title, followed by a blank line
        fprintf(fid, ['frame num,molecule num,x (nm),y (nm),z (nm),' ...
            'x fid-corrected (nm),y fid-corrected (nm), z fid-corrected (nm),'...
            'photons detected,mean background photons\n']);
        fclose(fid);
        totalPSFfits = f_concatSMfits(s.fitFilePrefix,s.useFids,s.fidFilePrefix);
        goodFit = totalPSFfits(:,17)>0;
        dlmwrite([csvPath csvFile],totalPSFfits(goodFit,[1 2 25 26 27 28 29 30 21 15]),...
            '-append');
        disp('Output file written successfully.');
        updateGUI;
    end
    function buttonOutScatter_Callback(~,~) 
        totalPSFfits = f_concatSMfits(s.fitFilePrefix,s.useFids,s.fidFilePrefix);
        f_scatter3(totalPSFfits,s.useFids);
        
        updateGUI;
    end
    function buttonOutHist_Callback(~,~) 
        totalPSFfits = f_concatSMfits(s.fitFilePrefix,s.useFids,s.fidFilePrefix);
        f_hist2(totalPSFfits,s.useFids);
        
        updateGUI;
    end


end