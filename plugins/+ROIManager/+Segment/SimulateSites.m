classdef SimulateSites<interfaces.DialogProcessor&interfaces.SEProcessor
    %SimulateSites is a localization based simulation engine for SMAP. It
    %uses as an input a list of localizations, a matlab function that
    %returns coordiantes or an image which defines a 2D structure. It
    %returns simulated localizations to SMAP using a realistic model for the
    %photophysics of the dye. Simulated structures are added to the
    %RoiManager
    properties
        lLocMoFitGUI_loaded = false;    % whether any LocMoFitGUI is loaded.
    end
    methods
        function obj=SimulateSites(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_sitefov','se_cellpixelsize','se_siteroi'};
            obj.history=true;
            obj.guiselector.show=true;
        end
        function initGui(obj)
            setvisibility(obj);
            % added by Yu-Le
            % check through the name of loaded eval plugins and find LocMoFit 
            if ~isempty(obj.locData.SE.processors.eval.guihandles.modules.Data)
                lFitterFound = strfind(obj.locData.SE.processors.eval.guihandles.modules.Data(:,2), 'LocMoFitGUI');
                lFitterFound = [lFitterFound{:}];
                lFitterFound = any(lFitterFound);
            else
                lFitterFound = 0;
            end
            if lFitterFound 
                fitterFound = 'on';
            else
                fitterFound = 'off';
            end
            obj.guihandles.useFitter_button.Visible=fitterFound;
            obj.guihandles.setModPars_button.Visible='off';
        end
        function out=run(obj,p)  
            [locst,possites,parameters]=simulatelocs(p, 1);
            
           if ~p.savez
               locst=rmfield(locst,{'znm','znm_gt'});
           end
           

               labels=num2str(p.labeling_efficiency*100,'%2.0f');
                phots=num2str(p.photons,'%3.0f');
                blinks=num2str(p.blinks,'%3.0f');
                bg=num2str(p.background,'%3.0f');
                lt=num2str(p.lifetime,'%3.0f');
                filename=[p.model.selection '_L' labels 'P' phots 'B' bg 'R' blinks 'L' lt ];
           
           
           obj.locData.addfile(['simulated_' num2str(obj.locData.files.filenumberEnd) '_' filename]);
           obj.locData.files.file(end).info.simulationParameters=obj.getGuiParameters;
           obj.locData.addLocData(locst);
           obj.locData.sort('filenumber','frame');
%            try
           initGuiAfterLoad(obj);
%            catch err
%                err
%            end
           se=obj.locData.SE;
           cell=interfaces.SEsites;
           cell.pos=[mean(locst.xnm) mean(locst.ynm)];
           cell.ID=0;
           cell.info.filenumber=obj.locData.files.filenumberEnd;
           se.addCell(cell);
           for k=1:length(possites)
               thissite=interfaces.SEsites;
               thissite.pos=[possites(k).x possites(k).y];
               thissite.info.cell=cell.ID;
               thissite.info.filenumber=cell.info.filenumber;
               thissite.evaluation.simulatesites=parameters(k);
                % thissite.cellnumber=sitepar.currentcell.number;
        %         thissite.number=sitepar.sitelist.cellnumber+1;
                se.addSite(thissite);
           end 
           
%            obj.setPar('SimulateSitesParameters',p);

            try
           se.currentsite=se.sites(1);
           se.currentcell=se.cells(1);
           se.currentfile=se.files(1);
           se.processors.preview.updateFilelist;
           se.processors.preview.updateCelllist;
           se.processors.preview.updateSitelist; 
            se.processors.preview.nextsite(1)
%            se.processors.plotsite(se.sites(1));
            catch err
                err
            end
           
            if p.savenow.Value==2
                obj.locData.loc.xnm=obj.locData.loc.xnm_gt;
                obj.locData.loc.ynm=obj.locData.loc.ynm_gt;
                obj.locData.loc.znm=obj.locData.loc.znm_gt;
                obj.locData.loc=rmfield(obj.locData.loc,{'xnm_gt','ynm_gt','znm_gt'});
            end
            if p.savenow.Value>1
                lastsml=obj.getPar('lastSMLFile');
                if ~isempty(lastsml)
                    path=fileparts(lastsml);
                else
                    path='';
                end
               [file,path]= uiputfile([path filesep obj.locData.files.file(1).name]);
               if file
                    obj.locData.savelocs([path file])
                    obj.setPar('lastSMLFile',[path file]);
               end
            end
          out=[];
        end
        
        function prepare_LocMoFit(obj, fitter)
            fitter_cp = copy(fitter);  % copy the fitter object to not overwrite it
            fitter_cp.allParsArg.fix = true(size(fitter_cp.allParsArg.fix));
            fitter_cp.rmConvertRules;
            for l = 1:length(fitter_cp.allModelLayer)
                fitter_cp.addPar({90+fitter_cp.allModelLayer(l),{'sim'},{'numOfMol'},0, inf,0,1,{''},0,inf}) % add this parameter to control the number of molecules.
            end
            if isempty(fitter_cp.dimension)
                fitter_cp.dimension = fitter_cp.dataDim;
            end
            obj.setPar('fitter',fitter_cp)
            obj.setPar('fitter_ori',fitter)
            obj.guihandles.setModPars_button.Visible='on';
            obj.guihandles.coordinatefile.String='-- Internal LocMoFit';
        end

        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end



function load_callback(a,b,obj)
f=obj.getSingleGuiParameter('coordinatefile'); 
[f,p]=uigetfile({'*.*';'*.tif';'*.png';'*.csv';'*.txt';'*.mat'},'Choose coordinate file',f);
if ~f
    return
end
obj.setGuiParameters(struct('coordinatefile',[p f]));
setvisibility(obj)
end

function loadpsf_callback(a,b,obj)
f=obj.getSingleGuiParameter('psf_file'); 
[f,p]=uigetfile('*_3dcal.mat','Choose bead calibration psf file',f);
if ~f
    return
end
obj.setGuiParameters(struct('psf_file',[p f]));
end


function useFitter_callback(a,b,obj)
    fig = figure(512);
    clf(fig);
    selectionTable = uitable(fig);
    nameEvalPlugins = obj.locData.SE.processors.eval.guihandles.modules.Data(:,2);
    
    lFitterFound = strfind(nameEvalPlugins, 'LocMoFitGUI');
    for k = 1:length(lFitterFound)
        lFitterFound{k} = ~isempty(lFitterFound{k});
    end
    lFitterFound = [lFitterFound{:}];
    
    nameLocMoFitGUI = nameEvalPlugins(logical(lFitterFound));
    selectionTable.Data = [num2cell(false(size(nameLocMoFitGUI))) nameLocMoFitGUI];
    selectionTable.ColumnEditable = [true false];
    selectionTable.CellEditCallback = {@selectionTable_CECallback};
    selectionTable.Position = [20 50 300 300];
    apply_button = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Apply');
    apply_button.Position = [20 20 60 30];
    apply_button.Callback = {@applySelecedFitter_callback, obj, selectionTable,fig};
end

function selectionTable_CECallback(a,b)
    a.Data(:,1) = num2cell(false);
    a.Data{b.Indices(1),b.Indices(2)} = true;
end

function applySelecedFitter_callback(a,b,obj, selectionTable, fig)
    % Added by Yu-Le
    % set the selected LocMoFit as a parameter of SimulateSites.
    selected = selectionTable.Data(:,1);
    idxSelected = find([selected{:}]);
    nameEvalPlugins = obj.locData.SE.processors.eval.guihandles.modules.Data(:,2);
    
    lFitterFound = strfind(nameEvalPlugins, 'LocMoFitGUI');
    for k = 1:length(lFitterFound)
        lFitterFound{k} = ~isempty(lFitterFound{k});
    end
    idxFitterFound = find([lFitterFound{:}]);
    
    idxSelected_final = idxFitterFound(idxSelected);
    fitter_ori = obj.locData.SE.processors.eval.processors{idxSelected_final}.fitter;
   
    obj.prepare_LocMoFit(fitter_ori);
    close(fig);
end

function setModPars_callback(a,b,obj)
    % Added by Yu-Le
    % Use the function of allParsArg to set the range for simulation.
    
    %% GUI
    fig = obj.getPar('parameter_handle');
    if isempty(fig)||~isgraphics(fig)
        fig = figure('Name','Model parameters');
        obj.setPar('parameter_handle',fig);
    end
    clf(fig);
    guihandles = [];
    guihandles.parArgTable = uitable(fig);
    guihandles.parArgTable.Position = [1 2 3 10];
    
    guihandles.convertTable = uitable(fig);
    guihandles.convertTable.Position = [4.2 2 2 10];
    
    figHeight = fig.Position(4);
    oneLine = 20;
    unitWidth = 50;
    
    
    % Model type
    guihandles.typeOption = uicontrol(fig, 'Style','popupmenu','String',{'Point','Image'},'Value',1);
    guihandles.typeOption.Position = [1 12 0.5 1];
    guihandles.typeOption.Callback = {@typeOption_callback,obj};
    typeOption_callback(guihandles.typeOption,[],obj);
    
    % Final ROI size
    % This is now replaced by depth.
    guihandles.t_FinalROISize = uicontrol(fig, 'Style','text','String','Final ROI size:');
    guihandles.t_FinalROISize.Position = [2 12 1 1];
    
    guihandles.t_Depth = uicontrol(fig, 'Style','text','String','Depth:');
    guihandles.t_Depth.Position = [2 12 1 1];
    
    % get settings
    finalROISize = obj.getPar('finalROISize');
    depth = obj.getPar('depth');
    useDepth = obj.getPar('useDepth');

    % initiate settings if required
    if isempty(useDepth)
        useDepth = 0;
        obj.setPar('useDepth', useDepth);
    end
    if isempty(finalROISize)
        obj.setPar('finalROISize','200')
        finalROISize = '200';
    end
    if isempty(depth)
        obj.setPar('depth','200')
        depth = '200';
    end
    guihandles.finalROISize = uicontrol(fig, 'Style','edit','String',finalROISize);
    guihandles.finalROISize.Position = [3 12 0.5 1];
    guihandles.finalROISize.Callback = {@finalROISize_callback,obj};
    
    guihandles.depth = uicontrol(fig, 'Style','edit','String',depth);
    guihandles.depth.Position = [3 12 0.5 1];
    guihandles.depth.Callback = {@depth_callback,obj};

    if useDepth
        guihandles.t_Depth.Visible = 'on';
        guihandles.depth.Visible = 'on';
        guihandles.t_FinalROISize.Visible = 'off';
        guihandles.finalROISize.Visible = 'off';
    else
        guihandles.t_Depth.Visible = 'off';
        guihandles.depth.Visible = 'off';
        guihandles.t_FinalROISize.Visible = 'on';
        guihandles.finalROISize.Visible = 'on';
    end
    
    guihandles.toggle = uicontrol(fig, 'Style','pushbutton','String','<->');
    guihandles.toggle.Position = [3.5 12 0.5 1];
    guihandles.toggle.Callback = {@toggle_callback,{obj,guihandles}};
    
    %% Data
    % Acquire the LocMoFit obj, and then display parameters based on the allParsArg
    fitter = obj.getPar('fitter');

    parName = fitter.allParsArg.name;
    parType = fitter.allParsArg.type;
    parModel = fitter.allParsArg.model;
    
    % If 'fix' is ticked, the corresponding parameters will be single
    % values. Otherwise a range.
    parUb = fitter.allParsArg.ub;
    parLb = fitter.allParsArg.lb;
    parFix = fitter.allParsArg.fix;
    parVal = cellstr(num2str(fitter.allParsArg.value));
    parRange = cellstr([num2str(parLb) '' num2str(parUb)]);
    parRange = regexprep(parRange,'\s+',' ');
    parVal(~parFix)=parRange(~parFix);
    parVal = regexprep(parVal,'^\s+','');
    
    guihandles.t_parTable = uicontrol(fig, 'Style','text','String','Parameters for simulation:');
    guihandles.t_parTable.Position = [1 1 2 1];
    
    guihandles.t_convertTable = uicontrol(fig, 'Style','text','String','User-defined variables:');
    guihandles.t_convertTable.Position = [4.2 1 2 1];
    
    % Table properties.
    guihandles.parArgTable.Data = [parName parType num2cell(parModel) parVal repmat({''},size(parVal,1),1)];
    guihandles.parArgTable.ColumnName = {'Name','Type','Model','Value','Convert'};
    guihandles.parArgTable.ColumnEditable = [false false false true true];
    guihandles.parArgTable.CellEditCallback = {@parArgTable_CellEditCallback, fitter};
    guihandles.parArgTable.ColumnWidth = {70 40 30 50 100};    
    
    guihandles.convertTable.Data = [];
    guihandles.convertTable.ColumnName = {'Name','Rule'};
    guihandles.convertTable.ColumnEditable = [true true];
    guihandles.convertTable.ColumnWidth = {70 100};
    
    guihandles.addRow = addRowButton(fig,guihandles.convertTable);
    guihandles.addRow.Position = [4.2 12 0.2 0.8];
    
    guihandles.rmRow = rmRowButton(fig,guihandles.convertTable);
    guihandles.rmRow.Position = [4.4 12 0.2 0.8];
    
    guihandles.apply = uicontrol(fig, 'Style','pushbutton','String','Apply','Callback',{@applyConvertRules,guihandles.convertTable,fitter});
    guihandles.apply.Position = [4.6 12 0.5 0.8];
    
    %% Editor
    guihandles.editor = uicontrol(fig, 'Style','pushbutton','String','Editor');
    guihandles.editor.Position = [3.5 1 0.5 1];
    guihandles.editor.Callback = {@editor_callback,{guihandles.parArgTable,guihandles.convertTable,obj}};

    %% save and load
    guihandles.saveParArg = uicontrol(fig, 'Style','pushbutton','String','Save');
    guihandles.saveParArg.Position = [2.5 1 0.5 1];
    guihandles.saveParArg.Callback = {@saveParArg_callback,{guihandles.parArgTable,guihandles.convertTable}};
    
    guihandles.loadParArg = uicontrol(fig, 'Style','pushbutton','String','Load');
    guihandles.loadParArg.Position = [3 1 0.5 1];
    guihandles.loadParArg.Callback = {@loadParArg_callback,{guihandles.parArgTable,guihandles.convertTable}};
    guiStyle(guihandles, fieldnames(guihandles))
end

function editor_callback(a,b,inputs)
%     fig = figure('Name','Editor');
%     text2show = sprintf(['%s\t%s\t%s\t%s\t%s\n'], string(parArgTable.Data'));
%     guihandles.editor = uicontrol(fig, 'Style','edit','String', text2show);
%     guihandles.editor.Max = 2;
%     guihandles.editor.Position = [1 1 5 12];
%     guihandles.editor.HorizontalAlignment = 'left';
%     
%     guihandles.save = uicontrol(fig, 'Style','pushbutton','String', 'save');
%     guihandles.save.Position = [1 13 1 1];
%     guihandles.save.Callback = {@editorSave_callback,guihandles.editor,parArgTable};
%     guiStyle(guihandles,fieldnames(guihandles));
    parArgTable = inputs{1};
    convertTable = inputs{2};
    obj = inputs{3};
    szConvertTable = size(convertTable.Data);
    blank = repmat({''},[szConvertTable(1) 3]);

    % Save parTable to the temp folder
    warning('off', 'MATLAB:MKDIR:DirectoryExists');
    mkdir('.\LocMoFit\temp')
    warning('on', 'MATLAB:MKDIR:DirectoryExists');
    fileName = ['temp_' datestr(datetime, 'yyyymmddHHMMSS') '.csv'];
    saveAs = ['.\LocMoFit\temp\' fileName];
    writecell([parArgTable.Data;[convertTable.Data blank]], saveAs);
    
    % Open it in the default editor of .csv files.
    winopen(saveAs)

    % monitor whether the .csv is changed
    % if so, load the updated values back to the GUI via
    % loadParArg_callback
    file = System.IO.FileSystemWatcher('LocMoFit\temp\');
    file.Filter = fileName;
    file.EnableRaisingEvents = true;
    addlistener(file,'Changed',@(x,y)loadParArg_callback(x,y,inputs(1:2),saveAs));
    obj.setPar('monitorFile',file);
end

function saveParArg_callback(a,b,tables)
    parArgTable = tables{1};
    convertTable = tables{2};
    szConvertTable = size(convertTable.Data);
    blank = repmat({''},[szConvertTable(1) 3]);
    [file,path] = uiputfile({'*.txt';'*.csv'},'Save parameter argument table');
    writecell([parArgTable.Data;[convertTable.Data blank]], [path file]);
end

function loadParArg_callback(a,b,tables,varargin)
    parArgTable = tables{1};
    convertTable = tables{2};
    
    if isempty(varargin)
        [file,path] = uigetfile('*.txt','Load parameter argument table');
    else
        [path,name,ext] = fileparts(varargin{1});
        path = [path filesep];
        file = [name ext];
    end
    opts = delimitedTextImportOptions('NumVariables',5);
    oneTable = readcell([path file],opts);
    oneTable = string(oneTable);
    oneTable = cellstr(oneTable);
    
    % Use the model column to separate the par and convert tables
    model = oneTable(:,3);
    col_convertTable = cellfun(@isempty,model);
    
    parArgTable.Data = oneTable(~col_convertTable,:);
    convertTable.Data = oneTable(col_convertTable,1:2);
    
    callBack = parArgTable.CellEditCallback{1};
    locMoFitObj = parArgTable.CellEditCallback{2};
    editable = [4 5];
    for k = 1:length(editable)
        for l = 1:size(parArgTable.Data,1)
            holder.Indices(1) = l;
            holder.Indices(2) = editable(k);
            holder.NewData = parArgTable.Data{l,editable(k)};
            callBack(parArgTable,holder,locMoFitObj);
        end
    end
end

function editorSave_callback(a,b,editor,parArgTable)
    t = regexprep(string(editor.String)', '\t$', '\t\r\n');
    t = textscan(char(t')','%s %s %s %s %s','delimiter',sprintf('\t'));
    t = [t{:}];
    parArgTable.Data = t;
end

function typeOption_callback(a,b,obj)
    obj.setPar('modelType',a.String{a.Value});
    setvisibility(obj)
end

function finalROISize_callback(a,b,obj)
    obj.setPar('finalROISize',a.String);
end

function depth_callback(a,b,obj)
    obj.setPar('depth',a.String);
end

function toggle_callback(a,b,allObj)
    obj = allObj{1};
    guihandles = allObj{2};
    currentUseDepth = obj.getPar('useDepth');
    if isempty(currentUseDepth)||currentUseDepth ==0
        guihandles.t_Depth.Visible = 'on';
        guihandles.depth.Visible = 'on';
        guihandles.t_FinalROISize.Visible = 'off';
        guihandles.finalROISize.Visible = 'off';
        obj.setPar('useDepth',1);
    else
        guihandles.t_Depth.Visible = 'off';
        guihandles.depth.Visible = 'off';
        guihandles.t_FinalROISize.Visible = 'on';
        guihandles.finalROISize.Visible = 'on';
        obj.setPar('useDepth',0);
    end
end

function parArgTable_CellEditCallback(a,b,obj)
    % Assign the change to the parArgTable
    % This function deal will all editings in the parArgTable. Only column
    % 4 and 5 are editable.
    indEdited = b.Indices(1);
    colId = b.Indices(2);
    elements = strsplit(b.NewData,' ');
    switch colId
        case 4
            % column 4: values
            if length(elements) == 2
                elements = str2double(elements);
                obj.allParsArg.fix(indEdited) = false;
                obj.allParsArg.value(indEdited) = 0;
                obj.allParsArg.lb(indEdited) = elements(1);
                obj.allParsArg.ub(indEdited) = elements(2);
            elseif length(elements) == 1
                obj.allParsArg.fix(indEdited) = true;
                obj.allParsArg.value(indEdited) = str2double(b.NewData);
            else
                warning('The input length is not acceptable. Please assign only 1 (fixed) to 2 (a range) elements.')
            end
        case 5
            % column 5: convert
            parId = ['m' num2str(obj.allParsArg.model(indEdited)), '.',obj.allParsArg.type{indEdited}, '.', obj.allParsArg.name{indEdited}];
            if isempty(b.NewData)
                obj.rmOneConvertRule(parId);
            else
                obj.converter(obj, b.NewData, parId);
            end
    end
end

function applyConvertRules(a,b,hTable,obj)
    data = hTable.Data;
    
    % Remove user defined targets that are not in the current table
    indUsr = startsWith(obj.converterRules.target,'usr_');
    allTarget2set = data(:,1);
    allTarget2set = join([cellstr(repmat('usr_',size(allTarget2set))) allTarget2set],'');
    indNot = ~ismember(obj.converterRules.target, allTarget2set);
    indRm = indUsr&indNot;
    obj.rmOneConvertRule(obj.converterRules.target(indRm))
    
    % Assign convert rules
    for r = 1:size(data,1)
        obj.converter([],data{r,2},allTarget2set{r})
    end
end

function setvisibility(obj)
f=obj.getSingleGuiParameter('coordinatefile');
% added by Yu-Le
if startsWith(f,'--')
    ext = 'LocMoFit';
else
    [p,fh,ext]=fileparts(f);
end
switch ext
    case {'.txt','.csv'}
        txt='on';
        tif='off';
    case {'.tif','.png','.jpg','.jpeg'}
        txt='off';
        tif='on';
    case '.mat'
        l=load([p f]);
        if isfield(l,'image')
            txt='off';
            tif='on';
        else
            txt='on';
            tif='off';            
        end
    case '.m'
        if isdeployed
            return
        end
                
        cf=pwd;
        cd(p)
        [~,fh]=fileparts(f);
        l=eval(fh);
        cd(cf);
        
        if isfield(l,'image')
            txt='off';
            tif='on';
        else
            txt='on';
            tif='off';            
        end
    case 'LocMoFit'
        modelType = obj.getPar('modelType');
        if ~isempty(modelType)
            switch modelType
                case 'Image'
                    txt='off';
                    tif='on';
                case 'Point'
                    txt='on';
                    tif='off';
                case ''
                    txt='off';
                    tif='on';
            end
        else
            txt='off';
            tif='on';
        end
end
obj.guihandles.labeling_efficiency.Visible=txt;
obj.guihandles.t_labelingefficiency.Visible=txt;
obj.guihandles.tif_density.Visible=tif;
obj.guihandles.tif_numbermode.Visible=tif;
obj.guihandles.tif_imagesizet.Visible=tif;
obj.guihandles.tif_imagesize.Visible=tif;
obj.guihandles.linkageerrort.Visible=txt;
obj.guihandles.linkageerrorfix.Visible=txt;
obj.guihandles.linkageerrorft.Visible=txt;
obj.guihandles.linkageerrorfree.Visible=txt;
end


function pard=guidef(obj)

% 
pard.coordinatefile.object=struct('String','plugins/+ROIManager/+Segment/hidden/MakeNPCCoordinates.m','Style','edit');
pard.coordinatefile.position=[1,1];
pard.coordinatefile.Width=3;
pard.coordinatefile.TooltipString=sprintf('.txt or .csv file with coordinates, .tif or .png file in which the pixel values are a \n measure for the concentration of the labels or a matlab function that \n returns the position of the lables in form of a structure with the fields .x .y .z');

pard.load_button.object=struct('String','Load','Style','pushbutton','Callback',{{@load_callback,obj}});
pard.load_button.position=[1,4];
pard.load_button.TooltipString=pard.coordinatefile.TooltipString;

% Added by Yu-Le
pard.useFitter_button.object = struct('Style','pushbutton','String', 'Use LocMoFIt', 'Callback', {{@useFitter_callback,obj}});
pard.useFitter_button.position = [2 4];

pard.setModPars_button.object = struct('Style','pushbutton','String', 'Set model pars', 'Callback', {{@setModPars_callback,obj}});
pard.setModPars_button.position = [2 3];


pard.tif_numbermode.object=struct('String',{{'Density (labels/um^2)','Number of labels'}},'Style','popupmenu');
pard.tif_numbermode.Width=1.5;
pard.tif_numbermode.position=[3,1];

pard.tif_density.object=struct('String','100','Style','edit');
pard.tif_density.position=[3,2.5];
pard.tif_density.Width=0.5;

pard.tif_imagesizet.object=struct('String','Image width (nm)','Style','text');
pard.tif_imagesizet.Width=1.5;
pard.tif_imagesizet.position=[3,3];

pard.tif_imagesize.object=struct('String','300','Style','edit');
pard.tif_imagesize.position=[3,4.5];
pard.tif_imagesize.Width=0.5;

pard.t_labelingefficiency.object=struct('String','Labeling efficiency','Style','text');
pard.t_labelingefficiency.position=[3,1];
pard.t_labelingefficiency.Width=1.5;
pard.t_labelingefficiency.TooltipString=sprintf('fraction of target molecules that carry a fluorophore');

pard.labeling_efficiency.object=struct('String','.5','Style','edit');
pard.labeling_efficiency.Width=.5;
pard.labeling_efficiency.position=[3,2.5];
pard.labeling_efficiency.TooltipString=pard.t_labelingefficiency.TooltipString;

pard.linkageerrort.object=struct('String','Linkage err fix:','Style','text');
pard.linkageerrort.position=[3,3];
pard.linkageerrort.Width=0.8;

pard.linkageerrorfix.object=struct('String','0','Style','edit');
pard.linkageerrorfix.Width=.4;
pard.linkageerrorfix.position=[3,3.8];

pard.linkageerrorft.object=struct('String','free:','Style','text');
pard.linkageerrorft.position=[3,4.2];
pard.linkageerrorft.Width=0.4;

pard.linkageerrorfree.object=struct('String','0','Style','edit');
pard.linkageerrorfree.Width=.4;
pard.linkageerrorfree.position=[3,4.6];


pard.modelt.object=struct('String','Model:','Style','text');
pard.modelt.Width=.5;
pard.modelt.position=[4,1];
pard.modelt.TooltipString=sprintf('Select model for blinking: \n simple: re-activations appear randomly in any frame \n PAFP: re-activations appear shortly after previous activation \n Dye: if a dye is not bleached, a re-activation occurs a random time point after the previous one (exponential distribution)');

pard.model.object=struct('String',{{'simple','PAFP','Dye'}},'Style','popupmenu');
pard.model.Width=.75;
pard.model.position=[4,1.5];
pard.model.TooltipString=pard.modelt.TooltipString;

pard.t2.object=struct('String','re-activations','Style','text');
pard.t2.position=[4,2.25];
pard.t2.Width=.75;

pard.blinks.object=struct('String','1','Style','edit');
pard.blinks.Width=.5;
pard.blinks.position=[4,3];
pard.blinks.TooltipString=sprintf('Number of re-activations. Zero means: only one actvation per fluorophore');
pard.t2.TooltipString=pard.blinks.TooltipString;

pard.EMon.object=struct('String','EM on','Style','checkbox');
pard.EMon.position=[4,4];
pard.EMon.Width=1;

pard.t3.object=struct('String','lifetime (fr)','Style','text');
pard.t3.position=[5,4];
pard.t3.Width=0.65;

pard.lifetime.object=struct('String','1','Style','edit');
pard.lifetime.Width=.35;
pard.lifetime.position=[5,4.65];
pard.lifetime.TooltipString=sprintf('average on-time (in frames) of an activated fluorophore');
pard.t3.TooltipString=pard.lifetime.TooltipString;

pard.t4.object=struct('String','mean photons','Style','text');
pard.t4.position=[5,1];
pard.t4.Width=.75;

pard.photons.object=struct('String','2000','Style','edit');
pard.photons.Width=.5;
pard.photons.position=[5,1.75];
pard.photons.TooltipString=sprintf('mean number of photons emitted by an activated fluorophore before off-switching \n (this number is distributed among the frames the fluorophore is on).');
pard.t4.TooltipString=pard.photons.TooltipString;

pard.photonsigmat.object=struct('String','+/-','Style','text');
pard.photonsigmat.position=[5,2.25];
pard.photonsigmat.Width=.25;
pard.photonsigma.object=struct('String','0','Style','edit');
pard.photonsigma.position=[5,2.45];
pard.photonsigma.Width=.35;

pard.t5.object=struct('String','BG/pixel','Style','text');
pard.t5.position=[5,3];
pard.t5.Width=0.5;

pard.background.object=struct('String','20','Style','edit');
pard.background.Width=.35;
pard.background.position=[5,3.5];
pard.background.TooltipString=sprintf('Background in photons/pixel/frame');
pard.t5.TooltipString=pard.background.TooltipString;

% PSF model
p(1).value=0;p(1).on={};p(1).off={'psf_file','load_button_psf'};
p(2).value=1;p(2).on=p(1).off; p(2).off=p(1).on;
pard.use_psf.object=struct('String','Experimental PSF:','Style','checkbox','Value',0,'Callback',{{@obj.switchvisible,p}});
pard.use_psf.position=[6,1];
pard.use_psf.Width=1.5;

pard.psf_file.object=struct('String','*_3dcal.mat','Style','edit','Visible','off');
pard.psf_file.position=[6,2.5];
pard.psf_file.Width=1.5;

pard.load_button_psf.object=struct('String','Load','Style','pushbutton','Callback',{{@loadpsf_callback,obj}},'Visible','off');
pard.load_button_psf.position=[6,4];



pard.t6.object=struct('String','Number of sites','Style','text');
pard.t6.position=[8,1];
pard.t6.Width=1.5;

pard.numberofsites.object=struct('String','6 5','Style','edit');
pard.numberofsites.Width=.5;
pard.numberofsites.position=[8,2.5];
pard.numberofsites.TooltipString=sprintf('Number of simulated structures. Pass on two numbers [M N] to generate a grid of MxN structures.');
pard.t6.TooltipString=pard.numberofsites.TooltipString;

pard.randomrot.object=struct('String','Random rot theta (deg):','Style','checkbox');
pard.randomrot.Width=1.75;
pard.randomrot.position=[7,1];

pard.randomrotangle.object=struct('String','15','Style','edit');
pard.randomrotangle.Width=.5;
pard.randomrotangle.position=[7,2.5];
pard.randomrot.TooltipString=sprintf('randomly rotate the structures by a maximum angle (degree) defined here. \n The structures are rotated additionally by a random angle around the z-axis.');
pard.randomrotangle.TooltipString=pard.randomrot.TooltipString;

pard.randomxy.object=struct('String','Random position (nm):','Style','checkbox');
pard.randomxy.Width=1.5;
pard.randomxy.position=[7,3];
pard.randomxy.TooltipString=sprintf('Structures are displaced in x,y,z by a random distance');

pard.randomxyd.object=struct('String','20','Style','edit');
pard.randomxyd.Width=.5;
pard.randomxyd.position=[7,4.5];
pard.randomxyd.TooltipString=pard.randomxy.TooltipString;

pard.savez.object=struct('String','save z','Style','checkbox','Value',1);
pard.savez.Width=0.5;
pard.savez.position=[2,2.2];
pard.savez.TooltipString=sprintf('Also simulate z-coordinate');
pard.savez.Optional=true;

pard.t7.object=struct('String','number of frames','Style','text');
pard.t7.position=[8,3];
pard.t7.Width=1.5;

pard.maxframes.object=struct('String','100000','Style','edit');
pard.maxframes.Width=.5;
pard.maxframes.position=[8,4.5];
pard.maxframes.TooltipString=sprintf('Maximum number of frames that contain all localizations');
pard.t7.TooltipString=pard.maxframes.TooltipString;

pard.savenow.object=struct('String',{{'No saving','Save ground truth','Save simulated locs'}},'Style','popupmenu');
pard.savenow.Width=1.2;
pard.savenow.position=[2,1];
pard.savenow.Optional=true;
pard.plugininfo.type='ROI_Analyze';

pard.plugininfo.description=sprintf(['SimulateSites is a localization based simulation engine for SMAP. It \n'...
    'uses as an input a list of localizations, a matlab function that\n'...
    'returns coordiantes or an image which defines a 2D structure. It \n'...
    'returns simulated localizations to SMAP using a realistic model for the \n'...
    'photophysics of the dye. Simulated structures are added to the RoiManager']);
end
