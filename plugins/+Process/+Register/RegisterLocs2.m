classdef RegisterLocs2<interfaces.DialogProcessor
%     calculates transformation (global or local) based on localizations
%     (e.g. multi-color beads or fluorophores in ratiometric imaging)
    properties
        isz=0;
        transformation=[];
        register_parameters;
    end
    methods
        function obj=RegisterLocs2(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.showresults=true;
            obj.history=true;
        end
        
        function out=run(obj,p)
            if p.useT
                transformation=loadtransformation(obj,p.Tfile,p.dataselect.Value);
                if isempty(transformation)
                    out.error='selected transformation file does not have a valid transformation. Uncheck use inital T';
                    return 
                end
            end
            p.isz=obj.isz;
            p.register_parameters=obj.register_parameters;
            fn=fieldnames(p.register_parameters);
            nmax=1;
            for k=1:length(fn)
                nmax=max(nmax,length(p.register_parameters.(fn{k})));
            end
            
%             rp=p.register_parameters;
            for n=1:nmax
                for k=1:length(fn)
                    ind=min(n,length(p.register_parameters.(fn{k})));
                    rp(n).(fn{k})=p.register_parameters.(fn{k})(ind);
                end
            end
            
%             if p.saveoldformat
%                 obj.transformation=transform_locs(obj.locData,p);
%             else
            for n=1:nmax
                p.register_parameters=rp(n);
                 p.repetition=num2str(n);
                obj.transformation=transform_locsN(obj.locData,p);
                p.Tfile=obj.transformation;
                p.useT=true;
            end
%             end
            fv=p.dataselect.Value;
            obj.locData.files.file(fv).transformation=obj.transformation;
            if obj.processorgui==false %run from WF
                f=obj.locData.files.file(fv).name;
                fn=strrep(f,'_sml.mat','_T.mat');
                   obj.guihandles.Tfile.String=fn;
                   transformation=obj.transformation; 
                   save(fn,'transformation');
                   obj.setPar('transformationfile',fn);
            end
            out=[];
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function makeGui(obj)
            makeGui@interfaces.DialogProcessor(obj);
            obj.guihandles.save.Callback=@obj.save_callback;
            obj.guihandles.browse.Callback=@obj.browse_callback;
            obj.guihandles.parameters.Callback=@obj.parameters_callback;
            obj.guihandles.useT.Callback=@obj.useT_callback; %gray out unused fields
            obj.parameters_callback(0,0)
        end
        function initGui(obj)
            obj.addSynchronization('transformationfile',obj.guihandles.Tfile,'String');
        end
        function save_callback(obj,a,b)
            if isempty(obj.transformation)
                errordlg('first calculate a transformation')
            else
                 f=obj.locData.files.file(1).name;
                fn=strrep(f,'_sml.mat','_T.mat');
%                 fn=obj.guihandles.Tfile.String;
                [f,path]=uiputfile(fn,'Save last transformation as transformation file _T.mat');
                if f
                    obj.guihandles.Tfile.String=[path f];
                    transformation=obj.transformation; 
                    save([path f],'transformation');
                    obj.setPar('transformationfile',[path f]);
                end 
                if obj.getSingleGuiParameter('updatesmlfile')
                    filenumber=obj.getSingleGuiParameter('dataselect').Value;
                    obj.locData.files.file(filenumber).transformation=obj.transformation;
                    obj.locData.savelocs(obj.locData.files.file(filenumber).name,[],[],[],[],filenumber);
                end
            end
        end
        function browse_callback(obj,a,b)
            fn=obj.guihandles.Tfile.String;
            [f,path]=uigetfile(fn,'Open transformation file _T.mat');
            if f
                Tload=load([path f]);
                if ~isfield(Tload,'transformation')
                    msgbox('could not find transformation in file. Load other file?')
                end
                obj.guihandles.Tfile.String=[path f];
                obj.guihandles.useT.Value=1;
                obj.setPar('transformationfile',[path f]);
            end
        end
        function parameters_callback(obj,callobj,b)
            if isempty(obj.register_parameters)        
                par.pixelsizenm=[500 50];
                par.maxshift_corr=15000;
                par.maxlocsused=[500000 1e9];
                par.maxshift_match=[1250 150];
                par.initial_mag=1;
                par.initialshiftx=0;
                 par.initialshifty=0;
                 obj.register_parameters=par;
                 return
            else 
                par=obj.register_parameters;
            end
            
            fn=fieldnames(par);
            for k=1:length(fn)
                if length(par.(fn{k}))>1 %vector
                    par.(fn{k})=num2str(par.(fn{k}));
                end
            end
            if isa(callobj,'matlab.ui.control.UIControl')

                [settings, button] = settingsdlg(...
                    'Description', 'Parameters for registration',... 
                    'title'      , 'Parameters',...  
                    {'Pixelsize (nm) for correlation';'pixelsizenm'}, par.pixelsizenm,...
                    {'Max shift for correlation (nm)';'maxshift_corr'}, par.maxshift_corr,...
                    {'Max locs for matching';'maxlocsused'}, par.maxlocsused,...
                    {'Initial magnification';'initial_mag'}, par.initial_mag,...
                    {'Initial shift X (nm)';'initialshiftx'}, par.initialshiftx,...
                    {'Initial shift Y (nm)';'initialshifty'}, par.initialshifty,...
                    {'Max shift matching (nm)';'maxshift_match'}, par.maxshift_match);

                if strcmpi(button,'ok')
                    par=addFields(par,settings);
                    for k=1:length(fn)
                        if ischar(par.(fn{k})) %vector
                            par.(fn{k})=str2num(par.(fn{k}));
                        end
                    end
                    
                    obj.register_parameters=par;
                end
            else
                obj.register_parameters=par;
            end
        end
        function useT_callback(obj,a,b)
        end
    end
end



function pard=guidef
pard.texta.object=struct('String','target:','Style','text');
pard.texta.position=[2,2.05];
% pard.textb.object=struct('String','reference','Style','text');
% pard.textb.position=[1,1];
pard.dataselect.object=struct('Style','popupmenu','String','File');
pard.dataselect.position=[2,1];
pard.dataselect.object.TooltipString='file for target localizations';

pard.allfiles.object=struct('String','use all files','Style','checkbox');
pard.allfiles.position=[3,1];

% 
% pard.refselect.object=struct('Style','popupmenu','String','File');
% pard.refselect.position=[2,1];
% pard.refselect.object.TooltipString='file for reference localizations';

% pard.refpos.object=struct('Style','popupmenu','String','center|up|down|left|right');
% pard.refpos.position=[3,1];
% pard.refpos.object.TooltipString='center position of reference';

pard.targetpos.object=struct('Style','popupmenu','String','center|top|bottom|left|right');
pard.targetpos.position=[3,2];
pard.targetpos.object.TooltipString='center position of target';

pard.targetmirror.object=struct('Style','popupmenu','String','no mirror|left-right|up-down|both');
pard.targetmirror.position=[4,2];
pard.targetmirror.object.TooltipString='mirror target part';

pard.uselayers.object=struct('Style','checkbox','String','use layers');
pard.uselayers.position=[4,1];
pard.uselayers.object.TooltipString='if checked, use only localizations in ROI and selected layer, make sure to select both halves';


pard.texttl.object=struct('String','T:','Style','text');
pard.texttl.position=[5,1];
pard.texttl.Width=0.15;

pard.targetlayer.object=struct('Style','popupmenu','String','layer1|layer2|layer3|layer4|','Value',2);
pard.targetlayer.position=[5,1.15];
pard.targetlayer.object.TooltipString='layer';
pard.targetlayer.Width=0.85;

pard.texttr.object=struct('String','R:','Style','text');
pard.texttr.position=[6,1];
pard.texttr.Width=0.15;

pard.reflayer.object=struct('Style','popupmenu','String','layer1|layer2|layer3|layer4|');
pard.reflayer.position=[6,1.15];
pard.reflayer.object.TooltipString='layer';
pard.reflayer.Width=0.85;

pard.Tfile.object=struct('Style','edit','String','settings/temp/temp_T.mat');
pard.Tfile.position=[8,1];
pard.Tfile.Width=3;
pard.Tfile.object.TooltipString='default file for transformation matrix. You can select new file after transformation has been calculated.';

pard.browse.object=struct('Style','pushbutton','String','load T');
pard.browse.position=[8,4];
pard.browse.object.TooltipString='Save the newly calculated transformation matrix.';

pard.parameters.object=struct('Style','pushbutton','String','Parameters');
pard.parameters.position=[2.5,4];
pard.parameters.object.TooltipString='additional parameters';
pard.parameters.Height=1.5;

pard.texttt.object=struct('String','Transformation:','Style','text');
pard.texttt.position=[2,3];
pard.transform.object=struct('Style','popupmenu','String','projective|affine|similarity|polynomial|lwm|pwl');
pard.transform.position=[3,3];
pard.transform.object.TooltipString='select one of Matlabs transformations. Not all might work.';

pard.transformparam.object=struct('Style','edit','String','3');
pard.transformparam.position=[4,3];
pard.transformparam.object.TooltipString='Parameter for lwm and polynomial';

pard.useT.object=struct('Style','checkbox','String','use inital T');
pard.useT.position=[4,4];
pard.useT.object.TooltipString='';
pard.Width=2;

% pard.saveoldformat.object=struct('Style','checkbox','String','old  format');
% pard.saveoldformat.position=[6,3];
% pard.saveoldformat.object.TooltipString='';


pard.save.object=struct('Style','pushbutton','String','save T');
pard.save.position=[6.5,4];
pard.save.object.TooltipString='';
pard.save.Height=1.5;

pard.updatesmlfile.object=struct('Style','checkbox','String','add T to current .sml','Value',1);
pard.updatesmlfile.position=[6,2.5];
pard.updatesmlfile.object.TooltipString='If checked, the transformation file is appended to the .sml file and saved there as well when you click save T';
pard.updatesmlfile.Height=1.5;

pard.syncParameters={{'filelist_short','dataselect',{'String'}}};
pard.inputParameters={'currentfileinfo'};
pard.plugininfo.description='calculates transformation (global or local) based on localizations (e.g. multi-color beads or fluorophores in ratiometric imaging)';
pard.plugininfo.name='Register Localizations';
pard.plugininfo.type='ProcessorPlugin';
end