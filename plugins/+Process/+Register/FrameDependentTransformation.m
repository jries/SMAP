classdef FrameDependentTransformation<interfaces.DialogProcessor
%     calculates transformation (global or local) based on localizations
%     (e.g. multi-color beads or fluorophores in ratiometric imaging)
    properties
        isz=0;
        transformation=[];
    end
    methods
        function obj=FrameDependentTransformation(varargin)        
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
            pard=guidef(obj);
        end
        function initGui(obj)
%             obj.addSynchronization('transformationfile',obj.guihandles.Tfile,'String');
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
%         function browse_callback(obj,a,b)
%             fn=obj.guihandles.Tfile.String;
%             [f,path]=uigetfile(fn,'Open transformation file _T.mat');
%             if f
%                 Tload=load([path f]);
%                 if ~isfield(Tload,'transformation')
%                     msgbox('could not find transformation in file. Load other file?')
%                 end
%                 obj.guihandles.Tfile.String=[path f];
%                 obj.guihandles.useT.Value=1;
%                 obj.setPar('transformationfile',[path f]);
%             end
%         end
    end
end



function pard=guidef(obj)
% pard.Tfile.object=struct('Style','edit','String','');
% pard.Tfile.position=[8,1];
% pard.Tfile.Width=3;
% pard.Tfile.object.TooltipString='default file for transformation matrix. You can select new file after transformation has been calculated.';
% 
% pard.browse.object=struct('Style','pushbutton','String','load T','Callback',@obj.browse_callback);
% pard.browse.position=[8,4];
% pard.browse.object.TooltipString='Save the newly calculated transformation matrix.';


pard.texttt.object=struct('String','Transformation:','Style','text');
pard.texttt.position=[2,3];
pard.transform.object=struct('Style','popupmenu','String','projective|affine|similarity|polynomial|lwm|pwl');
pard.transform.position=[3,3];
pard.transform.object.TooltipString='select one of Matlabs transformations. Not all might work.';

pard.transformparam.object=struct('Style','edit','String','3');
pard.transformparam.position=[4,3];
pard.transformparam.object.TooltipString='Parameter for lwm and polynomial';


pard.updatesmlfile.object=struct('Style','checkbox','String','write T to .sml','Value',1);
pard.updatesmlfile.position=[6,3];
pard.updatesmlfile.object.TooltipString='If checked, the transformation file is appended to the .sml file and saved there as well when you click save T';


pard.inputParameters={'currentfileinfo'};
pard.plugininfo.description='calculates transformation (global or local) based on localizations (e.g. multi-color beads or fluorophores in ratiometric imaging)';
pard.plugininfo.type='ProcessorPlugin';
end