classdef Get2CIntLoc<interfaces.DialogProcessor
        % Get2CIntLoc gets intensities corresponding localizations in two channels.
    methods
        function obj=Get2CIntLoc(varargin)        
                obj@interfaces.DialogProcessor(varargin{:}) ;   
            obj.history=true;
            obj.showresults=false;
        end
        
        function out=run(obj,p)
           out=[];
           transformation=loadtransformation(obj,p.Tfile,p.dataselect.Value);
%             transformation=loadtransformation(p.Tfile,'transformation',p.dataselect.Value);
            if isempty(transformation)
                out.error='selected transformation file does not have a valid transformation';
                return
            end
            obj.locData.files.file(p.dataselect.Value).transformation=transformation;
            file=obj.locData.files.file(p.dataselect.Value);
            loc=get2Clocintensities(obj.locData.loc,transformation,file,p);
            obj.locData.loc=copyfields(obj.locData.loc,loc);
            obj.locData.regroup;
             obj.setPar('locFields',fieldnames(obj.locData.loc))
        end
        
        function pard=guidef(obj)
            pard=guidef;
        end

        function initGui(obj)
            obj.addSynchronization('filelist_short',obj.guihandles.dataselect,'String')
            obj.addSynchronization('transformationfile',obj.guihandles.Tfile,'String');
        end

        function makeGui(obj)
            makeGui@interfaces.DialogProcessor(obj);
            obj.guihandles.loadbutton.Callback=@obj.loadbutton;
        end
        function loadbutton(obj,a,b)
            fn=obj.guihandles.Tfile.String;
            [f,path]=uigetfile(fn,'Select transformation file _T.mat');
            if f
                Tload=load([path f]);
                if ~isfield(Tload,'transformation')
                    msgbox('could not find transformation in file. Load other file?')
                end
                obj.guihandles.Tfile.String=[path f];
            end      
        end        
    end
end




function pard=guidef
pard.texta.object=struct('String','dataset','Style','text');
pard.texta.position=[1,1];

pard.dataselect.object=struct('Style','popupmenu','String','File');
pard.dataselect.position=[2,1];

% pard.datapart.object=struct('Style','popupmenu','String','all|left|right|top|bottom');
% pard.datapart.position=[3,1];


pard.Tfile.object=struct('Style','edit','String','*_T.mat');
pard.Tfile.position=[8,1];
pard.Tfile.Width=3;

pard.loadbutton.object=struct('Style','pushbutton','String','load');
pard.loadbutton.position=[8,4];
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='Get2CIntLoc gets intensities corresponding localizations in two channels.';
pard.plugininfo.name='2C intensities from Localizations';
end