classdef ApplyTransform<recgui.DialogProcessor
    methods
        function obj=ApplyTransform(varargin)        
                obj@recgui.DialogProcessor;
%                 if ~isempty(parameter)
%                     obj.setParameters(parameter)
%                 end
            if nargin>0
                obj.handle=varargin{1};
            
%             if nargin>1
%                 obj.makeGui(guidef);
            end
%             end   

        end
        
        function out=run(obj)
            p=obj.getGuiParameters.par;
            if strcmp(obj.resultshandle.Visible,'on')
                p.ploton=1;
            else
                p.ploton=0;
            end
%             locs=obj.locData.getlocRoi('frame','xnm','ynm',0);
%             p.resultstabs(1)=obj.guihandles.resultstab1;
%             p.resultstabs(2)=obj.guihandles.resultstab2;
%             p.maxframeall=max(obj.locData.loc.frame);
%             drift=finddriftfeature(locs,p);
%             locsall=copyfields([],obj.locData.loc,{'xnm','ynm','frame'});
%             grouplocsall=copyfields([],obj.locData.grouploc,{'xnm','ynm','frame'});
%             locsnew=applydriftcorrection(drift,locsall);
%             grouplocsnew=applydriftcorrection(drift,grouplocsall);
%             obj.locData.loc=copyfields(obj.locData.loc,locsnew,{'xnm','ynm'});
%             obj.locData.grouploc=copyfields(obj.locData.grouploc,grouplocsnew,{'xnm','ynm'});
%             
%             fn=obj.locData.files(obj.locData.loc.filenumber(1)).file.name;
%             fnn=[fn(1:end-7) '_driftc_sml.mat'];
%             obj.locData.savelocs(fnn);
%             
%             out=1;
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function attachLocData(obj,locData)
            attachLocData@recgui.GuiProcessor(obj,locData);
%             addlistener(obj.locData,'synchronizeGui',@obj.synchronizeGui);
            addlistener(obj.locData,'loaded',@obj.updateGui);
        end
        function updateGui(obj,event,data)
            ff=obj.locData.files.file;
            str={};
            for k=1:length(ff)
                str{end+1}=ff(k).name;
            end
            obj.guihandles.targetselect.String=str;
            obj.guihandles.refselect.String=str;
        end

    end
    methods(Static)
        function info=info(obj)
            info.name='ApplyTransform';
            info.class=@ApplyTransform;
            info.tag='ApplyTransform';
%             obj.info=info;
        end

    end
end




function pard=guidef
pard.texta.object=struct('String','target','Style','text');
pard.texta.position=[1,1];
pard.textb.object=struct('String','reference','Style','text');
pard.textb.position=[1,2];
pard.targetselect.object=struct('Style','popupmenu','String','File');
pard.targetselect.position=[2,1];
pard.refselect.object=struct('Style','popupmenu','String','layer1|layer2|layer3|layer4|layer5');
pard.refselect.position=[2,2];

pard.targetpart.object=struct('Style','popupmenu','String','all|left|right|top|bottom');
pard.targetpart.position=[3,1];
pard.refpart.object=struct('Style','popupmenu','String','all|left|right|top|bottom');
pard.refpart.position=[3,2];

pard.targetmirror.object=struct('Style','popupmenu','String','no mirror|left-right|up-down');
pard.targetmirror.position=[4,1];
pard.refmirror.object=struct('Style','popupmenu','String','no mirror|left-right|up-down');
pard.refmirror.position=[4,2];

pard.savefile.object=struct('Style','edit','String','settings/temp/temp_T.mat');
pard.savefile.position=[8,1];
pard.savefile.Width=3;

pard.savebutton.object=struct('Style','pushbutton','String','save');
pard.savebutton.position=[8,4];

pard.registeronlocs.object=struct('Style','checkbox','String','register locs');
pard.registeronlocs.position=[1,3.5];

pard.transform.object=struct('Style','popupmenu','String','translate|affine|LWM');
pard.transform.position=[2,3.5];
end