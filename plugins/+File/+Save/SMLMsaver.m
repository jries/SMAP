classdef SMLMsaver<interfaces.DialogProcessor
%     saves localizations in SMAP proprietary _sml.mat MATLAB format
    properties
        excludesavefields={'groupindex','numberInGroup','colorfield'};
    end
    methods
        function obj=SMLMsaver(varargin)  
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'filelist_long','mainfile','mainGui','numberOfLayers','sr_layerson'};
        end
        
        function out=save(obj,p)
            obj.status('save sml file')
            lastfile=obj.getPar('lastSMLFile');
            
            if isempty(lastfile)
                fn=p.filelist_long.selection;
                ind=strfind(fn,'_sml');
                if isempty(ind)
                    ind=strfind(fn,'_fitpos');
                end
                if isempty(ind)
                    ind=length(fn)-3;
                end
                of=[fn(1:ind-1) '_sml.mat'];
            else
                of=lastfile;
            end
              
            
            [f,path]=uiputfile(of);
            if f
                if isempty(strfind(f,'_sml'))
                    f(end-3:end)=[];
                    f=[f '_sml.mat'];
                end   
            par=obj.getAllParameters;
            
            par.saveroi=par.savevisible;
            savelocData=obj.locData.copy;
%             savelocData.loc=rmfield(savelocData.loc,obj.excludesavefields);
            savesml(savelocData,[path f],par,obj.excludesavefields)
            
            end
            obj.status('save done')
          
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function run(obj,p)
            obj.save(p)
        end
    end
end

function outputfields_callback(a,b,obj)
fn=fieldnames(obj.locData.loc);
tosave=setdiff(fn,obj.excludesavefields);
[save,inds]=checknames(fn,tosave);
if ~isempty(save)
    obj.excludesavefields=fn(~inds);
end
end


function pard=guidef(obj)
pard.savevisible.object=struct('Style','checkbox','Visible','on','String','only save visible','Value',0);
pard.savevisible.position=[1,1];
pard.savevisible.Width=4;
pard.savevisible.object.TooltipString='save only those filtered localizations that have been used to render the image';

pard.savefile.object=struct('Style','checkbox','Visible','on','String','Only: ','Value',0);
pard.savefile.position=[2,1];
pard.savefile.Width=1;
pard.savefile.object.TooltipString='save only selected file';

pard.dataselect.object=struct('Style','popupmenu','Visible','on','String',{{'empty'}});
pard.dataselect.position=[2,2];
pard.dataselect.Width=3;
pard.dataselect.object.TooltipString='save only selected file';

pard.selectfields.object=struct('Style','pushbutton','String','Fields to save','Callback',{{@outputfields_callback,obj}});
pard.selectfields.object.TooltipString='Select which fields to save. Use preview before.';
pard.selectfields.position=[1,3];
pard.selectfields.Width=2;
pard.selectfields.Optional=true;
            
pard.syncParameters={{'filelist_short','dataselect',{'String'}}};
pard.plugininfo.type='SaverPlugin';
pard.plugininfo.description='saves localizations in SMAP proprietary _sml.mat MATLAB format';

end