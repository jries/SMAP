classdef storerois<interfaces.DialogProcessor
    properties
        rois
    end
    methods
        function obj=storerois(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'sr_layerson','linewidth_roi','znm_min','znm_max','sr_pixrec','layernames'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)           

        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function pushbutton(obj,a,b,action)
            switch a.String
                case 'add'
                    formathandle=obj.getPar('guiFormat');
                    roih=obj.getPar('sr_roihandle');
                    roisave.roimode=class(roih);
                    roisave.position=roih.getPosition;
                    linewidth=formathandle.getSingleGuiParameter('linewidth_roi');
                    roisave.linewidth=linewidth;
                    obj.rois{end+1}=roisave;
                    setroilist(obj)
                case 'remove'
                    line=obj.guihandles.roilist.Value;
                    obj.rois(line)=[];
                    setroilist(obj);
                case 'apply'
                    formathandle=obj.getPar('guiFormat');
                    line=obj.guihandles.roilist.Value;
                    roih=obj.rois{line};
                    roih.isvalid=true;
                    if strcmp(roih.roimode,'imline')
                        formathandle.setGuiParameters(struct('linewidth_roi',roih.linewidth))
                    end
                    formathandle.roiset(roih); 
                case 'save'
                    fn=obj.getPar('lastSMLFile');
                    fn=strrep(fn,'_sml.mat','_rois.mat');
                    [f,path]=uiputfile(fn);
                    rois=obj.rois;
                    if f
                        save([path f],'rois')
                    end
                case 'load'
                    path=fileparts(obj.getPar('lastSMLFile'));
                    [f,path]=uigetfile([path filesep '.mat']);
                    if f
                        l=load([path filesep f]);
                        obj.rois=l.rois;
                        setroilist(obj);
                    end
                case 'overview'
                    f=obj.getPar('sr_figurehandle');
                    f2=f.copy;hold on
                    for k=1:length(obj.rois)
                        roih=obj.rois{k};
                        plot(roih.position(:,1),roih.position(:,2),'m','LineWidth',3)
                        text(mean(roih.position(:,1)),mean(roih.position(:,2)),num2str(k),"BackgroundColor",[1,1,1])                
                    end
                    %  for k=1:length(obj.rois)
                    %     text(mean(roih.position(:,1)),mean(roih.position(:,2)),num2str(k))                
                    % end                  

            end
        end
    end
end

function setroilist(obj)
ff='%2.1f';
for k=1:length(obj.rois)
    pos=obj.rois{k}.position;
    names{k}=k+ ". "+(obj.rois{k}.roimode)+ " " + num2str(pos(1),ff)+ "," + num2str(pos(1),ff);
end
obj.guihandles.roilist.String=names;
obj.guihandles.roilist.Value=k;
end




function pard=guidef(obj)
buttonwidth=0.5;
pard.add.object=struct('String','add','Style','pushbutton','Callback',{{@obj.pushbutton,'add'}});
pard.add.position=[1,1];
pard.add.Width=buttonwidth;
pard.apply.object=struct('String','apply','Style','pushbutton','Callback',{{@obj.pushbutton,'apply'}});
pard.apply.position=[2,1];
pard.apply.Width=buttonwidth;
pard.remove.object=struct('String','remove','Style','pushbutton','Callback',{{@obj.pushbutton,'remove'}});
pard.remove.position=[3,1];
pard.remove.Width=buttonwidth;
pard.save.object=struct('String','save','Style','pushbutton','Callback',{{@obj.pushbutton,'save'}});
pard.save.position=[4,1];
pard.save.Width=buttonwidth;
pard.load.object=struct('String','load','Style','pushbutton','Callback',{{@obj.pushbutton,'load'}});
pard.load.position=[5,1];
pard.load.Width=buttonwidth;
pard.ov.object=struct('String','overview','Style','pushbutton','Callback',{{@obj.pushbutton,'load'}});
pard.ov.position=[6,1];
pard.ov.Width=buttonwidth;

pard.roilist.object=struct('String','','Style','listbox');
pard.roilist.position=[6,buttonwidth+1];
pard.roilist.Width=5-buttonwidth;
pard.roilist.Height=6;

pard.plugininfo.description=sprintf('stores and retrieves ROIs');
pard.plugininfo.type='ProcessorPlugin';
end