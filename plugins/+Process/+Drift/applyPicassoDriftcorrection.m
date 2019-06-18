classdef applyPicassoDriftcorrection<interfaces.DialogProcessor
%     Applies drift correction file from Picasso to current localization
%     data
    methods
        function obj=applyPicassoDriftcorrection(varargin)        
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.history=true;
                obj.showresults=true;
                obj.guiselector.show=true;
        end
        
        function out=run(obj,p)
            out=[];
            obj.setPar('undoModule','applydriftcorrection');
            notify(obj.P,'backup4undo');
            drift=table2array(readtable(p.file));
            ax=initaxis(p.resultstabgroup,'drift');
            frame=1:size(drift,1);
            plot(ax,frame,drift(:,1),frame,drift(:,2));
            if size(drift,2)>2
                ax=initaxis(p.resultstabgroup,'drift z');
                plot(ax,frame,drift(:,3));
            end
            
            pixelsize=obj.getPar('cam_pixelsize_um')*1000;

            if p.correctxy
                dr.xy.mirror='none';
                
                dr.xy.x=drift(:,1)*pixelsize(1);
                dr.xy.y=drift(:,2)*pixelsize(end);
            end
            if p.correctz&& size(drift,2)>2
%                 dr.z=driftinfo.z;
                dr.z.z=drift(:,3);
            end
            locsnew=applydriftcorrection(dr,obj.locData.loc);
            obj.locData.loc.xnm=locsnew.xnm;
            obj.locData.loc.ynm=locsnew.ynm;

            obj.locData.files.file(1).driftinfo=dr;
            fn=obj.locData.files.file(1).name;
            if strfind(fn,'_sml')
                fnn=strrep(fn,'_sml','_adriftc_sml');
            else
                fnn=strrep(fn,'fitpos','adriftc_sml');
            end
            if p.save_dc
                obj.locData.savelocs(fnn); 
            end
            obj.locData.regroup;
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function loadfile(a,b,obj)
oldf=obj.getSingleGuiParameter('file');
if ~exist(oldf,'file')
    oldf=obj.getPar('lastSMLFile');
end
if isempty(oldf)
    directory='';
else
    directory=fileparts(oldf);
end
[f,path]=uigetfile([directory filesep '*.txt']);
if f
    obj.setGuiParameters(struct('file',[path f]));
end
end


function pard=guidef(obj)

pard.correctxy.object=struct('String','Correct xy-drift','Style','checkbox','Value',1);
pard.correctxy.position=[1,1];

pard.correctz.object=struct('String','Correct z-drift','Style','checkbox','Value',0);
pard.correctz.position=[2,1];


pard.file.object=struct('String','','Style','edit');
pard.file.position=[3,1];
pard.file.Width=3.5;

pard.loadfile.object=struct('String','load','Style','pushbutton','Callback',{{@loadfile,obj}});
pard.loadfile.position=[3,4.5];
pard.loadfile.Width=0.5;

pard.save_dc.object=struct('String','Save driftcorrected SML','Style','checkbox','Value',0);
pard.save_dc.position=[4,1];
pard.save_dc.Width=2;
pard.save_dc.Optional=true;

pard.plugininfo.name='apply drift correction from Picasso';
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='Applies drift correction file from Picasso to current localization data';

end