classdef TransformLocsManually<interfaces.DialogProcessor
%     Manually change position, magnification and rotation of localization
%     data
    methods
        function obj=TransformLocsManually(varargin)    
            obj@interfaces.DialogProcessor(varargin{:}) ;  
        end
        
        function out=run(obj,p) 
            out=[];
            obj.setPar('undoModule','Mirror');
            notify(obj.P,'backup4undo');
            filenumber=p.dataselect.Value;
            file=obj.locData.files.file(filenumber);
            loc=obj.locData.loc;
            gloc=obj.locData.grouploc;
            indfile=loc.filenumber==filenumber;
            indfileg=gloc.filenumber==filenumber;
            
            loc.xnm(indfile)=loc.xnm(indfile)+p.shiftx;
            gloc.xnm(indfileg)=gloc.xnm(indfileg)+p.shiftx;
            loc.ynm(indfile)=loc.ynm(indfile)+p.shifty;
            gloc.ynm(indfileg)=gloc.ynm(indfileg)+p.shifty;           
            if isfield(loc,'znm')
                loc.znm(indfile)=loc.znm(indfile)+p.shiftz;
                gloc.znm(indfileg)=gloc.znm(indfileg)+p.shiftz;
            end
            
            if p.rotcset
                cx=p.rotc(1);cy=p.rotc(2);
            else
                cx=p.sr_pos(1);cy=p.sr_pos(2);
            end
            
            xh=(loc.xnm(indfile)-cx)*p.mag(1);
            yh=(loc.ynm(indfile)-cy)*p.mag(end);
            xhg=(gloc.xnm(indfileg)-cx)*p.mag(1);
            yhg=(gloc.ynm(indfileg)-cy)*p.mag(end);
            loc.xnm(indfile)=cosd(p.rot)*(xh)+sind(p.rot)*yh+cx;
            loc.ynm(indfile)=-sind(p.rot)*(xh)+cosd(p.rot)*yh+cy;
            gloc.xnm(indfileg)=cosd(p.rot)*(xhg)+sind(p.rot)*yhg+cx;
            gloc.ynm(indfileg)=-sind(p.rot)*xhg+cosd(p.rot)*yhg+cy;
            
            
            obj.locData.loc.xnm=loc.xnm;
            obj.locData.loc.ynm=loc.ynm;
            
            obj.locData.grouploc.xnm=gloc.xnm;
            obj.locData.grouploc.ynm=gloc.ynm;
            if isfield(obj.locData.loc,'znm')
                obj.locData.grouploc.znm=gloc.znm;
                obj.locData.loc.znm=loc.znm;
            end
%             
        end
        function pard=guidef(obj)
            pard=guidef;
        end

        function updateGui(obj,event,data)
            ff=obj.locData.files.file;
            str={};
            for k=1:length(ff)
                if 0% ff(k).istiff
                    str{end+1}=' ';
                else
                    str{end+1}=ff(k).name;
                end
            end
            obj.guihandles.dataselect.Value=min(obj.guihandles.dataselect.Value,length(ff));
            obj.guihandles.dataselect.String=str; 
        end
    end
end



function pard=guidef
pard.texta.object=struct('String','dataset','Style','text');
pard.texta.position=[1,1];

pard.dataselect.object=struct('Style','popupmenu','String','File');
pard.dataselect.position=[2,1];
pard.dataselect.Width=3;
pard.dataselect.object.TooltipString='choose localization file data set';


pard.textb.object=struct('String','shift x','Style','text');
pard.textb.position=[3,1];
pard.textb.Width=0.5;
pard.shiftx.object=struct('String','0','Style','edit');
pard.shiftx.position=[3,1.5];
pard.shiftx.Width=0.5;
pard.shiftyt.object=struct('String','shift y','Style','text');
pard.shiftyt.position=[3,2];
pard.shiftyt.Width=0.5;
pard.shifty.object=struct('String','0','Style','edit');
pard.shifty.position=[3,2.5];
pard.shifty.Width=0.5;

pard.shiftzt.object=struct('String','shift z','Style','text');
pard.shiftzt.position=[3,3];
pard.shiftzt.Width=0.5;
pard.shiftz.object=struct('String','0','Style','edit');
pard.shiftz.position=[3,3.5];
pard.shiftz.Width=0.5;

pard.rott.object=struct('String','rotate (deg)','Style','text');
pard.rott.position=[4,1];
pard.rott.Width=1;
pard.rot.object=struct('String','0','Style','edit');
pard.rot.position=[4,2];
pard.rot.Width=0.5;

pard.magt.object=struct('String','magnification','Style','text');
pard.magt.position=[4,3];
pard.magt.Width=1;
pard.mag.object=struct('String','1','Style','edit');
pard.mag.position=[4,4];
pard.mag.Width=0.5;


pard.rotcset.object=struct('String','rotation center x,y (otherwise FoV)','Style','checkbox');
pard.rotcset.position=[5,1];
pard.rotcset.Width=2;
pard.rotc.object=struct('String','12500 12500','Style','edit');
pard.rotc.position=[5,3];
pard.rotc.Width=1;

pard.syncParameters={{'filelist_short','dataselect',{'String'}}};
pard.plugininfo.type='ProcessorPlugin';
pard.inputParameters={'sr_pos','sr_size'};
pard.plugininfo.description='Manually change position, magnification and rotation of localization data.';
end