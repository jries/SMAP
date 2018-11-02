classdef Mirror<interfaces.DialogProcessor
    methods
        function obj=Mirror(varargin)    
            obj@interfaces.DialogProcessor(varargin{:}) ;  
        end
        
        function out=run(obj,p) 
            out=[];
            obj.setPar('undoModule','Mirror');
            notify(obj.P,'backup4undo');
            filenumber=p.dataselect.Value;
            file=obj.locData.files.file(filenumber);
            loc=obj.locData.loc;
            indfile=loc.filenumber==filenumber;
            pixelsize=file.info.cam_pixelsize_um;
%             roi=file.info.roi*pixelsize*1000;
            
            midpx=256*pixelsize(1)*1000;
            midpy=256*pixelsize(2)*1000;
%             midpy=roi(2)+roi(4)/2;
            
            switch p.mirrorpart.selection
                case 'all'
                    ind=true(size(loc.xnm));
                    midmx=midpx;
                    midmy=midpy;
                case 'top'
                    ind=loc.ynm<midpy;
                    midmx=midpx;
                    midmy=midpy*.5;
                case 'bottom'
                    ind=loc.ynm>=midpy;
                    midmx=midpx;
                    midmy=midpy*1.5;
                case 'left'
                    ind=loc.xnm<midpx;
                    midmx=midpx*.5;
                    midmy=midpy;
                case 'right'
                    ind=loc.xnm>midpx;
                    midmx=midpx*1.5;
                    midmy=midpy;
            end
            ind=ind&indfile;
            switch p.mirrormode.selection
                case 'top-bottom'
                    obj.locData.loc.ynm(ind)=2*midmy-loc.ynm(ind);
                case 'left-right'
                    obj.locData.loc.xnm(ind)=2*midmx-loc.xnm(ind);
                case 'both'
                    obj.locData.loc.ynm(ind)=2*midmy-loc.ynm(ind);
                    obj.locData.loc.xnm(ind)=2*midmx-loc.xnm(ind);
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

pard.dataselect.object.TooltipString='choose localization file data set';


pard.textb.object=struct('String','Part to mirror','Style','text');
pard.textb.position=[3,1];
pard.mirrorpart.object=struct('String','top|bottom|left|right|','Style','popupmenu','Value',2);
pard.mirrorpart.position=[3,2];


pard.textc.object=struct('String','mirror operation','Style','text');
pard.textc.position=[4,1];
pard.mirrormode.object=struct('String','top-bottom|left-right|both','Style','popupmenu','Value',1);
pard.mirrormode.position=[4,2];

pard.syncParameters={{'filelist_short','dataselect',{'String'}}};
pard.plugininfo.type='ProcessorPlugin';
end