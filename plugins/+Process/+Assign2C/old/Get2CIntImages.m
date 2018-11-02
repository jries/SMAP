classdef Get2CIntImages<interfaces.DialogProcessor
    methods
        function obj=Get2CIntImages(varargin)   
            obj@interfaces.DialogProcessor(varargin) ;
        end
        
        function out=run(obj,p)
            load(p.Tfile)
            file=obj.locData.files.file(p.dataselect.Value);
            fo=strrep(file.name,'_sml.mat','_dc_sml.mat');
%             fo=[file.name(1:end-7) 'dc_sml.mat'];
            loc=get2CIntfromImages(obj.locData.loc,transformation,file,p);
            obj.locData.loc=copyfields(obj.locData.loc,loc);
            obj.locData.savelocs(fo);
            obj.locData.regroup;
            obj.setPar('locFields',fieldnames(obj.locData.loc))
%             notify(obj.locData,'updateGui')
        end
        function pard=guidef(obj)
            pard=guidef;
        end
%         function attachLocData(obj,locData)
%             attachLocData@recgui.GuiProcessor(obj,locData);
% %             addlistener(obj.locData,'synchronizeGui',@obj.synchronizeGui);
%             addlistener(obj.locData,'loaded',@obj.updateGui);
%         end
        function makeGui(obj)
            makeGui@interfaces.DialogProcessor(obj);
            obj.guihandles.loadbutton.Callback=@obj.loadbutton;
        end
        function updateGui(obj,event,data)
            ff=obj.locData.files.file;
            str={};
            for k=1:length(ff)
                    str{end+1}=ff(k).name;
            end
            obj.guihandles.dataselect.String=str; 
            
        end
        function loadbutton(obj,a,b)
            fn=obj.guihandles.Tfile.String;
            [f,path]=uigetfile(fn,'Select transformation file _T.mat');
            if f
                obj.guihandles.Tfile.String=[path f];
            end      
        end

    end
    methods(Static)
        function info=info(obj)
            info.name='Get2CIntImages';
            info.class=@Get2CIntImages;
            info.tag='Get2CIntImages';
        end

    end
end




function pard=guidef

pard.dataselect.object=struct('Style','popupmenu','String','File');
pard.dataselect.position=[1,1];

%sum
pard.checksum.object=struct('Style','checkbox','String','sum','Value',1);
pard.checksum.position=[1,2];

pard.t2.object=struct('Style','text','String','BG \xi');
pard.t2.position=[3,1];
pard.t2.Width=0.5;
pard.filterx.object=struct('Style','edit','String','3');
pard.filterx.position=[3,1.5];
pard.filterx.Width=0.5;

pard.t3.object=struct('Style','text','String','BG \tau');
pard.t3.position=[4,1];
pard.t3.Width=0.5;
pard.filtert.object=struct('Style','edit','String','100');
pard.filtert.position=[4,1.5];
pard.filtert.Width=0.5;

pard.t4.object=struct('Style','text','String','kernelSize');
pard.t4.position=[2,2];
pard.t4.Width=0.5;
pard.sizeRoiSum.object=struct('Style','edit','String','3');
pard.sizeRoiSum.position=[2,2.5];
pard.sizeRoiSum.Width=0.5;
% int
% pard.checkint.object=struct('Style','checkbox','String','intensity');
% pard.checkint.position=[1,3];

% fit
pard.checkfit.object=struct('Style','checkbox','String','fit','Value',1);
pard.checkfit.position=[1,3.5];

pard.t9.object=struct('Style','text','String','Roi size fit (pix)');
pard.t9.position=[2,3.5];
% pard.t9.Width=0.5;
pard.sizeRoiFit.object=struct('Style','edit','String','7');
pard.sizeRoiFit.position=[2,4.5];
pard.sizeRoiFit.Width=0.5;

pard.fitUsePSF.object=struct('Style','checkbox','String','use fitted PSF','Value',1);
pard.fitUsePSF.position=[3,3.5];

pard.t10.object=struct('Style','text','String','Size PSF (nm)');
pard.t10.position=[4,3.5];
% pard.t10.Width=0.5;
pard.PSFxnm.object=struct('Style','edit','String','150');
pard.PSFxnm.position=[4,4.5];
pard.PSFxnm.Width=0.5;

pard.fitOnBg.object=struct('Style','checkbox','String','fit on BG','Value',1);
pard.fitOnBg.position=[5,3.5];

% 
% % 2roi
% pard.check2roi.object=struct('Style','checkbox','String','2 ROI');
% pard.check2roi.position=[1,4];
% pard.t5.object=struct('Style','text','String','size in');
% pard.t5.position=[2,4];
% pard.t5.Width=0.5;
% pard.size2roi1.object=struct('Style','edit','String',3);
% pard.size2roi1.position=[2,4.5];
% pard.size2roi1.Width=0.5;
% pard.t6.object=struct('Style','text','String','size out');
% pard.t6.position=[3,4];
% pard.t6.Width=0.5;
% pard.size2roi2.object=struct('Style','edit','String',5);
% pard.size2roi2.position=[3,4.5];
% pard.size2roi2.Width=0.5;

pard.Tfile.object=struct('Style','edit','String','settings/temp/temp_T.mat');
pard.Tfile.position=[8,1];
pard.Tfile.Width=3;

pard.loadbutton.object=struct('Style','pushbutton','String','load');
pard.loadbutton.position=[8,4];


end