classdef GrabFijiStacks<interfaces.WorkflowModule
% Opens an instance of Fiji in which you can open any image stack (if
% large, use virtual stack). This image stack can then be selected in the
% GUI and used for fitting.As metadata is not parsed, set it manually in
% the Camera Converter
    properties     
        framestop
        framestart
        numberOfFrames
        timerfitstart
        windows
        currentImage
    end
    methods
        function obj=GrabFijiStacks(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=0;
            obj.isstartmodule=true;
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
             obj.inputParameters={'loc_subtractbg','loc_blocksize_frames'};            
        end
        function prerun(obj,p)
            p=obj.getAllParameters;
            
            if p.locdata_empty
                obj.locData.clear;
            end
            

            obj.framestop=p.framestop;
            obj.framestart=p.framestart;
            if obj.getPar('loc_preview')
                previewframe=obj.getPar('loc_previewframe');               
                if p.loc_subtractbg
                    dt=p.loc_blocksize_frames;
                    frameload=max(1,previewframe-floor(dt/2));
                    obj.framestart=frameload;
                    obj.currentImage=frameload;
                    obj.framestop=frameload+dt;
                else
                   obj.currentImage=previewframe;
                    obj.framestop=previewframe;
                end               
            else
                obj.currentImage=obj.framestart;
            end
            obj.timerfitstart=tic;            
        end
        function run(obj,data,p)
            global SMAP_stopnow        
            if p.filelist.Value>length(obj.windows)
                disp('update filelist first')
                return
            end
            stack=obj.windows(p.filelist.Value).stack;
            id=1;
            sizeim=obj.windows(p.filelist.Value).size;
            image=readstack(stack,obj.currentImage,sizeim);   
            
            while ~isempty(image)&&obj.currentImage<=obj.framestop&&~SMAP_stopnow                
                datout=interfaces.WorkflowData;
                datout.data=image;
                datout.frame=obj.currentImage;
                datout.ID=id;
                id=id+1;
                obj.output(datout)
                obj.currentImage=obj.currentImage+1;
                image=readstack(stack,obj.currentImage,sizeim);   
                %display
                if mod(datout.frame,10)==0
                    totalf=min(sizeim(3),obj.framestop)-obj.framestart;
                    elapsed=toc(obj.timerfitstart);
                    totaltime=elapsed/(datout.frame-obj.framestart+1)*totalf;
                    statuss=['frame ' int2str(datout.frame-obj.framestart) ' of ' int2str(totalf) ...
                        '. Time: ' num2str(elapsed,'%4.0f') ' of ' num2str(totaltime,'%4.0f') 's'];
                    obj.status(statuss);
                end
            end
            dateof=interfaces.WorkflowData;
            dateof.frame=obj.currentImage+1;
            dateof.ID=id;
            dateof.eof=true;
            obj.output(dateof)
            disp('fitting done')
        end

        function fijibutton_callback(obj,a,b)
            openfiji(obj);
            obj.status('import images in fiji as (virtual) stack');
        end
        function refreshfilelist_callback(obj,a,b)
            ijh=obj.getPar('IJ');
            ij=ijh.getInstance();          
            frames=ij.getFrames;
            obj.windows=[];
            ind=1;
            for k=1:length(frames)
                if strcmp(frames(k).class,'ij.gui.StackWindow')&&~isempty(frames(k).getImagePlus)
                    obj.windows(ind).name=char(frames(k).getTitle);
                    obj.windows(ind).number=k;
                    stack=frames(k).getImagePlus.getStack;
                    obj.windows(ind).stack=stack;
                    obj.windows(ind).imagePlus=frames(k).getImagePlus;
                    obj.windows(ind).size=[stack.getWidth stack.getHeight stack.getSize];
                    ind=ind+1;
                end
            end
            obj.selectfilelist_callback(0,0);
        end 
        function info=makefilestruct(obj,window)
            ip=window.imagePlus;
            of=ip.getOriginalFileInfo;
            fdir=char(of.directory);
            totalfile=[fdir char(of.fileName)];
            finf=obj.getPar('loc_fileinfo');
            il=imageloaderAll(totalfile,finf,obj.P);
            info=il.metadata;
            info.numberOfFrames=obj.framestop;
            info.Width=window.size(1);
            info.Height=window.size(2);
            obj.setPar('loc_fileinfo',info);        
        end
        
        function selectfilelist_callback(obj,a,b)
            if isempty(obj.windows)
                stacklist={'no open stacks found'};
            else
                stacklist={obj.windows(:).name};
            end
            obj.guihandles.filelist.String=stacklist;
            obj.guihandles.filelist.Value=max(1,min(length(stacklist),obj.guihandles.filelist.Value));
            vs=obj.guihandles.filelist.Value;
            if ~isempty(obj.windows)
                obj.framestop=obj.windows(vs).size(3);
                obj.guihandles.framestop.String=num2str(obj.framestop);
                fi=obj.makefilestruct(obj.windows(vs));
                obj.setPar('loc_fileinfo',fi);
            end
        end
    end
end



function image=readstack(stack,number,ss)
if nargin<3||isempty(ss)
    ss=[stack.getWidth stack.getHeight stack.getSize];
end
if number>0&&number<=ss(3)
    pixel=stack.getPixels(number);
    image=reshape(pixel,ss(1),ss(2))';
else
    image=[];
end
end
function pard=guidef(obj)
pard.text.object=struct('Style','text','String','and open images as virtual stack');
pard.text.position=[1,2];
pard.text.Width=2;

pard.fijibutton.object=struct('Style','pushbutton','String','Open Fiji','Callback',@obj.fijibutton_callback);
pard.fijibutton.position=[1,1];

pard.refreshfiles.object=struct('Style','pushbutton','String','Refresh list: ','Callback',@obj.refreshfilelist_callback);
pard.refreshfiles.position=[2,1];

pard.filelist.object=struct('Style','popupmenu','String',{'empty'},'Callback',@obj.selectfilelist_callback);
pard.filelist.position=[2,2];
pard.filelist.Width=2;

pard.fitall.object=struct('Style','checkbox','String','Fit all','Value',0);
pard.fitall.position=[2,4];

pard.textf.object=struct('Style','text','String','Frame range: ');
pard.textf.position=[4.2,1.25];
pard.textf.Width=0.75;
pard.framestart.object=struct('Style','edit','String','1');
pard.framestart.position=[4.2,2];
pard.framestart.Width=0.7;
pard.framestop.object=struct('Style','edit','String','1000000');
pard.framestop.position=[4.2,2.7];
pard.framestop.Width=0.7;

pard.locdata_empty.object=struct('Style','checkbox','String','Empty localizations','Value',1);
pard.locdata_empty.position=[4.2,3.5];
pard.locdata_empty.Width=1.5;
pard.locdata_empty.TooltipString='empty localization data before fitting. Important if post-processing (eg drift correction) is perfromed as part of workflow';
pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='Opens an instance of Fiji in which you can open any image stack (if large, use virtual stack). This image stack can then be selected in the GUI and used for fitting.';
end