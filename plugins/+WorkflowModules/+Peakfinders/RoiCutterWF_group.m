classdef RoiCutterWF_group<interfaces.WorkflowModule
    properties
        loc_ROIsize
        preview
        temprois
        tempinfo
        dX
        dT
        dn
    end
    methods
        function obj=RoiCutterWF_group(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=2; 
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
            obj.addSynchronization('loc_ROIsize',obj.guihandles.loc_ROIsize,'String')
            obj.setInputChannels(2,'frame');

        end
        function prerun(obj,p) 
%              global tempinfo temprois
            p=obj.getAllParameters;
            obj.loc_ROIsize=p.loc_ROIsize;
            obj.preview=obj.getPar('loc_preview');
            obj.setPar('loc_ROIsize',p.loc_ROIsize);
            obj.dX=p.loc_dTdx(2);
            obj.dT=p.loc_dTdx(1);
            kernelSize=obj.loc_ROIsize;
            obj.dn=ceil((kernelSize-1)/2);
            obj.run; %initialize
%             ninit=100;
%             init=zeros(ninit,1);
%             tempinfo=struct('inuse',false(size(init)),'x',init,'y',init,'dT',init,'numrois',init,'inframes',{{}});
%             temprois=zeros(obj.loc_ROIsize,obj.loc_ROIsize,ninit,'single');
           
        end
        function outputdat=run(obj,data,p)
            persistent tempinfo temprois tempdT numrois tempx tempy inuse inframes
            if nargin==1 %reset
                ninit=100;
                init=zeros(ninit,1);
%                 tempinfo=struct('inuse',false(size(init)),'x',init,'y',init,'dT',init,'numrois',init,'inframes',{{}});
                temprois=zeros(obj.loc_ROIsize,obj.loc_ROIsize,ninit,'single');
                inframes=zeros(ninit,3);
                tempdT=init;
                numrois=init; tempx=init; tempy=init; inuse=false(size(init));
                return
            end
            if ~p.loc_group
                outputdat=run_nogroup(obj,data,p);
                return
            end
        

            %for challange: fit two times (with/wo this) or include info
            %about frames 'on' and then resdistribute.
            %this is anyways needed e.g. for color assignment
            outputdat=[];
            image=data{1}.data;%get;
            frameh=data{1}.frame;
            if ~isempty(image)     
            maxima=data{2}.data;%get;
            if isempty(maxima.x)
                return;
            end
            

            dX=obj.dX;dT=obj.dT;dn=obj.dn;
            for k=1:length(maxima.x)
                addroi(image,maxima.x(k),maxima.y(k),frameh,dX,dT,dn)
            end 
            
          
            
            if obj.preview 
                outputfig=obj.getPar('loc_outputfig');
                if ~isvalid(outputfig)
                    outputfig=figure(209);
                    obj.setPar('loc_outputfig',outputfig);
                end
                outputfig.Visible='on';
                figure(outputfig)
                hold on
                col=[0.3 0.3 0.3];
                    ax=gca;
                for k=1:length(maxima.x)
                    pos=[maxima.x(k)-obj.dn maxima.y(k)-obj.dn maxima.x(k)+obj.dn maxima.y(k)+obj.dn ];
                    plotrect(ax,pos,col);
                end
                [cutoutimages,maximap]=purgeall();
            else
                [cutoutimages,maximap]=purgerois();
            end 

            if ~isempty(cutoutimages)
                info=maximap;
                info.frame=maximap.x*0+frameh;
                outs.info=info;
                outs.img=cutoutimages;
                dato=data{1};%.copy;
                dato.data=outs;%set(outs);
                outputdat=dato; 
            end
            
            else
                outputdat=data{1};
            end
        
        function addroi(image,x,y,frame,dX,dT,dn)
            maxgroup=100;
%             global tempinfo temprois
%             tempinfo=obj.tempinfo;
            %later: keep x0 as cut out position (has to be same), but
            %update xsearch
%             dX=obj.dX;
%             inxy=find((tempinfo.x(inuse)-x).^2+(tempinfo.y(inuse)-y).^2<dX^2);
            
            indtemp= find(inuse & tempx>x-dX & tempx<x+dX & tempy>y-dX & tempy<y+dX,1,'first');
            if ~isempty(indtemp) && numrois(indtemp)<maxgroup %already there
%                 finuse=find(inuse);
%                 indtemp=finuse(inxy(1)); %later: choose closest

                temprois(:,:,indtemp)=temprois(:,:,indtemp)+cutoutimage(image,tempx(indtemp),tempy(indtemp),dn);
                %fill info
%                 tempinfo.dT(indtemp)=dT; %reset
               
                numrois(indtemp)=numrois(indtemp)+1;
%                 tempinfo.inframes{indtemp}(end+1)=frame;
                inframes(indtemp,numrois(indtemp))=frame; 
                tempdT(indtemp)=dT;
              
            else %new ROI, this would be standard
                newind=find(inuse==false,1,'first');
                if isempty(newind)
                    newind=length(inuse)+1;
                end
                temprois(:,:,newind)=cutoutimage(image,x,y,dn);
%                 tempinfo.dT(newind)=dT; %reset
                numrois(newind)=1;
                tempx(newind)=x;tempy(newind)=y;  
                inuse(newind)=true;
%                 tempinfo.inframes{newind}=frame;
                inframes(newind,1)=frame;
                tempdT(newind)=dT;
            end
        end

        function [cutoutimages,maximap]=purgerois()
            
%               global tempinfo temprois

            finuse=find(inuse);
            indout=tempdT(finuse)<0;
            tempdT(finuse)=tempdT(finuse)-1;
            fout=finuse(indout);
            if isempty(fout)
                cutoutimages=[];maximap=[];
                return
            end
            
            cutoutimages=temprois(:,:,fout);
            maximap.x=tempx(fout);
            maximap.y=tempy(fout);
            maximap.inframes=mat2arrayhere(inframes,numrois,fout);
%             maximap.inframes=inframes(fout,:); %XXX
            inuse(fout)=false;
%             tempinfo.dT(finuse)=tempinfo.dT(finuse)-1;% count down dark frames
            
        end
         function [cutoutimages,maximap]=purgeall()
%              global tempinfo temprois
            finuse=find(inuse);
            fout=finuse;
            if isempty(fout)
                cutoutimages=[];maximap=[];
                return
            end
            cutoutimages=temprois(:,:,fout);
            maximap.x=tempx(fout);
            maximap.y=tempy(fout);
            maximap.inframes=mat2arrayhere(inframes,numrois,fout);
%             maximap.inframes=tempinfo.inframes(fout);
            
            inuse(fout)=false;
%             tempinfo.dT(finuse)=tempinfo.dT(finuse)-1;% count down dark frames
         end
        end
    end
end
function out=mat2arrayhere(inframes,numrois,fout)
out={};
for k=length(fout):-1:1
    out{k}=inframes(fout(k),1:numrois(k));
end
end

        function out=cutoutimage(image, x,y,dn)
            sim=size(image);
                if  (y<=dn||y>sim(1)-dn||x<=dn||x>sim(2)-dn)
                    y=max(dn+1,min(y,sim(1)-dn));
                    x=max(dn+1,min(x,sim(2)-dn));
                end  
                out=image(y-dn:y+dn,x-dn:x+dn);
        end
        
        
        
function outputdat=run_nogroup(obj,data,p) %from RoiCutterWF
            outputdat=[];
            image=data{1}.data;
            if ~isempty(image)     
            maxima=data{2}.data;
            if isempty(maxima.x)
                return;
            end
            
            kernelSize=obj.loc_ROIsize;
            dn=ceil((kernelSize-1)/2);
            sim=size(image);
            
%             if p.loc_filterforfit>0 && length(p.loc_filterforfit)==1
%                 %test: filter image befor fitting for z from 2D
%                 %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%                 h=fspecial('gaussian',3,p.loc_filterforfit);
%                 image=filter2(h,image);
%             %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%             end
            
            cutoutimages=zeros(kernelSize,kernelSize,length(maxima.x),'single');
            ind=0;
            goodind=~(maxima.y<=dn|maxima.y>sim(1)-dn|maxima.x<=dn|maxima.x>sim(2)-dn);

            for k=1:length(maxima.x)
                 ind=ind+1;
                if goodind(k)
                    cutoutimages(:,:,ind)=image(maxima.y(k)-dn:maxima.y(k)+dn,maxima.x(k)-dn:maxima.x(k)+dn); %coordinates exchanged.
                else
                    maxima.y(k)=max(dn+1,min(maxima.y(k),sim(1)-dn));
                    maxima.x(k)=max(dn+1,min(maxima.x(k),sim(2)-dn));
                    cutoutimages(:,:,ind)=image(maxima.y(k)-dn:maxima.y(k)+dn,maxima.x(k)-dn:maxima.x(k)+dn); %coordinates exchanged.
                end  
            end 
            info=maxima;
            frameh=data{1}.frame;
            info.frame=maxima.x*0+frameh;
            outs.info=info;
            outs.img=cutoutimages(:,:,1:ind);
            dato=data{1};%.copy;
            dato.data=outs;%set(outs);
%             obj.output(dato)
            outputdat=dato;

            
            
            
            if obj.preview 
                outputfig=obj.getPar('loc_outputfig');
                if ~isvalid(outputfig)
                    outputfig=figure(209);
                    obj.setPar('loc_outputfig',outputfig);
                end
                outputfig.Visible='on';
                figure(outputfig)
                hold on
                    col=[0.3 0.3 0.3];
                    ax=gca;
                for k=1:length(maxima.x)
                    pos=[maxima.x(k)-dn maxima.y(k)-dn maxima.x(k)+dn maxima.y(k)+dn ];
                    
                    plotrect(ax,pos,col);
                end
            end 
            else
                outputdat=data{1};
            end
end
            
function pard=guidef
pard.text.object=struct('Style','text','String','ROI (pix):');
pard.text.position=[1,1];
pard.text.Width=0.7;

pard.loc_ROIsize.object=struct('Style','edit','String','7');
pard.loc_ROIsize.position=[1,1.5];
pard.loc_ROIsize.Width=0.5;
pard.loc_ROIsize.TooltipString=sprintf('Size (pixels) of regions around each peak candidate which are used for fitting. \n Depends on fitter. Use larger ROIs for 3D data.');



pard.loc_group.object=struct('Style','checkbox','String','group');
pard.loc_group.position=[2,1];
pard.loc_group.TooltipString=sprintf('Groupping parameters. dT (frames), dX (pixels)');
pard.loc_group.Width=0.7;

pard.loc_dTdx.object=struct('Style','edit','String','1 1.5');
pard.loc_dTdx.position=[2,1.5];
pard.loc_dTdx.TooltipString=sprintf('Groupping parameters. dT (frames), dX (pixels)');
pard.loc_dTdx.Width=0.5;

pard.syncParameters={{'loc_ROIsize','loc_ROIsize',{'String'}}};

pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='This plugin cuts out regions of interest of a defined size around the candidate positions and passes these on to the fitter';
end