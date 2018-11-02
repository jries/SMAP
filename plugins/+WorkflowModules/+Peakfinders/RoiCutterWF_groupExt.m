classdef RoiCutterWF_groupExt<interfaces.WorkflowModule
    properties
        loc_ROIsize
        preview
        temprois
        tempinfo
        dn
    end
    methods
        function obj=RoiCutterWF_groupExt(varargin)
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
             global tempinfo temprois
            p=obj.getAllParameters;
            obj.loc_ROIsize=p.loc_ROIsize;
            obj.preview=obj.getPar('loc_preview');
            obj.setPar('loc_ROIsize',p.loc_ROIsize);

            kernelSize=obj.loc_ROIsize;
            obj.dn=ceil((kernelSize-1)/2);
            ninit=1000;
            init=zeros(ninit,1);
            tempinfo=struct('inuse',false(size(init)),'x',init,'y',init,'dT',init,'numrois',init,'groupindex',init,'numberInGroup',init,'phot',init,'bg',init);
            temprois=zeros(obj.loc_ROIsize,obj.loc_ROIsize,ninit,'single');
           
        end
        function outputdat=run(obj,data,p)
            outputdat=[];
            image=data{1}.data;%get;
            if ~isempty(image)     
            maxima=data{2}.data;%get;
            if isempty(maxima.x)
                return;
            end
            for k=1:length(maxima.x)
                obj.addroi(image,copystructReduce(maxima,k))
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
                [cutoutimages,maximap]=obj.purgeall;
            else
                [cutoutimages,maximap]=obj.purgerois;
            end 
             info=maximap;
            frameh=data{1}.frame;
            info.frame=maximap.x*0+frameh;
            if ~isempty(cutoutimages)
                outs.info=info;
                outs.img=cutoutimages;
                dato=data{1};%.copy;
                dato.data=outs;%set(outs);
                outputdat=dato; 
            end
            
            else
                outputdat=data{1};
            end
        end
        function addroi(obj,image,maxima)
            %   use groupnumber and numberingroup
%   if same groupnumber: add. If numrois>=numberingroup: purge
% maxima: needs fields: .groupindex, .numberInGroup

            global tempinfo temprois
            inuse=tempinfo.inuse;
            x=maxima.x;y=maxima.y;
            %later: keep x0 as cut out position (has to be same), but
            %update xsearch
            inxy=find(tempinfo.groupindex(inuse)==maxima.groupindex);
            if ~isempty(inxy) %already there
                finuse=find(inuse);
                indtemp=finuse(inxy(1)); %later: choose closest
                temprois(:,:,indtemp)=temprois(:,:,indtemp)+cutoutimage(image,tempinfo.x(indtemp),tempinfo.y(indtemp),obj.dn);
                %fill info
                tempinfo.numrois(indtemp)=tempinfo.numrois(indtemp)+1;
                tempinfo.phot(indtemp)=tempinfo.phot(indtemp)+maxima.phot;
                tempinfo.bg(indtemp)=tempinfo.bg(indtemp)+maxima.bg;
            else %new ROI, this would be standard
                newind=find(tempinfo.inuse==false,1,'first');
                if isempty(newind)
                    newind=length(tempinfo.inuse)+1;
                end
                temprois(:,:,newind)=cutoutimage(image,x,y,obj.dn);
                tempinfo.numrois(newind)=1;
%                 tempinfo.x(newind)=x;tempinfo.y(newind)=y;  
                tempinfo.inuse(newind)=true;
                fn=fieldnames(maxima);
                for k=1:length(fn)
                    tempinfo.(fn{k})(newind)=maxima.(fn{k});
                end
                
%                 tempinfo.numberInGroup(newind)=maxima.numberInGroup;
%                 tempinfo.groupindex(newind)=maxima.groupindex;
            end
        end

        function [cutoutimages,maximap]=purgerois(obj)
              global tempinfo temprois

            finuse=find(tempinfo.inuse);
            indout=tempinfo.numrois(finuse)>=tempinfo.numberInGroup(finuse);
            fout=finuse(indout);
            
            cutoutimages=temprois(:,:,fout);
            maximap=copystructReduce(tempinfo,fout);
%             maximap.x=tempinfo.x(fout);
%             maximap.y=tempinfo.y(fout);
            
            tempinfo.inuse(fout)=false;
        end
         function [cutoutimages,maximap]=purgeall(obj)
             global tempinfo temprois
            finuse=find(tempinfo.inuse);
            fout=finuse;
            cutoutimages=temprois(:,:,fout);
            maximap.x=tempinfo.x(fout);
            maximap.y=tempinfo.y(fout);
            
            tempinfo.inuse(fout)=false;
%             tempinfo.dT(finuse)=tempinfo.dT(finuse)-1;% count down dark frames
        end
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

function pard=guidef
pard.text.object=struct('Style','text','String','Size ROI (pix)');
pard.text.position=[1,1];


pard.loc_ROIsize.object=struct('Style','edit','String','7');
pard.loc_ROIsize.position=[2,1];
pard.loc_ROIsize.Width=0.7;
pard.loc_ROIsize.TooltipString=sprintf('Size (pixels) of regions around each peak candidate which are used for fitting. \n Depends on fitter. Use larger ROIs for 3D data.');

pard.loc_fitgrouped.object=struct('Style','checkbox','String','group');
pard.loc_fitgrouped.position=[2,1];
pard.loc_fitgrouped.TooltipString=sprintf('Fit on combined images based on pre-calculated grouping');
pard.loc_fitgrouped.Width=1;

pard.syncParameters={{'loc_ROIsize','loc_ROIsize',{'String'}}};

pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='This plugin cuts out regions of interest of a defined size around the candidate positions and passes these on to the fitter';
end