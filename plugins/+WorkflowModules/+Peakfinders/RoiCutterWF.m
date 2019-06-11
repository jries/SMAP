classdef RoiCutterWF<interfaces.WorkflowModule
%     This plugin cuts out regions of interest of a defined size around the
%     candidate positions and passes these on to the fitter.
    properties
        loc_ROIsize
        preview
        disppreview
    end
    methods
        function obj=RoiCutterWF(varargin)
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
            
            p=obj.getAllParameters;
            obj.loc_ROIsize=p.loc_ROIsize;
            obj.preview=obj.getPar('loc_preview');
            obj.setPar('loc_ROIsize',p.loc_ROIsize);
            obj.disppreview=false;
           
        end
        function outputdat=run(obj,data,p)
            outputdat=[];
            image=data{1}.data;%get;
            if ~isempty(image)     
            maxima=data{2}.data;%get;
            if isempty(maxima.x)
                if obj.preview
                    obj.status('no localizations found');drawnow
                    error ('no localizations found')

                end
%                  outs.info=info;
%                     outs.img=cutoutimages(:,:,1:ind);
                dato=data{1};%.copy;
                dato.data=[];
                outputdat=dato;
                return;
            end
            
            kernelSize=obj.loc_ROIsize;
            dn=ceil((kernelSize-1)/2);
            sim=size(image);
            
            if p.loc_filterforfit>0
                %test: filter image befor fitting for z from 2D
                %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                h=fspecial('gaussian',3,p.loc_filterforfit);

                image=filter2(h,image);
            %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
            end
            
            cutoutimages=zeros(kernelSize,kernelSize,length(maxima.x),'single');
            ind=0;
            goodind=~(maxima.y<=dn|maxima.y>sim(1)-dn|maxima.x<=dn|maxima.x>sim(2)-dn);
            outside=(maxima.y<=1|maxima.y>sim(1)-1|maxima.x<=1|maxima.x>sim(2)-1);

            for k=1:length(maxima.x)
                 ind=ind+1;
                if goodind(k)
                    cutoutimages(:,:,ind)=image(maxima.y(k)-dn:maxima.y(k)+dn,maxima.x(k)-dn:maxima.x(k)+dn); %coordinates exchanged.
                elseif outside(k)
                    maxima.y(k)=max(dn+1,min(maxima.y(k),sim(1)-dn));
                    maxima.x(k)=max(dn+1,min(maxima.x(k),sim(2)-dn));  
                    cutoutimages(:,:,ind)=zeros(2*dn+1,'single');
                else
                    maxima.y(k)=max(dn+1,min(maxima.y(k),sim(1)-dn));
                    maxima.x(k)=max(dn+1,min(maxima.x(k),sim(2)-dn));
                    cutoutimages(:,:,ind)=image(maxima.y(k)-dn:maxima.y(k)+dn,maxima.x(k)-dn:maxima.x(k)+dn); %coordinates exchanged.
%                     cutoutimages(:,:,ind)=zeros(kernelSize,kernelSize,1,'single');
                end  
            end 
            info=maxima;
            frameh=data{1}.frame;
            info.frame=maxima.x*0+frameh;
           

            outs.info=info;
            outs.img=cutoutimages(:,:,1:ind);
            dato=data{1};%.copy;
            dato.data=outs;%set(outs);
            outputdat=dato;

            if obj.preview 
                obj.disppreview=true;
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
                if obj.preview && ~obj.disppreview
                    obj.status('image could not be loaded');drawnow
                    error('no image loaded')
                end
                outputdat=data{1};
            end
        end
    end
end



function pard=guidef
pard.text.object=struct('Style','text','String','Size ROI (pix)');
pard.text.position=[1,1];
pard.text.Width=0.8;

pard.loc_ROIsize.object=struct('Style','edit','String','7');
pard.loc_ROIsize.position=[1,1.8];
pard.loc_ROIsize.Width=0.4;
pard.loc_ROIsize.TooltipString=sprintf('Size (pixels) of regions around each peak candidate which are used for fitting. \n Depends on fitter. Use larger ROIs for 3D data.');

pard.loc_filterforfit.object=struct('Style','edit','String','0');
pard.loc_filterforfit.position=[1,2.2];
pard.loc_filterforfit.TooltipString=sprintf('Filter before fit. Sigma of Gaussian kernel in pixels (0: no filter).');
pard.loc_filterforfit.Width=0.3;

pard.syncParameters={{'loc_ROIsize','loc_ROIsize',{'String'}},{'loc_filterforfit','loc_filterforfit',{'String'}}};

pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='This plugin cuts out regions of interest of a defined size around the candidate positions and passes these on to the fitter';
end