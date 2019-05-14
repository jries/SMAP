classdef Roi_bg<interfaces.WorkflowModule
    properties
%         imbuffer
        imindex=1;
        bufferlength=100;
    end
    methods
        function obj=Roi_bg(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=2; 
            %tiff loader 
            % ROI
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
        end
        function prerun(obj,p)
            obj.imindex=1;
           
        end
        function outputdat=run(obj,data,p)
            if ~p.calculatebg
                outputdat=data{1};
                return
            end
            persistent imbufferlocal            
            roi=data{1}.data;%get;
            camimg=data{2}.data;
            scamimg=size(camimg);
            imindex=obj.imindex;
            if imindex==1
                imbufferlocal=zeros(scamimg(1),scamimg(2),obj.bufferlength,'single');
            end
            indxhere=mod(imindex-1,obj.bufferlength)+1;
            bufferfilled=indxhere>imindex;
            if isempty(camimg)
                outputdat=data{1};
                return
            end
            obj.imindex=obj.imindex+1;
            imbufferlocal(:,:,indxhere)=camimg;
            if ~isempty(roi)     
                sroi=size(roi.img);
                wx=(sroi(1)-1)/2;
                xpos=roi.info.x;ypos=roi.info.y;
                if ~bufferfilled
                    framerange=max(1,indxhere-p.numframes_bg-1):max(indxhere-1,1);
                else
                    framerange=mod(indxhere-p.numframes_bg-1:indxhere-1,obj.bufferlength);
                end
                for k=length(xpos):-1:1
                    roih=imbufferlocal(ypos(k)-wx:ypos(k)+wx,xpos(k)-wx:xpos(k)+wx,framerange);
                    roi.info.bgim(k,1)=single(getbackground(roih,p));
                end
               
                dato=data{1};
                dato.data=roi;
                outputdat=dato;
            else
                outputdat=data{1};
            end
        end
        

    end
end



function pard=guidef
pard.calculatebg.object=struct('Style','checkbox','String','calculate BG');
pard.calculatebg.position=[1,1];

pard.bgfunction.object=struct('Style','popupmenu','String',{{'quantile'}});
pard.bgfunction.position=[1,2];
pard.bgfunctionpar.object=struct('Style','edit','String','0.5');
pard.bgfunctionpar.position=[1,3];
pard.bgfunctionpar.Width=0.5;

pard.t1.object=struct('Style','text','String','# frames used for BG calculation');
pard.t1.position=[2,1];
pard.t1.Width=2;
pard.numframes_bg.object=struct('Style','edit','String','20');
pard.numframes_bg.position=[2,3];
pard.numframes_bg.Width=0.5;

pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='This plugin cuts out regions of interest of a defined size around the candidate positions and passes these on to the fitter';
end

function bg=getbackground(roi,p)
    bg=myquantilefast(roi(:),p.bgfunctionpar);
end