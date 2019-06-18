classdef WorkflowFitter<interfaces.WorkflowModule
    properties
%         imagestack
%         stackinfo
        stackind
        numberInBlock=0;
        newID=1;
        fittedlocs=0;
        spatial3Dcal=false;
        spatialXrange={[-inf inf]};
        spatialYrange={[-inf inf]};
        infofields={'x','y'};
    end
    methods

        function obj= WorkflowFitter(varargin)
            obj@interfaces.WorkflowModule(varargin{:});
        end
        function pard=guidef(obj)
            pard=guidef;
        end

        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
            obj.inputParameters={'loc_ROIsize','loc_cameraSettings'};
            obj.stackind=0;

        end
        function fitinit(obj) %dummy function
        end
        function prerun(obj,p)
            global fitterstackinfo fitterimagestack fitterbgstack
            fitterstackinfo=[];
%             obj.numberInBlock=1; %round(5500*100/roisize^2);
            obj.fitinit;
%             roisize=obj.getPar('loc_ROIsize');
% %             obj.numberInBlock=round(5500*100/roisize^2);
%             numx=length(obj.spatialXrange)-1;
%             numy=length(obj.spatialYrange)-1;
% %             disp(['number in block: ' num2str(obj.numberInBlock)]);
%             obj.stackind=zeros(numx,numy);
%             obj.fittedlocs=0;
% %             if obj.spatial3Dcal
%             fitterstackinfo=[];
%                 for X=numx:-1:1
%                     for Y=numy:-1:1
%                         fitterimagestack{X,Y}=zeros(roisize,roisize,obj.numberInBlock,'single');
%                         fitterstackinfo{X,Y}.x=zeros(obj.numberInBlock,1,'single');
%                         fitterstackinfo{X,Y}.y=zeros(obj.numberInBlock,1,'single');
%                         fitterstackinfo{X,Y}.frame=zeros(obj.numberInBlock,1,'double');
% %                         fitterstackinfo{X,Y}.X=X;
% %                         fitterstackinfo{X,Y}.Y=Y;
%                     end
%                 end
%                        
% %             else
% %                 fitterimagestack{1,1}=zeros(roisize,roisize,obj.numberInBlock,'single');
% %             end
%             if obj.inputChannels==2
%                 fitterbgstack=fitterimagestack;
%             else
%                 fitterbgstack=[];
%             end
            
%             infos=struct('x',0,'y',0,'frame',0);
%             fitterstackinfo=infos;

%             obj.fitinit;
            
            
            
%             obj.setInputChannels(2,'frame');
%             p=obj.getGuiParameters.par;
%             obj.loc_ROIsize=p.loc_ROIsize;
           
        end
        function out=run(obj,data,p)  
            
            global fitterimagestack fitterstackinfo fitterbgstack 
            
            if ~iscell(data) 
                dstruc=data.data;
                eof=data.eof;
            else
                dstruc=data{1}.data;%get;
                eof=data{1}.eof;
                
            end
            if isempty(fitterstackinfo)
                if isempty(dstruc)
%                     obj.status('error: image empty');drawnow
%                     error ('image empty')
                    out=[];
                    return
                else
                initstacks(obj,dstruc.info)
                end
            end
            
            persistent reporttimer
            out=[];
            passbg=(obj.inputChannels==2);%~isempty(fitterbgstack);
            fninfo=fieldnames(fitterstackinfo{1});

            
            if isempty(reporttimer)
                reporttimer=tic;
            end
            
            xrange=obj.spatialXrange;
            yrange=obj.spatialYrange;
            sxy=size(xrange);
            if ~isempty(dstruc)&&~isempty(dstruc.img)       
                 imgstack=dstruc.img;
                 
                 stackinf=dstruc.info;
%                  fninfo=fieldnames(fitterstackinfo{1});
                 s=size(imgstack);
                 if length(s)==2
                     s(3)=1;
                 end                 
                 
                 obj.fittedlocs=s(3)+obj.fittedlocs;
                 if toc(reporttimer)>2
                     obj.setPar('fittedLocs',obj.fittedlocs);
                     reporttimer=tic;
                 end
                
                if passbg
                    bgstack=data{2}.data.img;
                end
                
                numberInBlockh=obj.numberInBlock;
                
                stackindh=obj.stackind;%pointer to last element
                stackindh=stackindh+1; %new pointer
                 
                
                for X=1:sxy(1)
                    for Y=1:sxy(2)
                        if obj.spatial3Dcal&&numberInBlockh>1  %later: dont rearrange, but use instack pointer.
                            inblock=find(stackinf.x>=xrange{X,Y}(1) & stackinf.x<=xrange{X,Y}(2) & stackinf.y>=yrange{X,Y}(1) & stackinf.y<=yrange{X,Y}(2)); 
%                             inblock=find(stackinf.x>=xrange(X) & stackinf.x<=xrange(X+1) & stackinf.y>=yrange(Y) & stackinf.y<=yrange(Y+1)); 
                        else
                            inblock=1:length(stackinf.x);
                        end
                        fitstack
%                         if eof
%                             fiteof
%                         end
                    end
                end
            end
            if eof
             for X=1:sxy(1)
                for Y=1:sxy(2)
                    fiteof
                end
             end
            
            outputlocs(obj,[],[],obj.newID,true);
            end
            
            function fitstack
                if numberInBlockh>1 %make blocks
%                  sh=size(imgstack);
%                  if length(sh)==2
%                      sh(3)=1;
%                  end

                 %avoid loop
                 imagesleft=length(inblock);
                 startinstack=1;

                 newstackind=stackindh(X,Y)+imagesleft-1;
                 while newstackind>numberInBlockh
                     imagestowrite=numberInBlockh-stackindh(X,Y)+1;

                     %images
                     fitterimagestack{X,Y}(:,:,stackindh(X,Y):end)=imgstack(:,:,inblock(startinstack:startinstack+imagestowrite-1));
                     if passbg %background
                        fitterbgstack{X,Y}(:,:,stackindh(X,Y):end)=bgstack(:,:,inblock(startinstack:startinstack+imagestowrite-1));
                     end
                     for fsi=1:length(fninfo) %stackinfor
                        fitterstackinfo{X,Y}.(fninfo{fsi})(stackindh(X,Y):end)=stackinf.(fninfo{fsi})(inblock(startinstack:startinstack+imagestowrite-1));
                     end

                     %do fitting
                     stackinfoout=copyfields(fitterstackinfo{X,Y},struct('X',X,'Y',Y));
                     if passbg
                        locs=obj.fit(fitterimagestack{X,Y},fitterbgstack{X,Y},stackinfoout);
                     else
                         locs=obj.fit(fitterimagestack{X,Y},stackinfoout);
                     end
%                      [locs.asymmetry,locs.asymmdiag,locs.asymangle]=asymmetry(fitterimagestack{X,Y},true);
                     outputlocs(obj,locs,fitterstackinfo{X,Y},obj.newID,false);
                     obj.newID=obj.newID+1;

                     imagesleft=imagesleft-imagestowrite;                   
                     startinstack=startinstack+imagestowrite;         
                     stackindh(X,Y)=1;
                     newstackind=imagesleft;
                 end
                 fitterimagestack{X,Y}(:,:,stackindh(X,Y):stackindh(X,Y)+imagesleft-1)=imgstack(:,:,inblock(startinstack:end));
                 if passbg
                    fitterbgstack{X,Y}(:,:,stackindh(X,Y):stackindh(X,Y)+imagesleft-1)=bgstack(:,:,inblock(startinstack:end));
                 end
                 for fsi=1:length(fninfo) %stackinfor
                        fitterstackinfo{X,Y}.(fninfo{fsi})(stackindh(X,Y):stackindh(X,Y)+imagesleft-1)=stackinf.(fninfo{fsi})(inblock(startinstack:end));
                 end
                 obj.stackind(X,Y)=stackindh(X,Y)+imagesleft-1;

                else %go framewise
                     if passbg
                        locs=obj.fit(imgstack,bgstack,stackinf);
                     else
                        locs=obj.fit(imgstack,stackinf);
                     end
%                      [locs.asymmetry,locs.asymmdiag,locs.asymangle]=asymmetry(imgstack,true);
                     outputlocs(obj,locs,stackinf,obj.newID,false);
                     obj.newID=obj.newID+1;
                end
                
            end
            
            function fiteof
%                 if eof
                    if obj.numberInBlock>1
                        for fsi=1:length(fninfo) %stackinfor
                            fitterstackinfo{X,Y}.(fninfo{fsi})=fitterstackinfo{X,Y}.(fninfo{fsi})(1:obj.stackind(X,Y));
                        end
                        stackinfoout=copyfields(fitterstackinfo{X,Y},struct('X',X,'Y',Y));
                        if passbg
                            locs=obj.fit(fitterimagestack{X,Y}(:,:,1:obj.stackind(X,Y)),fitterbgstack{X,Y}(:,:,1:obj.stackind(X,Y)),stackinfoout);
                        else
                            locs=obj.fit(fitterimagestack{X,Y}(:,:,1:obj.stackind(X,Y)),stackinfoout);
                        end
%                         if obj.stackind(X,Y)>0
%                         [locs.asymmetry,locs.asymmdiag,locs.asymangle]=asymmetry(fitterimagestack{X,Y}(:,:,1:obj.stackind(X,Y)),true);
%                         end
                        outputlocs(obj,locs,fitterstackinfo{X,Y},obj.newID,false);
                        obj.newID=obj.newID+1;
%                     else
%                         outputlocs(obj,[],[],obj.newID,true);
                    end
%                 end
            end

        end
       
    end
end

function initstacks(obj,info)
global fitterstackinfo fitterimagestack fitterbgstack

 roisize=obj.getPar('loc_ROIsize');
%             obj.numberInBlock=round(5500*100/roisize^2);
            sxy=size(obj.spatialXrange);
            numx=sxy(1);numy=sxy(2);
%             numx=length(obj.spatialXrange)-1;
%             numy=length(obj.spatialYrange)-1;
%             disp(['number in block: ' num2str(obj.numberInBlock)]);
            obj.stackind=zeros(numx,numy);
            obj.fittedlocs=0;
%             if obj.spatial3Dcal
            fitterstackinfo=[];
            addfields=setdiff(fieldnames(info),[obj.infofields,{'frame'}]);
                for X=numx:-1:1
                    for Y=numy:-1:1
                        fitterimagestack{X,Y}=zeros(roisize,roisize,obj.numberInBlock,'single');
                        fitterstackinfo{X,Y}.frame=zeros(obj.numberInBlock,1,'double');
                        for k=1:length(obj.infofields)
                            fitterstackinfo{X,Y}.(obj.infofields{k})=zeros(obj.numberInBlock,1,'single');
                        end
                        for k=1:length(addfields)
                            if iscell(info.(addfields{k}))
                                fitterstackinfo{X,Y}.(addfields{k}){obj.numberInBlock,1}=[];
                            else
                            fitterstackinfo{X,Y}.(addfields{k})=zeros(obj.numberInBlock,1,'single');
                            end
                        end
         
                    end
                end
                       
%             else
%                 fitterimagestack{1,1}=zeros(roisize,roisize,obj.numberInBlock,'single');
%             end
            if obj.inputChannels==2
                fitterbgstack=fitterimagestack;
            else
                fitterbgstack=[];
            end
end

function outputlocs(obj,locs,stackinfo,tag,eof)
dato=interfaces.WorkflowData;
dato.ID=tag;
dato.eof=eof;
dato.data=(locs);
obj.output(dato);

end

function pard=guidef

pard=[];

end