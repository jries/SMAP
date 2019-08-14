classdef deepSMLM<interfaces.WorkflowModule
    properties (Access=private)
        deepobj
        mask
        imagebuffer
        bufferinfo
        buffercounter
        Xgrid
        Ygrid
        buffersize
        preview
    end
    methods
        function obj=deepSMLM(varargin)
            obj@interfaces.WorkflowModule(varargin{:});
            obj.inputParameters={'loc_loc_filter_sigma','EMon'};
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
             initGui@interfaces.WorkflowModule(obj);
%             obj.setPar('loc_subtractbg',true) %hack to have 3 images in preview
%             obj.setPar('loc_blocksize_frames',3);
        end
        function prerun(obj,p)
            gitdir=fileparts(pwd);
            deeppath= [gitdir '/ries-private/InferenceSMLM/Inference_SMLM'];
            addpath([deeppath  '/src/matlab'])
            obj.deepobj=mex_interface(str2fun([deeppath '/build/deepsmlm_interface']), p.modelfile);
            
            obj.imagebuffer=[];
            if obj.getPar('loc_preview')
                obj.buffersize=1;
            else
                obj.buffersize=p.buffersize;
            end
            obj.preview=obj.getPar('loc_preview');
            

        end
        function dato=run(obj,data,p)
            %
            dato=[];
            bufferfactor=1.2;
            xfactor=.6*bufferfactor;
            yfactor=.6*bufferfactor;
            zfactor=750*bufferfactor;
            photfactor=50000*bufferfactor;
            image=data.data;%get;
            image=single(image); %now we use directly camera frames
            if isempty(image)
                dato=data;
                dato.data=[];
                return
            end
            
            
            if isempty(obj.imagebuffer)
                obj.imagebuffer=zeros(obj.buffersize,size(image,1),size(image,2),'single');
                obj.bufferinfo.frame=zeros(obj.buffersize,1,'double');
                obj.buffercounter=0;
                nx=1:size(image,2);ny=1:size(image,1);
                [obj.Xgrid,obj.Ygrid]=meshgrid(nx,ny);
                obj.mask=obj.getPar('loc_roimask'); %run at first time
            end
            obj.buffercounter=obj.buffercounter+1;
%             bufferind=mod(obj.buffercounter-1,3)+1;
            obj.imagebuffer(obj.buffercounter,:,:)=image;
            obj.bufferinfo.frame(obj.buffercounter)=data.frame;
            if obj.buffercounter<obj.buffersize %we need to wait for buffersize images. check if eof, then process the rest.
                
                return
            end
            
%             indsort=[bufferind:3 1:bufferind-1];indsort=indsort([2 3 1]);
%             fitimage=obj.imagebuffer(1,indsort,:,:);
%             fitimage=permute(obj.imagebuffer,[1 3 2]);
            
            
            %XXXX later remove
%             fitimage=obj.imagebuffer;
%             camsettings=obj.getPar('loc_cameraSettings');
%             imbufferadu=fitimage/camsettings.pix2phot+camsettings.offset;

            imbufferadu=obj.imagebuffer; %already in ADU
            [xf, x_shape] = pseudo_col2row_major(imbufferadu);
            
            [outf, out_size] = obj.deepobj.forward(xf, x_shape);
            deepim = pseudo_row2_col_major(outf, out_size);
             deepim(:,6,:,:)=imbufferadu;
%             deepim=permute(deepim,[1 2 4 3]);
            
            for k=1:size(deepim,1) %process each frame individually
                probmap=squeeze(deepim(k,1,:,:));
                if all(size(obj.mask)==size(image)) %apply mask
                    probmap=probmap.*obj.mask;
                end
                maxima=maximumfindcall((probmap)); 
                cutoff=p.pcutoff;
                maxind= (maxima(:,3)>cutoff);
                y=maxima(maxind,1);
                x=maxima(maxind,2);
                linind=sub2ind(size(image),x,y);
                Ymap=squeeze(deepim(k,3,:,:))*yfactor+obj.Ygrid;
                Xmap=squeeze(deepim(k,4,:,:))*xfactor+obj.Xgrid;
                photmap=squeeze(deepim(k,2,:,:));
                Zmap=squeeze(deepim(k,5,:,:));
                dxmap=squeeze(deepim(k,4,:,:))*yfactor;
                dymap=squeeze(deepim(k,3,:,:))*yfactor;

                xfit=Xmap(linind);
                yfit=Ymap(linind);
                intensity=photmap(linind)*photfactor;
                zfit=Zmap(linind)*zfactor;
                pfit=probmap(linind);
                dx=dxmap(linind);
                dy=dymap(linind);

                %export for SMAP
                v1=0*xfit+1;
                locs.frame=v1*obj.bufferinfo.frame(k);
                locs.xcnn=xfit;
                locs.ycnn=yfit;
                locs.zcnn=zfit;
                locs.photcnn=intensity;
                locs.prob=pfit;
%                 locs.PSFxpix=v1;
%                 locs.bg=0*v1; %now set to zero, later determine from CNN or from photon converted image
%                 locs.xpixerr=sqrt((locs.PSFxpix.*locs.PSFxpix+1/12*v1)./( locs.phot)+8*pi*(locs.PSFxpix.*locs.PSFxpix).^2.* locs.bg./( locs.phot).^2);
%                 locs.ypixerr=locs.xpixerr;
                locs.dx=dx;
                locs.dy=dy;
                
                %to be used as a peak finder:
                locs.x=xfit;
                locs.y=yfit;
                locs.N=intensity;
                locs.znm=zfit;
                
                dato=data; 
                dato.frame=obj.bufferinfo.frame(k);
                dato.ID=dato.frame;
                dato.data=locs; 
                obj.output(dato);
            end
%             if obj.preview
%                 outputfig=obj.getPar('loc_outputfig');
%                 if ~isvalid(outputfig)
%                     outputfig=figure(209);
%                     obj.setPar('loc_outputfig',outputfig);
%                 end
%                 outputfig.Visible='on';
%                 figure(outputfig)
%                 hold off
%                 imagesc(image);
%                 colormap jet
%                 colorbar;
%                 axis equal
%                 hold on
%                 plot(xfit,yfit,'yo')
%             end
            obj.buffercounter=0;
        end
    end
end




function loadmodel_callback(a,b,obj)
file=obj.getSingleGuiParameter('modelfile');
if isempty(file)
    file='*.pt';
end
[file path]=uigetfile(file);
if file
    obj.setGuiParameters(struct('modelfile',[path  file]));
end
end

function pard=guidef(obj)
pard.txt.object=struct('Style','text','String','deepSMLM:');
pard.txt.position=[1,1];
pard.loadmodel.object=struct('Style','pushbutton','String','Load model','Callback',{{@loadmodel_callback,obj}});
pard.loadmodel.position=[1,2];
pard.modelfile.object=struct('Style','edit','String','','HorizontalAlignment','right');
pard.modelfile.position=[2,1];
pard.modelfile.Width=2;

pard.pcutofft.object=struct('Style','text','String','probability cutoff');
pard.pcutofft.position=[3,1];
pard.pcutofft.Width=1;
pard.pcutoff.object=struct('Style','edit','String','0.3');
pard.pcutoff.position=[3,2];
pard.pcutoff.Width=.35;

pard.buffersizet.object=struct('Style','text','String','buffer (frames)');
pard.buffersizet.position=[4,1];
pard.buffersizet.Width=1;
pard.buffersize.object=struct('Style','edit','String','100');
pard.buffersize.position=[4,2];
pard.buffersize.Width=.35;


pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='';
end