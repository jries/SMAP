classdef deepSMLM<interfaces.WorkflowModule
    properties (Access=private)
        deepobj
        mask
        imagebuffer
        buffercounter
        Xgrid
        Ygrid
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
            obj.setPar('loc_subtractbg',true) %hack to have 3 images in preview
            obj.setPar('loc_blocksize_frames',3);
        end
        function prerun(obj,p)
            gitdir=fileparts(pwd);
            deeppath= [gitdir '/ries-private/InferenceSMLM/Inference_SMLM'];
            addpath([deeppath  '/src/matlab'])
            obj.deepobj=mex_interface(str2fun([deeppath '/build/deepsmlm_interface']), p.modelfile);
            
            obj.imagebuffer=[];
            
            

        end
        function dato=run(obj,data,p)
            %
            bufferfactor=1.2;
            xfactor=.5*bufferfactor;
            yfactor=.5*bufferfactor;
            zfactor=750*bufferfactor;
            photfactor=10000*bufferfactor;
            image=data.data;%get;
            if isempty(image)
                dato=data;
                dato.data=[];
                return
            end
            
            
            if isempty(obj.imagebuffer)
                obj.imagebuffer=zeros(1,3,size(image,1),size(image,2),'single');
                obj.buffercounter=0;
                nx=1:size(image,1);ny=1:size(image,2);
                [obj.Xgrid,obj.Ygrid]=meshgrid(nx,ny);
                obj.mask=obj.getPar('loc_roimask'); %run at first time
            end
            obj.buffercounter=obj.buffercounter+1;
            bufferind=mod(obj.buffercounter-1,3)+1;
            obj.imagebuffer(1,bufferind,:,:)=image;
            if obj.buffercounter<3 %we need to wait for 3 images
                dato=[];
                return
            end
            
            indsort=[bufferind:3 1:bufferind-1];indsort=indsort([2 3 1]);
            fitimage=obj.imagebuffer(1,indsort,:,:);
%             fitimage=permute(fitimage,[1 2 4 3]);
            [xf, x_shape] = pseudo_col2row_major(fitimage);
            
            [outf, out_size] = obj.deepobj.forward(xf, x_shape);
            deepim = pseudo_row2_col_major(outf, out_size);
%             deepim=permute(deepim,[1 2 4 3]);
            probmap=squeeze(deepim(1,1,:,:));
            if all(size(obj.mask)==size(image)) %apply mask
                probmap=probmap.*obj.mask;
            end
            maxima=maximumfindcall((probmap)); 
            cutoff=p.pcutoff;
            maxind= (maxima(:,3)>cutoff);
            y=maxima(maxind,1);
            x=maxima(maxind,2);
            linind=sub2ind(size(image),x,y);
            Ymap=squeeze(deepim(1,3,:,:))*yfactor+obj.Ygrid;
            Xmap=squeeze(deepim(1,4,:,:))*xfactor+obj.Xgrid;
            photmap=squeeze(deepim(1,2,:,:));
            Zmap=squeeze(deepim(1,5,:,:));
            
            
            xfit=Xmap(linind);
            yfit=Ymap(linind);
            intensity=photmap(linind)*photfactor;
            zfit=Zmap(linind)*zfactor;
            pfit=probmap(linind);
            
            %export for SMAP
            v1=0*xfit+1;
            locs.frame=v1*data.frame;
            locs.xpix=xfit;
            locs.ypix=yfit;
            locs.znm=zfit;
            locs.phot=intensity;
            locs.prob=pfit;
            locs.PSFxpix=v1;
            locs.bg=median(image(:))*v1;
            locs.xpixerr=sqrt((locs.PSFxpix.*locs.PSFxpix+1/12*v1)./( locs.phot)+8*pi*(locs.PSFxpix.*locs.PSFxpix).^2.* locs.bg./( locs.phot).^2);
            locs.ypixerr=locs.xpixerr;
            dato=data;         
            dato.data=locs; 
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


pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='';
end