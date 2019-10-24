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
        scaling
        previewfigure
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
            addpath([deeppath  '/build'])
            
            
            obj.imagebuffer=[];
            if obj.getPar('loc_preview')
                obj.buffersize=1;
            else
                obj.buffersize=p.buffersize;
            end
            obj.preview=obj.getPar('loc_preview');
            ind=strfind(p.modelfile,'_');
%             [mpath mfile]=fileparts(p.modelfile);
            jsonfile=dir([p.modelfile(1:ind(2)) '*.json']);
%             jsonfile=[p.modelfile(1:ind(end-2)) 'param.json'];
            jsontxt=fileread([jsonfile.folder filesep jsonfile.name]);
            param=jsondecode(jsontxt);
            obj.scaling=param.Scaling;
            obj.scaling.channels=param.HyperParameter.channels_in;
            obj.deepobj=mex_interface(str2fun([deeppath '/build/deepsmlm_interface']), p.modelfile,obj.scaling.channels,'cpu',10); %waht is 10 ? XXXXXXXXXXX
        end
        function dato=run(obj,data,p)
            %
            timeblocks=obj.scaling.channels; %3 frames analyzed together
            dato=[];
            bufferfactor=obj.scaling.linearisation_buffer;
            xfactor=obj.scaling.dy_max*bufferfactor;
            yfactor=obj.scaling.dx_max*bufferfactor;
            zfactor=obj.scaling.z_max*bufferfactor;
            photfactor=obj.scaling.phot_max*bufferfactor;
            bgfactor=obj.scaling.bg_max*bufferfactor;
            
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
                firstframe=true;
            else
                firstframe=false;
            end
            obj.buffercounter=obj.buffercounter+1;
            obj.imagebuffer(obj.buffercounter,:,:)=image;
            obj.bufferinfo.frame(obj.buffercounter)=data.frame;
            
            lastframe=false;
            if obj.buffercounter<obj.buffersize %we need to wait for buffersize images. check if eof, then process the rest.
                if ~data.eof
                    return
                else
                    lastframe=true;
                end
            end
            
            %first block: only remove last frame in block
            %last block: only remove first frame in block
            %otherwise remove first and last
            %always create overlap of two frames
            

            imbufferadu=obj.imagebuffer; %already in ADU
            [xf, x_shape] = pseudo_col2row_major(imbufferadu);
            
            [outf, out_size] = obj.deepobj.forward(xf, x_shape);
            deepim = pseudo_row2_col_major(outf, out_size);
             deepim(:,7,:,:)=imbufferadu;
             
            
%             deepim=permute(deepim,[1 2 4 3]);
            if timeblocks ==1 
                if lastframe
                    frange=1:obj.buffercounter;
                else 
                    frange=1:obj.buffersize;
                end
            else %later extend to any number
                if obj.buffersize==1
                    frange=1;
                elseif firstframe
                    frange=1:obj.buffersize-1;
                elseif lastframe
                    frange=2:obj.buffercounter;
                else
                    frange=2:obj.buffersize-1;
                end
            end
                
            for k=frange %process each frame individually
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
                bgmap=squeeze(deepim(k,6,:,:))*bgfactor;

                xfit=Xmap(linind);
                yfit=Ymap(linind);
                intensity=photmap(linind)*photfactor;
                zfit=Zmap(linind)*zfactor;
                pfit=probmap(linind);
                dx=dxmap(linind);
                dy=dymap(linind);
                bg=bgmap(linind);

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
                locs.xpix=xfit;
                locs.ypix=yfit;
                locs.phot=intensity;
                locs.znm=zfit;
                locs.bg=bg;
                
                dato=data; 
                dato.frame=obj.bufferinfo.frame(k);
                dato.ID=dato.frame;
                dato.data=locs; 
                obj.output(dato);
                dato=[]; %not to output it again when returning from the function
            end
            if obj.buffersize>1 && timeblocks>1
            obj.imagebuffer(1:2,:,:)=obj.imagebuffer(obj.buffersize-1:end,:,:);
            obj.bufferinfo.frame(1:2)=obj.bufferinfo.frame(obj.buffersize-1:end,:,:);
            end
            if obj.preview
                obj.setPar('preview_peakfind',locs);
                if isempty(obj.previewfigure) || ~isvalid(obj.previewfigure)
                    obj.previewfigure=figure;
                end
                deepimplot=permute(deepim,[4 3 2 1]);
                deepimplot(:,:,2,:)=deepimplot(:,:,2,:)*photfactor;
                deepimplot(:,:,3,:)=deepimplot(:,:,3,:)*xfactor;
                deepimplot(:,:,4,:)=deepimplot(:,:,4,:)*yfactor;
                deepimplot(:,:,5,:)=deepimplot(:,:,5,:)*zfactor;
                deepimplot(:,:,6,:)=deepimplot(:,:,6,:)*bgfactor;
                tags={[],[],{'probability','photons','dx','dy','znm','bg','image'},[]};
                 imx(deepimplot,'Parent',obj.previewfigure,'Tags',tags)
            end
            obj.buffercounter=timeblocks-1; %before timeblocks: did we lose frames? XXXX
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

% pard.context_framest.object=struct('Style','text','String','Channels');
% pard.context_framest.position=[5,1];
% pard.context_framest.Width=1;
% pard.context_frames.object=struct('Style','edit','String','3');
% pard.context_frames.position=[5,2];
% pard.context_frames.Width=.35;
% pard.context_frames.TooltipString='Number of frames ananlyzed simultaneously';
% pard.context_framest.TooltipString=pard.context_frames.TooltipString;

pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='';
end