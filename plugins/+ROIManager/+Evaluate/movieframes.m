classdef movieframes<interfaces.SEEvaluationProcessor
    properties
        outmovie
    end
    methods
        function obj=movieframes(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj, p)
            %add save button (as tiff or movie)
            %average: optional
            out=[];
            layers=find(p.sr_layerson);
            ax=obj.setoutput(p.moviegallery.selection);
            locs=obj.getLocs({'xnmrot','ynmrot','frame'},'layer',layers,'size',p.se_siteroi*[1 1]);  
            maxf=max(locs.frame);
            numt=ceil(maxf/p.df);
            rows=ceil(numt/p.columns);
            nedge=-p.se_siteroi/2:p.se_sitepixelsize:p.se_siteroi/2;
            spix=length(nedge)-1;
            outim=zeros(spix*rows,spix*p.columns,3);
            allim=zeros(spix,spix,numt);
            h=fspecial('gauss',ceil(p.gausssize*6),p.gausssize);
%             maxv=max(h(:));
            ic=1;ir=1;
%             allout=zeros(spix,spix);
            for k=1:numt
                indh=locs.frame>(k-1)*p.df+1 & locs.frame<(k-1)*p.df+p.window;
                hc=histcounts2(locs.xnmrot(indh),locs.ynmrot(indh),nedge,nedge)';
                hf=filter2(h,hc);
                allim(:,:,k)=hf;
            end
            meanim=sum(allim,3);
            meanim=meanim/quantile(meanim(:),.9999);
            maxv=quantile(allim(:),0.9999);
            allim(allim>maxv)=maxv;
            allim=allim/maxv;
            allim(end,:,:)=1;
            allim(:,end,:)=1;

            cmap=hot(256); 
            
            scalebarlength=round(100/p.se_sitepixelsize);
            if isempty(p.exposuretime)
                dt=obj.locData.files.file(obj.site.info.filenumber).info.timediff;
            else
                dt=p.exposuretime;
            end
%             F(numt)=im2frame(ind2rgb(uint8(allim(:,:,1)*255),cmap));
            ismovie=p.moviegallery.Value==2;
            for k=1:numt
                imh=ind2rgb(uint8(allim(:,:,k)*255),cmap);
                imh=imh+(meanim)*0.7;
                imh(imh>1)=1;
                
                %time
                tf=(k-1)*dt/1000*p.df;
                t=[num2str(tf,'%2.1f') 's'];
                tim=text2im(t);
                s=size(tim);
                imh(1:s(1),1:s(2),:)=repmat(tim,1,1,3);
                
                
                %scalebar
                imh(end-6:end-4,end-scalebarlength-10:end-10-1,:)=1;
                if ismovie
                    F(k)=im2frame(imh);
                else
                    outim((ir-1)*spix+1:ir*spix,(ic-1)*spix+1:ic*spix,:)=imh;
                    ic=ic+1;
                    if ic>p.columns
                        ic=1;
                        ir=ir+1;
                    end
                end
                
            end
            if ismovie
                    F(numt)=im2frame(imh*0+0.5);
                    axis(ax,'off')
                    axis(ax,'equal')
                    movie(ax,F,2,6)  
                    outmovie=F;
            else

                    image(ax,outim)
                    colormap(ax,hot)
                    axis(ax,'equal')
                    axis(ax,'off')
                    outmovie=outim;


                    
            end
            
        end
     
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end

end



function pard=guidef(obj)

pard.moviegallery.object=struct('String',{{'Gallery','Movie'}},'Style','popupmenu');
pard.moviegallery.position=[1,1];
pard.moviegallery.Width=2;

pard.windowt.object=struct('String','Window (frames)','Style','text');
pard.windowt.position=[2,1];
pard.windowt.Width=2;
pard.window.object=struct('String','100','Style','edit');
pard.window.position=[2,3];
pard.window.Width=0.5;

pard.dft.object=struct('String','Stepsize (frames)','Style','text');
pard.dft.position=[3,1];
pard.dft.Width=2;
pard.df.object=struct('String','100','Style','edit');
pard.df.position=[3,3];
pard.df.Width=0.5;


pard.columnst.object=struct('String','Columns','Style','text');
pard.columnst.position=[4,1];
pard.columnst.Width=2;
pard.columns.object=struct('String','5','Style','edit');
pard.columns.position=[4,3];
pard.columns.Width=0.5;

pard.gausssizet.object=struct('String','Filter (pix)','Style','text');
pard.gausssizet.position=[5,1];
pard.gausssizet.Width=2;
pard.gausssize.object=struct('String','3','Style','edit');
pard.gausssize.position=[5,3];
pard.gausssize.Width=0.5;

pard.exposuretimet.object=struct('String','Exposure (ms)','Style','text');
pard.exposuretimet.position=[6,1];
pard.exposuretimet.Width=2;
pard.exposuretime.object=struct('String','','Style','edit');
pard.exposuretime.position=[6,3];
pard.exposuretime.Width=0.5;

% pard.dxt.Width=3;
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','se_sitepixelsize'};
pard.plugininfo.type='ROI_Evaluate';
pard.plugininfo.description='';
end

