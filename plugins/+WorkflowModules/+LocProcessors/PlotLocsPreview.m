classdef PlotLocsPreview<interfaces.WorkflowModule
%    Plots localizations of a single frame for optimization of fitting
%    parameters
    properties   
        outputfig
        previewdone
    end
    methods
       function obj=PlotLocsPreview(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1; 
       end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
        end
        function prerun(obj,p)
            obj.setPar('preview_locs',[]);
            obj.setPar('preview_image',[]);
            obj.setPar('preview_filtered',[]);
            obj.setPar('preview_background',[]);
            obj.setPar('preview_image_background',[]);
            obj.setPar('preview_peakfind',[]);
            obj.setPar('preview_mask',[]);
            obj.previewdone=false;
        end
        function pard=guidef(obj)
            pard.text.object=struct('Style','text','String','preview mode');
            pard.text.position=[1,1];
            pard.text.Width=1.0;
            pard.loc_previewmode.object=struct('Style','popupmenu','String','image|filtered|image-bg|bg','Value',2);
            pard.loc_previewmode.position=[1,1.7];
            pard.loc_previewmode.Width=.8;
            pard.loc_previewmode.TooltipString=sprintf('Determine which image to display in Preview mode. Peak finding is performed on norm(image)');

            pard.plugininfo.type='WorkflowModule'; 
            pard.plugininfo.description='Plots localizations of a single frame for optimization of fitting parameters';
        end
        function output=run(obj,data,p)
            output=[];
            if obj.getPar('loc_preview') %&& ~obj.previewdone
                locs=data.data;
                previewlocs=obj.getPar('preview_locs');
                if isempty(locs)
                    locs=previewlocs;
                elseif ~isempty(previewlocs) 
                    fn=fieldnames(previewlocs);
                    for k=1:length(fn)
                        locs.(fn{k})(end+1:end+length(previewlocs.(fn{k})))=previewlocs.(fn{k});
                    end
                end
                obj.setPar('preview_locs',locs);
                
                %make figure
                if isempty(obj.outputfig)|| ~isvalid(obj.outputfig)
                    obj.outputfig=figure;
                end
                fig=figure(obj.outputfig);
                obj.setPar('loc_outputfig',fig);
                ax=gca;
                hold(ax,'off')
                
                %plot image
                titletxt=p.loc_previewmode.selection;
                switch p.loc_previewmode.selection
                    case 'image-bg'
                        implot=obj.getPar('preview_image_background');
                    case 'image'
                        implot=obj.getPar('preview_image');
                    case 'filtered'
                        implot=obj.getPar('preview_filtered');
                    case 'bg'
                        implot=obj.getPar('preview_bg');
                    otherwise
                        implot=[];
                end
                if isempty(implot)
                    disp([p.loc_previewmode.selection ' not available. Use raw image instead for preview']);
                    implot=obj.getPar('preview_image');
                    titletxt=['image, as ' titletxt ' not available'];
                end
                
                if isempty(implot)
                    obj.setPar('status','image could not be loaded for preview')
                     error ('image could not be loaded for preview')
%                     return
                end
                imagesc(ax,implot);
                title(ax,titletxt);
                colorbar(ax);
%                 colormap(ax,'jet');
                hold(ax,'on')
                axis(ax,'equal')
                %mask
               maskim=obj.getPar('loc_roimask');
               if ~isempty(maskim)
                   maxv=ax.CLim(2);
                   imagesc(ax,'CData',maskim*0+maxv,'AlphaData',double(~maskim)*0.4) ;
               end

                
                % plot ROIs
                maxima=obj.getPar('preview_peakfind');
                if isempty(maxima)||isempty(maxima.xpix)
                     obj.setPar('status','no localizations found')
                     obj.setPar('errorindicator','no localizations found')
                     error ('no localizations found')
                end
                col=[0.3 0.3 0.];
                dn=floor(obj.getPar('loc_ROIsize')/2);
                for k=1:length(maxima.xpix)
                    pos=[maxima.xpix(k)-dn maxima.ypix(k)-dn maxima.xpix(k)+dn maxima.ypix(k)+dn ];
                    plotrect(ax,pos,col);
                end
                
                %plot locs
                if ~isempty(locs)&&~isempty(locs.xpix)
                    
                    plot(ax,locs.xpix,locs.ypix,'k.','Parent',ax,'MarkerSize',4)
                    dn=ceil((obj.getPar('loc_ROIsize')-1)/2);
                    if isempty(dn)
                        dn=3;
                    end
                    for k=1:length(locs.xpix)
                        pos=[locs.xpix(k)-dn locs.ypix(k)-dn locs.xpix(k)+dn locs.ypix(k)+dn ];
                        col=[1 0. 1];
                        plotrect(ax,pos,col);
                    end
      
                end
                obj.previewdone=true; %avoid calling it again with eof and overwriting stuff
                
            end 


        end

    end
end
