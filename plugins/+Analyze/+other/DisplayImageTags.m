classdef DisplayImageTags<interfaces.DialogProcessor&interfaces.SEProcessor
%    Filters a localization attribute. This is mainly used for the BatchAnalysis
%     plugin.
    methods
        function obj=DisplayImageTags(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
                obj.showresults=true;
        end
        
        function out=run(obj,p)  
            out=[];
            
            files=obj.locData.files.file;
            for f=1:length(files)
                imagetags=files(f).imagetags;
                for k=1:length(imagetags.tags)
                    name=imagetags.tags{k}(1:min(length(imagetags.tags{k}),7));
                    ax=obj.initaxis([name ':' num2str(f) num2str(k)]);
%                     tab=uitab(tg,'Title',imagetags.tags{k});
%                     ax=axes('Parent',tab);
                    dat=imagetags.data(k,:);
%                     if ischar(dat{1})
%                         datm=str2double(dat);
%                     else
%                         datm=cell2mat(dat);
%                     end
                    frames=(1:length(dat))';
                    frameind=dat~=0;
                    plot(ax,frames(frameind), dat(frameind))
                    xlabel('frame')
                    ylabel(imagetags.tags{k})
                    [~, fn]=fileparts(files(f).name);
                    title(fn,'Interpreter','none')
                    xlim([min(frames(frameind)) max(frames(frameind))])
                end
            end
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end




function pard=guidef(obj)

pard.t1.object=struct('String','Display image tags vs frame.','Style','text');
pard.t1.position=[1,1];
pard.t1.Width=4;

pard.plugininfo.type='ROI_Analyze';
pard.plugininfo.description='Display image tags vs. frame.';
end