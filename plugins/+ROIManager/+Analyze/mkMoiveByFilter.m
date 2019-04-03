classdef mkMoiveByFilter<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=mkMoiveByFilter(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)
%             p.sourceRange=[1 2430];
%             p.stepSize = 5;
%             p.windowSize = 30;
%             p.frameRate = 12;
            
            %% init
            se = obj.locData.SE;
            
            lowB = p.sourceRange(1); % the lower boundary
            framesOfMovie = zeros([obj.getPar('se_sitefov'),obj.getPar('se_sitefov'),3, floor(range(p.sourceRange)/p.stepSize)]); % the stack of frames. n*m*3*numOfFrames, where n*m is the dim of a ROI
            
            %% sliding window
            for k = 1:floor(range(p.sourceRange)/p.stepSize)
                nLayer = size(obj.locData.layer,2);
                for l = 1:nLayer
                    obj.setPar('selectedField', {'order',lowB, lowB+p.windowSize, 1,1 }, 'layer', l);
                end
                se.currentsite.image = [];
                framesOfMovie(:,:,:,k)= se.plotsite(se.currentsite, se.processors.preview.guihandles.siteax, se.processors.preview.guihandles.cellax).image;
                lowB = lowB+p.stepSize;
            end
            
            %% export the vedio
            if ~p.saveFrames%% preview
                m = immovie(framesOfMovie);
                if p.preview
                    implay(m,p.frameRate) % will initiate a pop-up asking the user saving or not
                end
                myVideo = VideoWriter([p.folderPath '\' p.videoPath]);
                myVideo.FrameRate = p.frameRate;  % Default 30
                myVideo.Quality = 100;    % Default 75
                open(myVideo);
                writeVideo(myVideo, m);
                close(myVideo);
            else
                mkdir([p.folderPath '\' p.videoPath]);
                for k = 1:size(framesOfMovie,4)
                    imwrite(framesOfMovie(:,:,:,k), [p.folderPath '\' p.videoPath '\' p.videoPath '_' num2str(k) '.tif']);
                end
            end
            out=[];
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end




function pard=guidef(obj)



pard.t1.object=struct('String','ROI number (from to)','Style','text');
pard.t1.position=[1,1];
pard.t1.Width=1;

pard.sourceRange.object=struct('String','1 100','Style','edit');
pard.sourceRange.position=[1,2];
pard.sourceRange.Width=1;

pard.preview.object=struct('String','Preview','Style','checkbox', 'Value', 1);
pard.preview.position=[1,3];
pard.preview.Width=1;

pard.saveFrames.object=struct('String','Save frames','Style','checkbox', 'Value', 0);
pard.saveFrames.position=[1,4];
pard.saveFrames.Width=1;

pard.t2.object=struct('String','Step size','Style','text');
pard.t2.position=[2,1];
pard.t2.Width=1;

pard.stepSize.object=struct('String','5','Style','edit');
pard.stepSize.position=[2,2];
pard.stepSize.Width=1;

pard.t3.object=struct('String','Window size','Style','text');
pard.t3.position=[3,1];
pard.t3.Width=1;

pard.windowSize.object=struct('String','30','Style','edit');
pard.windowSize.position=[3,2];
pard.windowSize.Width=1;

pard.t4.object=struct('String','Frame rate','Style','text');
pard.t4.position=[3,3];
pard.t4.Width=1;

pard.frameRate.object=struct('String','12','Style','edit');
pard.frameRate.position=[3,4];
pard.frameRate.Width=1;

pard.folder.object=struct('String','Choose folder','Style','pushbutton','Callback',{{@saveTo_callback,obj}});
pard.folder.position=[7,3];
pard.folder.Width=1;     

pard.plugininfo.type='ROI_Analyze';

pard.t_folderPath.object=struct('String','Save to','Style','text');
pard.t_folderPath.position=[7,1];
pard.t_folderPath.Width=1;

pard.folderPath.object=struct('String','/','Style','edit');
pard.folderPath.position=[7,2];
pard.folderPath.Width=1;

pard.videoPath.object=struct('String', 'my_video.avi', 'Style','edit');
pard.videoPath.position=[8,2:3];
pard.videoPath.Width=2;    
end

function saveTo_callback(a,b,obj)
f=obj.getSingleGuiParameter('folderPath');
f=uigetdir(f,'Choose folder for saving the result');
if ~f
    return
end
obj.setGuiParameters(struct('folderPath',f));
setvisibility(obj)
end
