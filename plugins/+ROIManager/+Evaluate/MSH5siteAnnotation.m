classdef MSH5siteAnnotation<interfaces.SEEvaluationProcessor
%     Calculates the line-profile along user-defined direction and fits it
%     with a Gaussian or double-Gaussian model
    properties
        %roicoordinates
        roihandles={[]};
        roisize
        roisizey
%         currentroi=0;
        %axis
        %images
%         roihandle
    end
    methods
        function obj=MSH5siteAnnotation(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj, p)
            out=[];
            display=obj.display;

            out.meioticStage=obj.getSingleGuiParameter('meioticStage').selection;
            out.rowstoDiplotene=str2double(obj.getSingleGuiParameter('rowstoDiplotene'));
            out.comment=obj.getSingleGuiParameter('comment');         
            markerlayer=1;
            clusterMarker=1;
            layerson = find(obj.locData.getPar('sr_layerson'));             % check the layers used
            visitFlag = false;                                      % a flag for the conditional loop for the first round
            fieldQ = {'locprecnm','locprecznm','PSFxnm','xnm','znm','ynm','frame','xnmrot','ynmrot','phot1','phot2','filenumber','channel'};    % fields will be used %%%%%%%%%%%%%%%%%%% need to add more for displayer
            
            %layerson=layerson(layerson~=cosa1layer);
            %Pick up data points
            fieldQNLayer = [fieldQ 'layer'];
            for k = layerson                                        % go through layers, collect all filtered locs -----from Yu-Le's LocMoFitGUI.m
                if ~visitFlag
                    locs=obj.getLocs(fieldQ,'layer',k,'size','freeroi');
                    locs.layer = ones(size(locs.(fieldQ{1})))*k;
                    fieldExact = fieldnames(locs);
                    lRm = ~ismember(fieldExact, fieldQNLayer);
                    locs = rmfield(locs,fieldExact(lRm));
                    visitFlag=true;
                else
                    locsNewLayer = obj.getLocs(fieldQ,'layer',k,'size','freeroi');
                    for l = length(fieldQNLayer):-1:1
                        if l == length(fieldQNLayer)
                            locsNewLayer.layer = ones(size(locsNewLayer.(fieldQNLayer{l-1})))*k;
                        end
                        if isfield(locs, fieldQNLayer{l})
                            locs.(fieldQNLayer{l})=[locs.(fieldQNLayer{l});locsNewLayer.(fieldQNLayer{l})];
                        end
                    end
                end
            end

            %Cluster marker layer Points
            warning('off', 'stats:pdist2:DataConversion');
            if clusterMarker==1
                epsilon=20;
                minpts=2;
                markerpoints=[locs.xnm(locs.layer==markerlayer) locs.ynm(locs.layer==markerlayer) locs.znm(locs.layer==markerlayer)];
%                 D = pdist2(markerpoints,markerpoints);
                idx=dbscan(markerpoints,40,50);
                [GC,GR] = groupcounts(idx);
                out.meanmarkerpoints=zeros(size(unique(idx),1),3);
                c=1;
                for cluster=unique(idx).'
                    out.meanmarkerpoints(c,:)=mean(markerpoints(idx==cluster,:),1);
                    c=c+1;
                end
                out.GC=GC;
                out.GR=GR;
                out.nclusters=max(GR);
                out.nlocs=size(markerpoints,2);
                out.dbscanparams=[epsilon minpts];
                if display
                    figure('Name','3dpreview - dbscan cluster of marker locs');
                    hold on
                    scatter3(markerpoints(:,1), markerpoints(:,2), markerpoints(:,3),15,idx,'filled'); % 
%                     scatter3(locs.xnm(locs.layer~=markerlayer), locs.ynm(locs.layer~=markerlayer), locs.znm(locs.layer~=markerlayer),5,'red','filled'); % 
                    
                    colormap("lines")
                    daspect([1 1 1])
                    hold off
                end
            end
            warning('on', 'stats:pdist2:DataConversion');
            
        end
     
        function pard=guidef(obj)
            pard=guidef(obj);
        end
%         
    end

end


function pard=guidef(obj)

pard.commentT.object=struct('String','Comment:','Style','text');
pard.commentT.position=[3.5,1];
pard.commentT.Width=2;

pard.comment.object=struct('Style','edit','String','NEWVAL_dye=AF647_refDye=CF680_protein=MSH-5::V5_ref=HTP-3::HA_orientation=frontal_shapedescriptor=irregular_use=True_imagedSecond=False_additionalcomment=NA');
pard.comment.position=[5.5,1];
pard.comment.Width=4;
pard.comment.Height=2;
pard.comment.TooltipString='Input Comment. If you need to add more values to outputtable add them in this format: NEWVAL_name1=value1_name2=value2';

pard.rowstoDiploteneT.object=struct('String','Rows to diplotene:','Style','text');
pard.rowstoDiploteneT.position=[2.5,1];
pard.rowstoDiploteneT.Width=3;
pard.rowstoDiploteneT.TooltipString='';

pard.rowstoDiplotene.object=struct('String','7','Style','edit');
pard.rowstoDiplotene.position=[2.5,4.5];
pard.rowstoDiplotene.Width=0.5;
pard.rowstoDiplotene.TooltipString='Site is this many rows away from diplotene:';

pard.meioticStageT.object=struct('String','Meiotic Stage:','Style','text');
pard.meioticStageT.position=[1,1];
pard.meioticStageT.Width=2;

pard.meioticStage.object=struct('Style','popupmenu','String',{["midP","lateP","earlyP","diplotene","premeiotic","TZ","diakinesis"]});
pard.meioticStage.position=[1,2.5];
pard.meioticStage.Width=1;
pard.meioticStage.TooltipString='Choose stage that best describes entire image.';



% pard.equaldistance.object=struct('Style','checkbox','String','equal distances');
% pard.equaldistance.position=[2,1];
% pard.equaldistance.Width=2;

% pard.fitpeaks.object=struct('Style','pushbutton','String','fit peaks','Callback',@obj.fit_callback);
% pard.fitpeaks.position=[1,1];
% pard.fitpeaks.Width=2;



% pard.dxt.Width=3;
% pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','se_sitepixelsize'};
pard.plugininfo.type='ROI_Evaluate';
pard.plugininfo.description='';
end
