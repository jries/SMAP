classdef ShuffelClusterIndex<interfaces.DialogProcessor
    % ShuffelClusterIndex randomizes the index of clusters. This is useful
    % for color-coded plotting of clusters.
    methods
        function obj=ShuffelClusterIndex(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters={'sr_pixrec','numberOfLayers','sr_pos','sr_size','layers','sr_layerson'};
            obj.history=true;
            obj.showresults=false;
        end

        function initGuiFinal(obj)
            locf=obj.getPar('locFields').String;
            indg=find(strcmp(locf,'groupindex'));
            if ~isempty(indg)
                obj.guihandles.assignfield1.Value=indg;
                obj.selectfield_callback;
            end
            
        end
        function out=run(obj,p)
            out=[];
            field=obj.getSingleGuiParameter('assignfield1').selection;
            outfield= obj.getSingleGuiParameter('resultfieldh');

            clusterindex=obj.locData.loc.(field);
            posind=clusterindex>0;
            
            mi=max(clusterindex);
            n=rand(mi,1);
            [~,indsort]=sort(n);
%             indsort=vertcat(1,indsort);
            clusterindexnew=zeros(size(clusterindex),'single');
            clusterindexnew(posind) = indsort(clusterindex(posind));
            obj.locData.setloc(outfield,clusterindexnew);
            % obj.locData.loc.(outfield)=clusterindexnew;
            obj.locData.regroup;
        end
        function selectfield_callback(obj,a,b)
            field=obj.getSingleGuiParameter('assignfield1').selection;
            obj.setGuiParameters(struct('resultfieldh',[field '_rnd']));
        end
        
        function pardef=guidef(obj)
            pardef.plugininfo.name='ShuffleClusterIndices';
            pardef.plugininfo.description='ShuffelClusterIndex randomizes the index of clusters. This is useful for color-coded plotting of clusters.';
            pardef.plugininfo.type='ProcessorPlugin';

            pardef.assignfield1.object=struct('Style','popupmenu','String','n1|n2','Callback',@obj.selectfield_callback);
            pardef.assignfield1.position=[1,3.5];

            pardef.resultfieldh.object=struct('String',{{''}},'Style','edit');
            pardef.resultfieldh.position=[1,1];
            pardef.resultfieldh.Width=1.5;

            pardef.txt.object=struct('String',' = random shuffle of ','Style','text');
            pardef.txt.position=[1,2.5];


            pardef.syncParameters={{'locFields','assignfield1',{'String'}}};

        end

    end
end

