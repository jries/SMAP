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
        function out=run(obj,p)
            out=[];
            
            if isfield(obj.locData.loc,'clusterindex')
                field='clusterindex';
            elseif isfield(obj.locData.loc,'track_id')
                field='track_id';
            else
                disp('no clusterindex found');
            end
            clusterindex=obj.locData.loc.(field);
            posind=clusterindex>0;
            
            mi=max(clusterindex);
            n=rand(mi,1);
            [~,indsort]=sort(n);
            indsort=vertcat(1,indsort);
            clusterindexnew=zeros(size(clusterindex),'single');
            clusterindexnew(posind) = indsort(clusterindex(posind));
            obj.locData.loc.(field)=clusterindexnew;
            obj.locData.regroup;
        end
        
        function pardef=guidef(obj)
            pardef.plugininfo.name='ShuffleClusterIndices';
            pardef.plugininfo.description='ShuffelClusterIndex randomizes the index of clusters. This is useful for color-coded plotting of clusters.';
            pardef.plugininfo.type='ProcessorPlugin';
        end

    end
end

