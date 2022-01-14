classdef clusterAnalysis<interfaces.SEEvaluationProcessor
%     performs clustering on data using dbscan, k-means, k-medion, Gaussian
%     mixture models or hierarchical clustering and returns positions of
%     clusters and statistics
    properties
        boundary
    end
    methods
        function obj=clusterAnalysis(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj, inp)
            % Import the kimograph
%             site = obj(1).locData
            layers=find(inp.sr_layerson);
            locs=obj.getLocs({'locprecnm','locprecznm','xnm','ynm','frame','znm','layer'},'layer',layers,'size',inp.se_siteroi/2);  
%             if isempty(locs.znm)
%                 locs.znm=0*locs.xnm;
%                 locs.locprecznm=0*locs.xnm+1;
%             end
            X=horzcat(locs.xnm,locs.ynm,locs.znm);
            Xerr=horzcat(locs.locprecnm,locs.locprecnm,locs.locprecznm);
            C=[];
            switch inp.method.selection
                case 'dbscan'
                    idx = dbscan(X,inp.eps_dbscan,inp.k_dbscan);
                    max(idx)
                case 'k-means'
                    [idx,C] = kmeans(X,inp.N_kmeans,'Replicates',10);
                  
                case 'k-medoid'
                    [idx,C] = kmedoids(X,inp.N_kmeans,'Replicates',10);
                case 'Gaussian mixture'
                    gm = fitgmdist(X,inp.N_kmeans,'Replicates',10);
                    idx = cluster(gm,X);
                    P = posterior(gm,X);
                    C=gm.mu;
                case 'hierarchical'
        
                    idx = clusterdata(X,inp.N_kmeans);
            end
            %evaluate positions
            Cdat=makeclusterpos(X,idx,inp.posmethod.selection,Xerr);
            
            %plot
            ax=obj.setoutput('scatter');
            delete(ax.Children)
            hold(ax,'off')
%             idxn=idx-min(idx)+1;
            
            nc=max(idx);
            hc=histcounts(idx(idx>0));
            [~,idr]=sort(hc,'descend');
            
            if size(idr,2)>1
                
                ddC.xyz=norm(Cdat(idr(2),:)-Cdat(idr(1),:));
                ddC.xy=norm(Cdat(idr(2),1:2)-Cdat(idr(1),1:2));
                if size(Cdat,2)>2
                ddC.z=abs(Cdat(idr(2),3)-Cdat(idr(1),2));
                end
            else
                ddC.xyz=NaN;
                ddC.xy=NaN;
                ddC.z=NaN;
            end
%             [~,idr]=sort(rand(nc));
            map=hsv(nc+1);
            idxn=idx;
            idxn(idx>0)=idr(idx(idx>0));
            incluster=idxn>0;
            if isempty(locs.znm)
                scatter(ax,X(incluster,1),X(incluster,2),5,map(idxn(incluster),:))
                hold (ax,'on')
                 scatter(ax,X(~incluster,1),X(~incluster,2),3,[0.5 0.5 0.5])
                 if size(idr,2)>1
                  plot(ax,Cdat(idr(1:2),1),Cdat(idr(1:2),2),'k-o')
                 end
            else
                scatter3(ax,X(incluster,1),X(incluster,2),X(incluster,3),5,map(idxn(incluster),:))
                hold (ax,'on')
                scatter3(ax,X(~incluster,1),X(~incluster,2),X(~incluster,3),3,[0.5 0.5 0.5])
                if size(idr,2)>1
                 plot3(ax,Cdat(idr(1:2),1),Cdat(idr(1:2),2),Cdat(idr(1:2),3),'k-o')
                end
            end
            hold (ax,'on')
%                 scatter3(ax,Cdat(idr(1),1),Cdat(idr(1),2),Cdat(idr(1),3),20,'k')
%                 scatter3(ax,Cdat(idr(2),1),Cdat(idr(2),2),Cdat(idr(2),3),20,'k')
               
            if ~isempty(C)
%                 scatter3(ax,C(idr(1),1),C(idr(1),2),C(idr(1),3),20,'b+')
%                 scatter3(ax,C(idr(2),1),C(idr(2),2),C(idr(2),3),20,'b+')
                if ~isempty(locs.znm)
                    plot3(ax,C(idr(1:2),1),C(idr(1:2),2),C(idr(1:2),3),'b-+')
                    dC.xyz=norm(C(idr(2),:)-C(idr(1),:));
                    dC.xy=norm(C(idr(2),1:2)-C(idr(1),1:2));
                    dC.z=abs(C(idr(2),3)-C(idr(1),3));
                else
                    dC=[];
                end
            else
                dC=[];
            end
            title(ax,['d: ' num2str(ddC.xyz,3)]);
%             f=figure;
%             h=gscatter(X(:,1),X(:,2),idx);
            
            out.method=inp.method.selection;
            out.clustermodel=C;
            out.clusterdata=Cdat;
            out.dmodel=dC;
            out.ddat=ddC;
            out.posmethod=inp.posmethod;
        end
     
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end

end


function C=makeclusterpos(X,idx,method,Xerr)
wm=false;
switch method
    case 'mean'
        fun=@mean;
    case 'robust Mean'
        fun=@robustMean;
    case 'median'
        fun=@median;
    case 'weighted Mean'
        wm=true;
end


C=zeros(max(idx),size(X,2));
for k=1:max(idx)
    indh=idx==k;
    if isempty(indh)
        continue
    end
    if wm
        C(k,:)=sum(X(indh,:)./Xerr(indh,:),1)./sum(1./Xerr(indh,:),1);
    else
        C(k,:)=fun(X(indh,:));
    end
end
end

function pard=guidef(obj)


p(1).value=1; p(1).on={'k_dbscant','k_dbscan','seteps_dbscan','eps_dbscan'}; p(1).off={'N_kmeans','N_kmeanst'};
p(2).value=2; p(2).on=p(1).off; p(2).off=p(1).on;
p(3)=p(2);p(3).value=3;p(4)=p(2);p(4).value=4;p(5)=p(2);p(5).value=5;

pard.method.object=struct('Style','popupmenu','String',{{'dbscan','k-means','k-medoid','Gaussian mixture','hierarchical'}},'Callback',{{@obj.switchvisible,p}});
pard.method.position=[1,1];
pard.method.Width=2;

pard.k_dbscant.object=struct('String','min objects in neighoburhood','Style','text');
pard.k_dbscant.position=[2,1];
pard.k_dbscant.Width=3;

pard.k_dbscan.object=struct('String','10','Style','edit');
pard.k_dbscan.position=[2,4];

pard.seteps_dbscan.object=struct('String','eps (neighbourhood radius): ','Style','text');
pard.seteps_dbscan.position=[3,1];
pard.seteps_dbscan.Width=3;

pard.eps_dbscan.object=struct('String','5','Style','edit');
pard.eps_dbscan.position=[3,4];

pard.N_kmeanst.object=struct('String','Number of clusters','Style','text');
pard.N_kmeanst.position=[2,1];
pard.N_kmeanst.Width=3;

pard.N_kmeans.object=struct('String','2','Style','edit');
pard.N_kmeans.position=[2,4];

pard.posmethodt.object=struct('String','position estimation','Style','text');
pard.posmethodt.position=[4,1];
pard.posmethodt.Width=2;
pard.posmethod.object=struct('String',{{'mean','robust Mean','median','weighted Mean'}},'Style','popupmenu');
pard.posmethod.position=[4,3];
pard.posmethod.Width=2;
% pard.dxt.Width=3;
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi'};
pard.plugininfo.type='ROI_Evaluate';
pard.plugininfo.description=' performs clustering on data using dbscan, k-means, k-medion, Gaussian mixture models or hierarchical clustering and returns positions of clusters and statistics ';
end
