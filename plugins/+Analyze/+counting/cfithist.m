classdef cfithist<interfaces.DialogProcessor
    %  Copyright (c)2017 Ries Lab, European Molecular Biology Laboratory,
    %  Heidelberg. This file is part of Single Molecule Analysis Platform (SMAP).
    
    % fits histogram of localizaitons/cluster with an analytical model to
    % determine number of proteins, maturation efficiency,  number of
    % re-activations
    
    methods
        function obj=cfithist(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
        end
        
        function out=run(obj,p)
            
            if p.bootstrap
                numberOfSubsets = length(obj.getResults('counting_histogram'));
            else
                numberOfSubsets = 1;
            end
            meanlocs = [];
            ticked = [];
            collection = [];
            if p.N0_fit
                ticked = [ticked string('N0_v')];
            end
            if p.pblink_fit
                ticked = [ticked string('pblink_v')];
            end
            if p.pmature_fit
                ticked = [ticked string('pmature_v')];
            end
            if p.monomer_fit
                ticked = [ticked string('monomer_v')];
            end
            if p.ncluster_fit
                ticked = [ticked string('ncluster_v')];
            end
            nTicked = length(ticked);
            
            for i = 1:numberOfSubsets
                if p.bootstrap
                    allHistogram=obj.getResults('counting_histogram');
                    histogram=allHistogram{i};
                else
                    histogram=obj.getResults('counting_histogram');
                end
%               histogram=obj.locData.guiData.counting.histogram;
                pout_ref=cluster_mmaple_fithist(p,histogram);pout=p_ref;
                rep = 0;
                while 0
                    pout=cluster_mmaple_fithist(pout_ref,histogram);
                    stopSignal = round(pout.N0_v ,3) == round(pout_ref.N0_v,3) & round(pout.pmature_v,3) == round(pout_ref.pmature_v,3) & round(pout.pblink_v,3) == round(pout_ref.pblink_v,3);
                    if stopSignal||(rep>10)
                        break
                    else
                        pout_ref=pout;
                        rep=rep+1;
                    end
                end
                
                meanlocs(i) = str2double(pout.meanlocs);
                for ii = 1:nTicked
                    tickedFlag = ticked(ii);
                    collection(i,ii) = pout_ref.(char(tickedFlag));
                end
%               locs=obj.locData.getloc({'frame','xnm','ynm','phot','bg','PSFxnm','locprecnm'},'layer',1,'position','roi');

%               pout.par=cluster_counting(locs,p);
            end
            
            
            if p.bootstrap
                outCollection = [mean(meanlocs) mean(collection);std(meanlocs) std(collection)];
                sOutCollection = size(outCollection);
                outCollection = mat2cell(outCollection, repelem(1, sOutCollection(1)), repelem(1, sOutCollection(2)));
                outCollection = [{'Mean';'Std'}, outCollection];
                obj.guihandles.bsResult = figure;
                t = uitable(obj.guihandles.bsResult, 'data', outCollection, 'columnName', [string('Type'), string('MeanLocs'), ticked]);
            end
            obj.setGuiParameters(pout);
            out.histogram=histogram;
            out.fit=pout;
        end
        
        function refit_callback(obj)
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        

    end
end




function pard=guidef(obj)
pard.text1.object=struct('String','parameters','Style','text');
pard.text1.position=[3,1];

pard.N0_fit.object=struct('String','N0','Style','radiobutton');
pard.N0_fit.position=[2,1];

pard.N0_v.object=struct('String','10','Style','edit');
pard.N0_v.position=[2,2];
pard.N0_v.isnumeric=1;


pard.pmature_fit.object=struct('String','p mature','Style','radiobutton');
pard.pmature_fit.position=[3,1];

pard.pmature_v.object=struct('String','.5','Style','edit');
pard.pmature_v.position=[3,2];
pard.pmature_v.isnumeric=1;


pard.pblink_fit.object=struct('String','p blink','Style','radiobutton');
pard.pblink_fit.position=[4,1];

pard.pblink_v.object=struct('String','.2','Style','edit');
pard.pblink_v.position=[4,2];
pard.pblink_v.isnumeric=1;


pard.monomer_fit.object=struct('String','monomer fraction','Style','radiobutton');
pard.monomer_fit.position=[5,1];

pard.monomer_v.object=struct('String','.2','Style','edit');
pard.monomer_v.position=[5,2];
pard.monomer_v.isnumeric=1;




pard.ncluster_fit.object=struct('String','n in cluster','Style','radiobutton');
pard.ncluster_fit.position=[6,1];

pard.ncluster_v.object=struct('String','0','Style','edit');
pard.ncluster_v.position=[6,2];
pard.ncluster_v.isnumeric=1;

pard.blinkmode.object=struct('Style','popupmenu','String','Poisson|Exponential|Noblink','Value',2);
pard.blinkmode.position=[7,1];
pard.blinkmode.Width=2;


pard.t1.object=struct('String','fitrange','Style','text');
pard.t1.position=[2,3];

pard.fitrange_min.object=struct('String','4','Style','edit');
pard.fitrange_min.position=[2,4];
pard.fitrange_min.Width=0.5;
pard.fitrange_max.object=struct('String','60','Style','edit');
pard.fitrange_max.position=[2,4.5];
pard.fitrange_max.Width=0.5;

pard.t2.object=struct('String','fit:','Style','text');
pard.t2.position=[4,3];
pard.t2.Width=0.25;
pard.fitselection.object=struct('String',{{'cumulative distribution','histogram','weighted histogram'}},'Style','popupmenu');
pard.fitselection.position=[4,3.25];
pard.fitselection.Width=1.75;

pard.bootstrap.object=struct('String','Bootstrap','Style','checkbox', 'Value', 0);
pard.bootstrap.position=[5,3];
pard.bootstrap.Width=1;

pard.bint.object=struct('String','binning','Style','text');
pard.bint.position=[6,3];
pard.bint.Width=1;
pard.bin.object=struct('String','1','Style','edit');
pard.bin.position=[6,4];
pard.bin.Width=1;

pard.plugininfo.name='fit brightness histogram';
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='fits histogram of localizaitons/cluster with an analytical model to determine number of proteins, maturation efficiency,  number of re-activations';
% pard.text2.object=struct('String','length scale nm','Style','text');
% pard.text2.position=[8,1];
% 
% pard.lengthscale.object=struct('String','15','Style','edit');
% pard.lengthscale.position=[8,2];
% pard.lengthscale.isnumeric=1;

% pard.results1='results 1';
% pard.results3='other results';


end