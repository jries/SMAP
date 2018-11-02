classdef Locstatistics<interfaces.DialogProcessor
    % Locstatistics calculates all kind of statistics for localization
    % data.
    methods
        function obj=Locstatistics(varargin)           
            obj@interfaces.DialogProcessor(varargin{:}) 
            obj.inputParameters={'sr_layerson','layernames'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)
        fields={'filenumber','frame','phot','locprecnm','znm','PSFxnm','locprecznm','numberInGroup','bg'};
        if p.useroi
            position='roi';
        else
            position='all';
        end

        layers=find(p.sr_layerson);
            if p.filter
                for m=length(layers):-1:1
                    locs{m}=obj.locData.getloc(fields,'layer',layers(m),'position',position);
                    modetxt{m}=['layer' num2str(layers(m))];
                end
            else
                locs{2}=obj.locData.getloc(fields,'position',position,'grouping','grouped');
                locs{1}=obj.locData.getloc(fields,'position',position,'grouping','ungrouped');
                modetxt{2}='grouped';
                modetxt{1}='ungrouped';
            end
            
            out=make_statistics2(locs,p,true);
            tcl=sprintf(' \t Nloc  \t muphot \t locprecmax  \t locprecmedian  \t locprecrise  \t mulifetime  \t bgmean ');
            if isfield(out,'PSFxnm')
                tcl=[tcl sprintf('\t PSFxnm')];
            end
            if isfield(out,'locprecznm')
                tcl=[tcl sprintf('\t locprecznm max') ];
            end
            tcl=[tcl 13];
            for k=1:length(out.photons.Nloc)
                th=sprintf([ modetxt{k} '\t' num2str(out.photons.Nloc(k)) '\t' num2str(out.photons.mu(k)) '\t'  num2str(out.locprec.max(k)) '\t'...
                    num2str(out.locprec.median(k)) '\t' num2str(out.locprec.rising(k)) '\t' num2str(out.lifetime.mu(k)) '\t'...
                    num2str(out.background.mean(k))]);
                if isfield(out,'PSFxnm')
                    th=[th 9 num2str(out.PSFxnm.max(k))];
                end
                if isfield(out,'locprecznm')
                    th=[th 9 num2str(out.locprecznm.max(k))];
                end
                tcl=[tcl th 13];
            end
            out.clipboard=tcl;

        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef
pard.useroi.object=struct('String','use Roi','Style','checkbox','Value',1);
pard.useroi.position=[1,1];

pard.filter.object=struct('String','use layers/filters','Style','checkbox','Value',1);
pard.filter.position=[1,2];

pard.overview.object=struct('String','plot overview','Style','checkbox','Value',0);
pard.overview.position=[1,3];

pard.tphot.object=struct('String','photon range:','Style','text');
pard.tphot.position=[3,1];
pard.tphot.Width=1.5;

pard.photrange.object=struct('String','800 10000','Style','edit');
pard.photrange.position=[3,2.5];

pard.tlt.object=struct('String','lifetime range (frames):','Style','text');
pard.tlt.position=[4,1];
pard.tlt.Width=1.5;

pard.lifetimerange.object=struct('String','1 30','Style','edit');
pard.lifetimerange.position=[4,2.5];

pard.plugininfo.name='Statistics';
pard.plugininfo.description=sprintf(['Locstatistics calculates all kind of statistics for localization data.\n'...
    'photons: N: number of localizations. <P>: mean.  mu: decay constant of exponential fit. \n'...
    'locprec: max, median and position of rising edge. \n'...
    'lifetime: how many frames does a fluorophore live (from grouping). mu: from exponential fit.\n'...
    'background: mean\n'...
    'either znm or PSFxnm.\n'...
    'locprecznm \n'...
    'frames: number of lcoalizations vs. frame. To see when localizations drop off.']);
pard.plugininfo.type='ProcessorPlugin';
end