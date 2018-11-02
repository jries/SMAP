classdef density_calculator<interfaces.DialogProcessor
    %  Copyright (c)2017 Ries Lab, European Molecular Biology Laboratory,
    %  Heidelberg. This file is part of Single Molecule Analysis Platform (SMAP).
    
    % density_calculator looks at the neighborhood and counts number of
    % neighbours. locData.clusterdensity=neighbours
    properties
        useind
    end
    methods
        function obj=density_calculator(varargin)    
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters={'sr_pixrec','numberOfLayers','sr_pos','sr_size','layers','sr_layerson'};
            obj.history=true;    
            obj.showresults=false;
        end
        
        function out=run(obj,p)
            out=[];
            if p.allfiles
                fieldremref={'filenumber'};
                fieldremall={'channel'};
                indl=1;
            else 
                fieldremref={};
                fieldremall={'filenumber','channel'};
            end
            
            
            %obj.locData.sort('filenumber','channel','xnm');
            %problem on filtered data, somehow it gets confused.
            %[locs,indin]=obj.locData.getloc({'xnm','ynm','channel','frame','znm','ingrouped','inungrouped'},'position','all');
            if contains(p.countwhat.selection,'layer')
                activelayers=find(p.sr_layerson);
                for k=1:length(activelayers)
                    fieldrem=setdiff(fieldnames(obj.locData.layer(k).filter),fieldremall)';
%                     fieldrem=setdiff(fieldrem,fieldremref);
                    locallh=obj.locData.getloc({'xnm','ynm','filenumber','channel','frame','znm','ingrouped','inungrouped'},'layer',activelayers(k),'position','all','removeFilter',fieldrem,'grouping','ungrouped');
                    locrefh=obj.locData.getloc({'xnm','ynm','filenumber','channel','frame','znm','ingrouped','inungrouped'},'layer',activelayers(k),'position','all','removeFilter',fieldremref);
                    fiu=find(locallh.inungrouped);
                    if p.allfiles
                        for f=1:max(locallh.filenumber)
                            infiles=locallh.filenumber==f;
                            locall(indl)=copystructReduce(locallh,infiles);
                            
                            
                            fiu2=fiu(infiles);
                            inung=false(length(locallh.inungrouped),1);
                            inung(fiu2)=true;
                            locall(indl).inungrouped=inung;
                            
                            infilesr=locrefh.filenumber==f;
                            locref(indl)=copystructReduce(locrefh,infilesr);
%                             locref(indl).inungrouped=locallh.inungrouped&infiles;                           
                            indl=indl+1;
                        end
                    else   
                        locall(k)=locallh;
                        locref(k)=locrefh;
                    end
                end
            else
                locall=obj.locData.getloc({'xnm','ynm','filenumber','channel','frame','znm','ingrouped','inungrouped'},'grouping','ungrouped');
                if ~isempty(obj.useind)
                    locref=copystructReduce(locall,obj.useind);
                else
                locref=locall;
                end
%                 locall=obj.locData.loc;
%                 locref=obj.locData.loc;
            end
            
           
                
            neighbourstotall=zeros(length(obj.locData.loc.xnm),1);
%            [locs,indin]=obj.locData.getloc({'xnm','ynm','filenumber','channel','frame','znm','ingrouped','inungrouped'},'layer',find(p.sr_layerson),'position','all');
            for lay=1:length(locall)
                neighbourstot=zeros(length(obj.locData.loc.xnm),1);
                sortmall=horzcat(locall(lay).xnm,(1:length(locall(lay).frame))');
                [sortedm,sortind]=sortrows(sortmall);
                xa=locall(lay).xnm(sortind);
                ya=locall(lay).ynm(sortind);
                
                sortmallr=horzcat(locref(lay).xnm,(1:length(locref(lay).frame))');
                [sortedmr,sortindr]=sortrows(sortmallr);
                xr=locref(lay).xnm(sortindr);
                yr=locref(lay).ynm(sortindr);

                dx=p.countingsize_xy;
                dz=p.countingsize_z;
                if isfield(locall(lay),'znm') && ~isempty(locall(lay).znm)&&~isempty(locref(lay).znm)
                    za=locall(lay).znm(sortind);
                    zr=locref(lay).znm(sortindr);
                    if p.countingregion.Value==1 %Gauss
                        neighbours=countneighbours3DGauss(double(xa),double(ya),double(za),double(dx),double(dz));
                        disp('not impolemented for Gauss')
                    else
                        neighbours=countneighbours3Dcirc2(double(xa),double(ya),double(za),double(xr),double(yr),double(zr),double(dx),double(dz));
                    end
                else
                    if p.countingregion.Value==1 %Gauss
                        neighbours=countneighbours2DGauss(double(xa),double(ya),double(dx));
                        disp('not impolemented for Gauss')
                    else
                        neighbours=countneighbours2Dcirc2(double(xa),double(ya),double(xr),double(yr),double(dx));
                    end
                end
                [~,sortbackind]=sort(sortedm(:,2));
                neighboursback=neighbours(sortbackind);
%                 if sum(locall(lay).-inungrouped)==length(locall(lay).xnm) %ungrouped data
                     neighbourstot(locall(lay).inungrouped)=neighboursback;
%                 else
%                     nbbackug=obj.locData.grouped2ungrouped(locall(lay).ingrouped,neighboursback);
%                     neighbourstot=nbbackug;
%                 end
                neighbourstotall=max(neighbourstotall,neighbourstot);
            end
            

            obj.locData.setloc('clusterdensity',single(neighbourstotall));
%             obj.locData.sort('filenumber','channel','frame');
            obj.locData.regroup;
            obj.setPar('locFields',fieldnames(obj.locData.loc));  
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        

    end
end




function pard=guidef
pard.countingregion.object=struct('String','Gauss wighted counting |circle counting','Style','popupmenu','Value',2);
pard.countingregion.object.TooltipString=sprintf('count Gauss-weighted locs (more accurate) or locs in circle/cylinder (faster)');
pard.countingregion.position=[1,1];
pard.countingregion.Width=2;

pard.countwhat.object=struct('String','all|layers','Style','popupmenu','Value',2);
pard.countwhat.object.TooltipString=sprintf('count all localizations or per layer (use visible ones as reference)');
pard.countwhat.position=[2,1];
pard.countwhat.Width=1;

pard.allfiles.object=struct('String','on all files','Style','checkbox','Value',0);
pard.allfiles.object.TooltipString=sprintf('Apply on all files together');
pard.allfiles.position=[2,2];
pard.allfiles.Width=1;

pard.texta.object=struct('String','size in x,y (nm)','Style','text');
pard.texta.position=[3,1];
pard.countingsize_xy.object=struct('String','12','Style','edit');
pard.countingsize_xy.position=[3,2];
pard.countingsize_xy.isnumeric=1;
pard.countingsize_xy.object.TooltipString=sprintf('radius of circle or sigma of gauss in lateral direction');


pard.text1.object=struct('String','size in z (nm)','Style','text');
pard.text1.position=[4,1];
pard.countingsize_z.object=struct('String','24','Style','edit');
pard.countingsize_z.position=[4,2];
pard.countingsize_z.isnumeric=1;
pard.countingsize_z.object.TooltipString=sprintf('size of cylinder or sigma of gauss in z direction');
pard.plugininfo.name='density calculator (number of neighbours)';
pard.plugininfo.description= 'density_calculator looks at the neighborhood and counts number of neighbours. If grouped or ungrouped data is used depends on setting in layers.';
pard.plugininfo.type='ProcessorPlugin';

end