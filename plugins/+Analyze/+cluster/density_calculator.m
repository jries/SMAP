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
            obj.inputParameters={'sr_pixrec','numberOfLayers','sr_pos','sr_size','sr_layerson'};
            obj.history=true;    
            obj.showresults=false;
        end
        
        function out=run(obj,p)
            out=[];
            
            % get locs to be counted
            activelayers=find(p.sr_layerson);
            if p.excludelayer
                targetlayers=setdiff(activelayers,p.excludelayernumber);
            else
                targetlayers=activelayers;
            end
            switch p.fileselector.Value
                case 1 %current
                    fieldremt={};
                case 2 %all files
                    fieldremt={'filenumber'};
                case 3 %individually
                    lf=obj.locData.getloc('filenumber');
                    allf=1:max(lf.filenumber);
                    ph=p;
                    ph.fileselector.Value=1;
                    for f=allf
                        ph.layer1_.ch_filelist.Value=f;
                        obj.locData.filter('filenumber',[],'inlist',f)
                        obj.run(ph); %recursive 
                    end
                    return
            end   
            targetlocs=obj.locData.getloc({'xnm','ynm','filenumber','channel','frame','znm','ingrouped','inungrouped'},'layer',targetlayers,'position','roi','removeFilter',fieldremt);
            % get locs for which to count: indices and position
            switch p.refselection.selection
                case 'filtered'
                    if p.filtermore
                        reflayers=intersect(activelayers,p.layern);
                    else
                        reflayers=activelayers;
                    end
                    [reflocs,indref]=obj.locData.getloc({'xnm','ynm','filenumber','channel','frame','znm','ingrouped','inungrouped',p.newfield},'layer',reflayers,'position','roi','removeFilter',fieldremt,'grouping','ungrouped');
                    indgood=true(size(reflocs.xnm));
                case 'all in ROI'
                    [reflocs,indref]=obj.locData.getloc({'xnm','ynm','filenumber','channel','frame','znm','ingrouped','inungrouped',p.newfield},'position','roi','grouping','ungrouped');
                    currentfilenumber=p.layer1_.ch_filelist.Value;
                    indgood=reflocs.filenumber==currentfilenumber;
                    %file number needs to be applied manually
            end
            if p.filtermore
                inchannel=any(reflocs.channel==p.channels,2);
                indgood=indgood&inchannel;
            end
            
            if ~isempty(reflocs.znm)&&~isempty(targetlocs.znm)  %3D
                refcoord=[reflocs.xnm(indgood),reflocs.ynm(indgood),reflocs.znm(indgood)];
                targetcoord=[targetlocs.xnm,targetlocs.ynm,targetlocs.znm];
                dxyz=[p.countingsize_xy, p.countingsize_z];
            else
                refcoord=[reflocs.xnm(indgood) reflocs.ynm(indgood)];
                targetcoord=[targetlocs.xnm,targetlocs.ynm];
                dxyz=[p.countingsize_xy];
            end
            
%             figure(88)
%             plot(refcoord(:,1),refcoord(:,2),'.',targetcoord(:,1),targetcoord(:,2),'.')
%             drawnow
            
            nn=countneighbours(refcoord,targetcoord,dxyz,p.countingregion.Value);
            findref=find(indref);
            if isempty(reflocs.(p.newfield)) %does not exist
                alln=zeros(size(indref),'single')-1;
                alln(findref(indgood))=nn;
                obj.locData.loc.(p.newfield)=alln;
            else
                obj.locData.loc.(p.newfield)(findref(indgood))=nn;
            end
            
            obj.locData.regroup;
            
%             return
%             if p.allfiles
%                 fieldremref={'filenumber'};
%                 fieldremall={'channel'};
%                 indl=1;
%             else 
%                 fieldremref={};
%                 fieldremall={'filenumber','channel'};
%             end
%             
%             
%             %obj.locData.sort('filenumber','channel','xnm');
%             %problem on filtered data, somehow it gets confused.
%             %[locs,indin]=obj.locData.getloc({'xnm','ynm','channel','frame','znm','ingrouped','inungrouped'},'position','all');
%             if contains(p.countwhat.selection,'layer')
%                 
%                 for k=1:length(activelayers)
%                     fieldrem=setdiff(fieldnames(obj.locData.layer(k).filter),fieldremall)';
% %                     fieldrem=setdiff(fieldrem,fieldremref);
%                     locallh=obj.locData.getloc({'xnm','ynm','filenumber','channel','frame','znm','ingrouped','inungrouped'},'layer',activelayers(k),'position','all','removeFilter',fieldrem,'grouping','ungrouped');
%                     locrefh=obj.locData.getloc({'xnm','ynm','filenumber','channel','frame','znm','ingrouped','inungrouped'},'layer',activelayers(k),'position','all','removeFilter',fieldremref);
%                     fiu=find(locallh.inungrouped);
%                     if p.allfiles
%                         for f=1:max(locallh.filenumber)
%                             infiles=locallh.filenumber==f;
%                             locall(indl)=copystructReduce(locallh,infiles);
%                             
%                             
%                             fiu2=fiu(infiles);
%                             inung=false(length(locallh.inungrouped),1);
%                             inung(fiu2)=true;
%                             locall(indl).inungrouped=inung;
%                             
%                             infilesr=locrefh.filenumber==f;
%                             locref(indl)=copystructReduce(locrefh,infilesr);
% %                             locref(indl).inungrouped=locallh.inungrouped&infiles;                           
%                             indl=indl+1;
%                         end
%                     else   
%                         locall(k)=locallh;
%                         locref(k)=locrefh;
%                     end
%                 end
%             else
%                 locall=obj.locData.getloc({'xnm','ynm','filenumber','channel','frame','znm','ingrouped','inungrouped'},'grouping','ungrouped');
%                 if ~isempty(obj.useind)
%                     locref=copystructReduce(locall,obj.useind);
%                 else
%                 locref=locall;
%                 end
% %                 locall=obj.locData.loc;
% %                 locref=obj.locData.loc;
%             end
%             
%            
%                 
%             neighbourstotall=zeros(length(obj.locData.loc.xnm),1);
% %            [locs,indin]=obj.locData.getloc({'xnm','ynm','filenumber','channel','frame','znm','ingrouped','inungrouped'},'layer',find(p.sr_layerson),'position','all');
%             for lay=1:length(locall)
%                 neighbourstot=zeros(length(obj.locData.loc.xnm),1);
%                 sortmall=horzcat(locall(lay).xnm,(1:length(locall(lay).frame))');
%                 [sortedm,sortind]=sortrows(sortmall);
%                 xa=locall(lay).xnm(sortind);
%                 ya=locall(lay).ynm(sortind);
%                 
%                 sortmallr=horzcat(locref(lay).xnm,(1:length(locref(lay).frame))');
%                 [sortedmr,sortindr]=sortrows(sortmallr);
%                 xr=locref(lay).xnm(sortindr);
%                 yr=locref(lay).ynm(sortindr);
% 
%                 dx=p.countingsize_xy;
%                 dz=p.countingsize_z;
%                 if isfield(locall(lay),'znm') && ~isempty(locall(lay).znm)&&~isempty(locref(lay).znm)
%                     za=locall(lay).znm(sortind);
%                     zr=locref(lay).znm(sortindr);
%                     if p.countingregion.Value==1 %Gauss
%                         countf=@countneighbours3DGauss2;
%                     else
%                         countf=@countneighbours3Dcirc2;
%                     end
%                     neighbours=countf(double(xa),double(ya),double(za),double(xr),double(yr),double(zr),double(dx),double(dz));
% 
%                 else
%                     if p.countingregion.Value==1 %Gauss
%                         countf=@countneighbours2DGauss2;
%                     else
%                         countf=@countneighbours2Dcirc2;
%                     end
%                     neighbours=countf(double(xa),double(ya),double(xr),double(yr),double(dx));
%                 end
%                 [~,sortbackind]=sort(sortedm(:,2));
%                 neighboursback=neighbours(sortbackind);
% %                 if sum(locall(lay).-inungrouped)==length(locall(lay).xnm) %ungrouped data
%                      neighbourstot(locall(lay).inungrouped)=neighboursback;
% %                 else
% %                     nbbackug=obj.locData.grouped2ungrouped(locall(lay).ingrouped,neighboursback);
% %                     neighbourstot=nbbackug;
% %                 end
%                 neighbourstotall=max(neighbourstotall,neighbourstot);
%             end
%             
% 
%             obj.locData.setloc('clusterdensity',single(neighbourstotall));
% %             obj.locData.sort('filenumber','channel','frame');
%             obj.locData.regroup;
%             obj.setPar('locFields',fieldnames(obj.locData.loc));  
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        

    end
end


function neighboursback=countneighbours(ref,target,dxyz,region)
sortmallr=horzcat(ref(:,1),(1:size(ref,1))');
[sortedmr,sortindr]=sortrows(sortmallr);
refsorted=ref(sortindr,:);

sortmallt=horzcat(target(:,1),(1:size(target,1))');
[sortedmt,sortindt]=sortrows(sortmallt);
targetsorted=target(sortindt,:);

if size(ref,2)==3 %3D
    if region==1 %Gauss
        countf=@countneighbours3DGauss2;
    else
        countf=@countneighbours3Dcirc2;
    end
    neighbours=countf(double(refsorted(:,1)),double(refsorted(:,2)),double(refsorted(:,3)),...
            double(targetsorted(:,1)),double(targetsorted(:,2)),double(targetsorted(:,3)),double(dxyz(1)),double(dxyz(end)));   
else
    if region==1 %Gauss
        countf=@countneighbours2DGauss2;
    else
        countf=@countneighbours2Dcirc2;
    end
    neighbours=countf(double(refsorted(:,1)),double(refsorted(:,2)),...
            double(targetsorted(:,1)),double(targetsorted(:,2)),double(dxyz(1)));
end
[~,sortbackind]=sort(sortedmr(:,2));
neighboursback=neighbours(sortbackind);
% xa:reference
                
%                 sortmallr=horzcat(locref(lay).xnm,(1:length(locref(lay).frame))');
%                 [sortedmr,sortindr]=sortrows(sortmallr);
%                 xr=locref(lay).xnm(sortindr);
%                 yr=locref(lay).ynm(sortindr);
% 
%                 dx=p.countingsize_xy;
%                 dz=p.countingsize_z;
%                 if isfield(locall(lay),'znm') && ~isempty(locall(lay).znm)&&~isempty(locref(lay).znm)
%                     za=locall(lay).znm(sortind);
%                     zr=locref(lay).znm(sortindr);
%                     if p.countingregion.Value==1 %Gauss
%                         neighbours=countneighbours3DGauss(double(xa),double(ya),double(za),double(dx),double(dz));
%                         disp('not impolemented for Gauss')
%                     else
%                         neighbours=countneighbours3Dcirc2(double(xa),double(ya),double(za),double(xr),double(yr),double(zr),double(dx),double(dz));
%                     end
%                 else
%                     if p.countingregion.Value==1 %Gauss
%                         neighbours=countneighbours2DGauss(double(xa),double(ya),double(dx));
%                         disp('not impolemented for Gauss')
%                     else
%                         neighbours=countneighbours2Dcirc2(double(xa),double(ya),double(xr),double(yr),double(dx));
%                     end
%                 end
%                 [~,sortbackind]=sort(sortedm(:,2));
%                 neighboursback=neighbours(sortbackind);
% %                 if sum(locall(lay).-inungrouped)==length(locall(lay).xnm) %ungrouped data
%                      neighbourstot(locall(lay).inungrouped)=neighboursback;
% %                 else
% %                     nbbackug=obj.locData.grouped2ungrouped(locall(lay).ingrouped,neighboursback);
% %                     neighbourstot=nbbackug;
% %                 end
%                 neighbourstotall=max(neighbourstotall,neighbourstot);
end

function pard=guidef(obj)
pard.countingregion.object=struct('String','Gauss wighted counting |circle counting','Style','popupmenu','Value',2);
pard.countingregion.object.TooltipString=sprintf('count Gauss-weighted locs (more accurate) or locs in circle/cylinder (faster)');
pard.countingregion.position=[1,1];
pard.countingregion.Width=1.5;

pard.texta.object=struct('String','size in x,y (nm)','Style','text');
pard.texta.position=[2,1];
pard.countingsize_xy.object=struct('String','12','Style','edit');
pard.countingsize_xy.position=[2,2];
pard.countingsize_xy.isnumeric=1;
pard.countingsize_xy.object.TooltipString=sprintf('radius of circle or sigma of gauss in lateral direction');
pard.countingsize_xy.Width=0.5;

pard.text1.object=struct('String','size in z (nm)','Style','text');
pard.text1.position=[3,1];
pard.countingsize_z.object=struct('String','24','Style','edit');
pard.countingsize_z.position=[3,2];
pard.countingsize_z.isnumeric=1;
pard.countingsize_z.object.TooltipString=sprintf('size of cylinder or sigma of gauss in z direction');
pard.countingsize_z.Width=0.5;

pard.newfieldt.object=struct('String','Field name','Style','text');
pard.newfieldt.position=[4,1];
pard.newfield.object=struct('String','clusterdensity','Style','edit');
pard.newfield.position=[4,2];


pard.targett.object=struct('String','Only filtered locs in ROI are counted','Style','text');
pard.targett.position=[1,3];
pard.targett.Width=2;

p(1).value=0; p(1).on={}; p(1).off={'excludelayernumber'};
p(2).value=1; p(2).on={'excludelayernumber'}; p(2).off={};
pard.excludelayer.object=struct('String','Exclude layer','Style','checkbox','Value',0,'Callback',{{@obj.switchvisible,p}});
pard.excludelayer.position=[2,3.5];
pard.excludelayernumber.object=struct('String','1','Style','edit','Visible','off');
pard.excludelayernumber.position=[2,4.5];
pard.excludelayernumber.Width=0.5;

pard.reft.object=struct('String','Count neighbours for these localizations:','Style','text');
pard.reft.position=[3,3];
pard.reft.Width=2;

pard.refselection.object=struct('String',{{'filtered','all in ROI'}},'Style','popupmenu');
pard.refselection.position=[4,3.5];
pard.refselection.Width=1;

p(1).value=0; p(1).on={}; p(1).off={'channelst','channels','layerst','layern'};
p(2).value=1; p(2).on={'channelst','channels','layerst','layern'}; p(2).off={};
pard.filtermore.object=struct('String','filter in addition (keep only):','Style','checkbox','Callback',{{@obj.switchvisible,p}});
pard.filtermore.position=[5,3.5];
pard.filtermore.Width=2;

pard.channelst.object=struct('String','Channels','Style','text','Visible','off');
pard.channelst.position=[6,3.5];
pard.channels.object=struct('String','0','Style','edit','Visible','off');
pard.channels.position=[6,4.5];
pard.channels.Width=0.5;

pard.layerst.object=struct('String','Layers (filtered only)','Style','text','Visible','off');
pard.layerst.position=[7,3.5];
pard.layern.object=struct('Style','edit','String','1','Visible','off');
pard.layern.position=[7,4.5];
pard.layern.Width=0.5;

pard.fileselector.object=struct('String',{{'current filtered','all files together','all files separately'}},'Style','popupmenu');
pard.fileselector.position=[6,1];
pard.fileselector.Width=1.5;

pard.plugininfo.name='density calculator (number of neighbours)';
pard.plugininfo.description= 'density_calculator looks at the neighborhood of each localizations and counts number of neighbours in a defined region If grouped or ungrouped data is used depends on setting in layers.';
pard.plugininfo.type='ProcessorPlugin';

end