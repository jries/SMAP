classdef LocalizationData<interfaces.GuiParameterInterface
    % Stores all localization related information and provides access to
    % localization related functionality
    properties
        loc  %ungrouped localizations. locs.(field) is vector of lenth(number of localizations)
        grouploc %grouped localizations. Not calculated on the fly for performance
        layer %stores layer-specific information = filters.
        files %localization files and file meta data
        SE %siteexplorer object, linked here why?
        history={};
%         iscopy=true;
    end
    
    methods

        function obj=LocalizationData
            %constructor. calls .clear to initialize 
            obj.clear;
        end
%         function set.loc(obj,value)
%             obj.loc=value;
%         end
%         function loc=get.loc(obj)
%             loc=obj.loc;
%         end
        function clear(obj,part)
            %initializes properties to empty fields, sets files property to
            %no files
            if nargin<2||strcmpi(part,'all')
            obj.loc=[];
            obj.grouploc=[];
            obj.layer=[];
             obj.layer.filter=[];
             obj.layer.groupfilter=[];
            obj.files.filenumberEnd=0;
            obj.files.file=[];    
            obj.history={};
            if ~isempty(obj.SE)
                obj.SE.clear;
            else
                obj.SE=interfaces.SiteExplorer(obj);
                obj.SE.attachPar(obj.P);
                obj.SE.attachLocData(obj);
%                 obj.SE.locData=obj;
            end
            elseif strcmpi(part,'filter')
                obj.layer=[];
            end
            obj.inputParameters={'numberOfLayers','filtertable','layer','group_dx','group_dt'};
        end
        function addfile(obj,filename)
            %initializes obj.files.file(end+1) to default. 
             if nargin <2
                 filename='';
             end
            obj.files.filenumberEnd=obj.files.filenumberEnd+1;
            filenew=initfile(filename);
            filenew.number= obj.files.filenumberEnd;
            if obj.files.filenumberEnd==1
                obj.files.file=filenew;
            else
                obj.files.file(obj.files.filenumberEnd)=obj.files.file(1);
            obj.files.file(obj.files.filenumberEnd)=copyfields(obj.files.file(obj.files.filenumberEnd),filenew,fieldnames(obj.files.file(obj.files.filenumberEnd)));
            end
            if isempty(obj.loc)
                emp=zeros(0,1);
                locs=struct('frame',double(emp),'xnm',single(emp),'ynm',single(emp),'channel',single(emp),'locprecnm',single(emp),'filenumber',single(emp));
                obj.loc=locs;
                obj.grouploc=locs;
            end
            obj.SE.addFile(filenew.name,obj.files.filenumberEnd,obj.files.file(obj.files.filenumberEnd).info);
        end
        function setloc(obj,name, value,indused)
            %Add a field to all localizations. setloc(field, value) sets obj.loc.(name)=value
            if nargin<4
                obj.loc.(name)=real(value);
            else %only part of locs are calculated. Fill rest with zeros. Can be grouped or ungrouped.
                nu=length(obj.loc.frame);
%                 ng=length(obj.grouploc.frame);
                if length(indused)==nu %ungrouped data
                    v=zeros(nu,1,'single');
                    v(indused)=value;
                else
                    v=obj.grouped2ungrouped(indused,value);
                end
                obj.loc.(name)=real(v);
                obj.regroup;    
            end
            
            locfields=fieldnames(obj.loc);
            if ~isempty(obj.P)
            obj.setPar('locFields',locfields,'String');
            end
        end
        function addloc(obj,name,value)
            %addloc(field, value) Appends localizations to exisiting field, sets localizations
            %if not exisitng. 
            %FIX: can result in field vecotrs of different length. Pad with
            %zeros?
            if isempty(obj.loc)
                obj.loc.(name)=real(value);
                return
            end
            if nargin<3
                fn=fieldnames(obj.loc);
                
                value=zeros(size(obj.loc.(fn{1})),'single');
            end
            if isfield(obj.loc,name)
                l1=length(obj.loc.(name));
                l2=length(value);
%                 v1=obj.loc;
%                 value=[obj.loc.(name);real(value)];
                obj.loc.(name)(l1+1:l1+l2,1)=value;
            else
                fn=fieldnames(obj.loc);
                l1=length(obj.loc.(fn{1}));
                l2=length(value);
                if l1==l2
                    obj.loc.(name)=value;
                else
                
                    obj.loc.(name)(1:l1,1)=0;
                    obj.loc.(name)(l1+1:l1+l2,1)=value;
                end
%                 obj.setPar('locFields',fieldnames(obj.loc));
%                 obj.setloc(name,value)
            end
        end
        function [locsout,indout,hroi]=getloc(obj,fields,varargin)
        %returns subset of localizations. 
        % call from the interfaces.LocalizationData object that is attached to all
        % plugins as obj.locData
        % [locsout,indcombined,hroio]=obj.locData.getloc({filed1, field2,...},'Name1','parameter1','Name2',...)

        %output: 
        %locsout:   structure with fields corresponding ot the input field names
        %           containing all localizations that are requested by the name, parameter
        %           pairs.
        %indcombined    indices of all localizations that are returned
        %hroio      handle to image ROI


        %fields: any locData field, in addition: 
        % 'ingrouped','inungrouped': return a logical array that is true if the
        %           localization is displayed. For grouped and ungrouped localizations,
        %           respectively.
        % 'within': return true if inside the ROI.
        % inlayeru, inlayerg: return cell array with indices for each layer

        %Properties: Name, value pairs
        %'grouping': If to return grouped or ungrouped localizations
        %    'grouped': return grouped
        %    'ungrouped': return ungrouped 
        %    'layer'(default). For each layer determine the grouping from the
        %               setting in the layer
        %'channels': double vector of channels. Then localizations in the channels
        %               are returned
        %'filenumber': double vector of filenumbers
        %'position': If to return localizations in a defined ROI
        %       'all' (default): do not filter by position (xnm, ynm)
        %       'roi': Return localizations in the user-defined image ROI. If no
        %               ROI is defined, use the fov
        %       'fov': return localizations in the FoV (area currently displayed in
        %               the reconstructed image
        %        numerical vector: [centerx, centery, widhtx widthy] for rectangular ROI or 
        %                   [centerx, centery,radius] for circular ROI
        %        a 'interfaces.SEsites object: localizations in site.
        %'layer': localizations from which layer to use
        %           double number or vector: use these layers and filter the
        %           localizations according to the layer. if parameter 'layer' is
        %           not set, use all localizations without filtering
        %'removeFilter': cell array of filter names to remove 
        %'within': pass on indicies of localizations, those are considered in the
        %           output
        %'shiftxy', shifts x, y  (and z) as specified in a numerical vector with 2
        %           (3) entries
            if any(strcmp(varargin,'locselector'))
                [locsout,indout,hroi]=getlocs2ind(obj,fields,varargin{:});
            elseif any(strcmp(varargin,'ingrouped')) || any(strcmp(varargin,'inungrouped'))
                [locsout,indout,hroi]=getlocs2ind(obj,fields,varargin{:});
            else
                [locsout,indout,hroi]=getlocs(obj,fields,varargin{:});
%                 [locsout,indout,hroi]=getlocs2(obj,fields,varargin{:});

            end
        end
        
        function im=getimage(obj,layer)
            if nargin<2
                im=obj.getPar('sr_image');
            else
                if layer>length(obj.layer)
                    im=[];
                else
                    im=obj.layer(layer).images.srimage;
                end
            end
        end
            
        function mask=getRoiMask(obj)
            %returns binary imge of roi: mask=obj.getRoiMask
            hroi=obj.parameters.roihandle;
            mask=createMask(hroi);
        end
        
        function removelocs(obj,ind,removepart)
            %removes localizations with specified indices from all fields
            %obj.removelocs(indices,removepart)
            %optional removepart='loc' (default) or 'grouploc', determine which part to remove
            if nargin<3
                removepart='loc';
            end
            fn=fieldnames(obj.(removepart));
            for k=1:length(fn)
                obj.(removepart).(fn{k})(ind)=[];
            end
        end
        
        function filter=getFilter(obj,layer, grouping)
            %get filter structure from layer, 
            %optional: specify grouping=true or false, if not specified, take from layer
            %filter=obj.getFilter(layer, grouping)
            if nargin<3||isempty(grouping)
                grouping=obj.isgrouped(layer);
            end
            
            if layer<=length(obj.layer)
            if grouping
                filter=obj.layer(layer).groupfilter;
            else
                filter=obj.layer(layer).filter;
            end  
            else
                filter=[];
            end
        end
        
        
        function setFilter(obj,filter,layer,grouping) 
            %set filter structure for specific layer (from getFilter)
            %obj.setFilter(filter,layer,grouping)
            if nargin<4||isempty(grouping)
                grouping=obj.isgrouped(layer);
            end
            if grouping
                obj.layer(layer).groupfilter=filter;
            else
                obj.layer(layer).filter=filter;
            end  
        end
        
        function addhistory(obj,p)
            obj.history{end+1}=p;
        end
        
        function removeFilter(obj,field,layer)
            %remove specified filter from specified layer
            % removeFilter(field,layer)
            if layer<=length(obj.layer)
            obj.layer(layer).filter=myrmfield(obj.layer(layer).filter,field);
            obj.layer(layer).groupfilter=myrmfield(obj.layer(layer).groupfilter,field);
            end
        end
        
        function indin=inFilter(obj,layer,grouping)
            %returns indices of filtered localizations in specific layer.
            %Grouping (t/f) can be specified or is taken from layer
            %indices=inFilter(layer,grouping)
            if nargin<3||isempty(grouping)
                grouping=obj.isgrouped(layer);
            end
            
            if grouping && ~isempty(obj.grouploc)
                fn=fieldnames(obj.grouploc);
                indin=true(length(obj.grouploc.(fn{1})),1);
            else
                fn=fieldnames(obj.loc);
                indin=true(length(obj.loc.(fn{1})),1);
            end
            filters=obj.getFilter(layer,grouping); 
            if ~isempty(filters)
            fields=fieldnames(filters);
                if ~isempty(fields)
                    for k=1:length(fields)     
                        if size(indin)==size(filters.(fields{k}))
                            indin=indin&filters.(fields{k});
                        else
                            disp('size of filters does not match. Call interfaces.LocData.filter before. Tried to fix it')
%                             obj.filter(fields{k},layer)
%                             indin=inFilter(obj,layer,grouping);
%                             filters
%                             break       
                        end
                    end
                end
            end
        end
        
        function filter(obj,fields,layers,filtermode,minmax)
            %recaluclates all filters for specified layers.
            %filter(fields,layers,filtermode,filtervalues)
            %fields: cell array of fieldnames which to filter. If empty:
            %take all fields, also filter for channel
            %layers: all layers which to filter on. if empty: use all
            %layers
            %filtermode (optional): 'inlist' (filtervalues: list of integers to
            %keep),'minmax' (filtervalues: [minfilter maxfilter]). If not
            %specified: take from filter table
            
            if isempty(obj.loc)
                return
            end
            f=fieldnames(obj.loc);     
            if isempty(obj.loc.(f{1}))
                return
            end
            if nargin<3||isempty(layers) %filter all layers
                layers=1:obj.getPar('numberOfLayers');
            end
            
            if nargin<2||isempty(fields)%filter everything
                fields=fieldnames(obj.loc);
                 fields{end+1}='channel';
                 for k=1:length(layers)
                    obj.layer(layers(k)).groupfilter=[];
                    obj.layer(layers(k)).filter=[];
                 end
                 filterall=true;
            else
                if ~iscell(fields)
                fields={fields};
                end
                filterall=false;
            end

            for k=1:length(layers)  
                layerh=layers(k);
                s=obj.getPar('filtertable','layer',layerh);
                if ~isempty(s)
                fn=s(:,1);
                else
                    fn={};
                end
                for f=1:length(fields)
                    if ~isfield(obj.loc,fields{f})
                        disp(['cannot filter ' fields{f} '. Field in locData.loc is missing.'])
                        continue
                    end
                    invert=false;
                    ind=find(strcmp(fn,fields{f}));
                    if (~isempty(ind) && (s{ind,7})) ||(~filterall&&isfield(obj.loc,fields{f})&&nargin>3) %if field specified: filter for sure
                        if nargin<4
                            minmax=s(ind,[2,6]);
                            filtermode='minmax';
                        end
                        if size(s,2)>=8 && s{ind,8}
                            invert=true;
                        end
                    elseif strcmp(fields{f},'channel')
                        filtermode='inlist';
                        if nargin <5
                            minmax=obj.getPar('channels','layer',k);
                        end
                    else 
                        continue
                    end   
                    if invert
                        obj.layer(layerh).filter.(fields{f})=~anyfilter(obj.loc.(fields{f}),filtermode,minmax);
                    else
                        obj.layer(layerh).filter.(fields{f})=anyfilter(obj.loc.(fields{f}),filtermode,minmax);
                    end
                    if ~isempty(obj.grouploc)
                        if invert
                            obj.layer(layerh).groupfilter.(fields{f})=~anyfilter(obj.grouploc.(fields{f}),filtermode,minmax);
                        else
                            obj.layer(layerh).groupfilter.(fields{f})=anyfilter(obj.grouploc.(fields{f}),filtermode,minmax);
                        end
                        
                    end
                end
            end            
        end
            
        function g=isgrouped(obj,layer)
            %returns grouping state of layer. Should belong to GuiChannel,
            %but link here is often useful
            %g=isgrouped(layer)
            g=zeros(length(layer),1);
            
            for k=1:length(layer)
                if layer(k)<=length(obj.layer)
             p=obj.getPar(['layer' num2str(layer(k)) '_']);
             g(k)=logical(p.groupcheck);  
                end
            end
            
            
        end
        
        function saveloc=savelocs(obj,filename,goodind,additionalsave,grouping,excludesavefields,filenumber)
            %saves localization data to specific filename.  Returns saved structure:
            %saveloc=savelocs(filename,goodind)
            %goodind (optional): indices which to save (applies to all
            %fields). Default: save all localizations.
            %grouping: if true, saves grouped data.
            %excludesavefields: those fields not to save
            %filenumber: only save specific filenumber
            
            if nargin<5||isempty(grouping)
                grouping=false;
            end
            if grouping
                locext='grouploc';
            else
                locext='loc';
            end
            saveloc.loc=obj.(locext);
            saveloc.file=obj.files.file;
            saveloc.history=obj.history;
            
            fieldsremove={'original_channel','groupindex','numberInGroup','colorfield'};
            if nargin>5 && ~isempty(fieldsremove)
                fieldsremove=union(fieldsremove,excludesavefields);
            end
            saveloc.loc=myrmfield(saveloc.loc,fieldsremove);
            
            if nargin>6 && ~isempty(filenumber)%filenumber
                    saveloc.file=saveloc.file(filenumber);
                    saveloc.history=saveloc.history(filenumber);
                    saveloc.loc.filenumber=ones(size(obj.loc.filenumber)); % yu-le mod
                    if isempty(goodind)
                        goodind=obj.loc.filenumber==filenumber;
                    else
                        goodind=goodind&obj.loc.filenumber==filenumber;
                    end       
            end
            
            if nargin>2&&~isempty(goodind)
                fields=fieldnames(saveloc.loc);
                for k=1:length(fields)
                    saveloc.loc.(fields{k})=saveloc.loc.(fields{k})(goodind);
                end
            end
            if nargin>3&&~isempty(additionalsave)
                saveloc=copyfields(saveloc,additionalsave);
            end
            if ~isempty(obj.SE)
                if nargin>6 && ~isempty(filenumber) && ~isempty(obj.SE.sites)%filenumber
                    SEtemp = obj.SE.copy;
                    siteFilenumber = getFieldAsVector(SEtemp.sites,'info.filenumber');
                    lsiteFilenumber = siteFilenumber == filenumber;
                    SEtemp.sites = SEtemp.sites(lsiteFilenumber);
                    
                    cellFilenumber = getFieldAsVector(SEtemp.cells,'info.filenumber');
                    lcellFilenumber = cellFilenumber == filenumber;
                    SEtemp.cells = SEtemp.cells(lcellFilenumber);
                    
                    allFilenumber = getFieldAsVector(SEtemp.files,'ID');
                    loneFilenumber = allFilenumber == filenumber;
                    SEtemp.files = SEtemp.files(loneFilenumber);
                    saveloc.siteexplorer=SEtemp.save;
                else
                    saveloc.siteexplorer=obj.SE.save;
                end
            end
            saveloc=concentratefilelist(saveloc);
            try
                rg=obj.getPar('mainGui');
                parameters=rg.saveParameters;
            catch err
                err
                parameters=[];
            end
            fileformat.name='sml';
            out=struct('saveloc',saveloc,'fileformat',fileformat,'parameters',parameters);
            if isfield(obj.files.file(1),'transformation')
                for k=length(obj.files.file):-1:1
                    if ~isempty(obj.files.file(k).transformation)
                        out.transformation=obj.files.file(k).transformation;
                        break
                    end
                end
            end
            
            if nargin>1&&~isempty(filename)
                if ~contains(filename,'_sml.mat')
                    [path,file]=fileparts(filename);
                    filename=[path filesep file '_sml.mat'];
                end

                if 1 %obj.getGlobalSetting('saveas73') %now I use this to save with the old saver: large files are saved as v7.3, not as small parts.
                    v=saverightversion(filename,out,'-v7.3');
                    disp(['saved as version ' v])
                else
                    savematparts(filename,out,{'saveloc','loc'});
                end
%                 out=struct('saveloc',saveloc);
%                 saverightversion(filename,out);
%                 save(filename,'saveloc','-v7.3')
            end
        end
        
        function setLocData(obj,locData)
            %sets localizations to those from interfaces.LocalizationData object
            %setLocData(locData:interfaces.LocalizationData)
            obj.clear
            obj.loc=locData.loc;
            obj.files=locData.files;
            if myisfield(locData,'SE')
                obj.SE=locData.SE;
                obj.SE.attachLocData(obj);
                obj.SE.attachPar(obj.P);
            end
        end
        
        function addLocData(obj,locData)  
            %adds localization from interfaces.LocalizationData object
            %addLocData(locData:interfaces.LocalizationData). Pad fields
            %which are not shared by both objects with zeros
            if myisfield(locData,'loc')
                loch=locData.loc;
            else
                loch=locData;
            end
            if ~isfield(loch,'filenumber')||all(loch.filenumber<1)
                fn=fieldnames(loch);
                loch.filenumber=obj.files.filenumberEnd+0*loch.(fn{1});
            end
            
            fn2=fieldnames(loch);
            lennew=length(loch.(fn2{1}));
            if ~isempty(obj.loc)
                fn1=fieldnames(obj.loc);
                lenold=length(obj.loc.(fn1{1}));
                
            else
                fn1={};
                lenold=0;
            end           
            if lenold>0 %fill missing channels with 0
            fd=setdiff(fn2,fn1);
            for k=1:length(fd)
                if iscell(loch.(fd{k})) 
                    obj.loc.(fd{k}){lenold}=[];
                else
                    obj.loc.(fd{k})=zeros(lenold,1,'like',loch.(fd{k}));
                end
            end  
            end
            %add new data
            for k=1:length(fn2)
                obj.addloc(fn2{k},loch.(fn2{k}))
            end     
            %add zeros for the missing data
            fd=setdiff(fn1,fn2);
            for k=1:length(fd)
                if iscell(obj.loc.(fd{k})) 
                    obj.loc.(fd{k}){lennew}=[];
                else
                obj.addloc(fd{k},zeros(lennew,1,'like',obj.loc.(fd{k})));
                end
            end                     
        end
        
        function sort(obj,varargin)
            %sorts localizations according to field1, field2,...
            %sort(field1,field2,...)
            
            %ungrouped
            sortm=zeros(length(obj.loc.xnm),length(varargin));
            for k=1:length(varargin)
                if isfield(obj.loc,varargin{k})
                    sortm(:,k)=(obj.loc.(varargin{k}));
                end
            end
            [~,sortind]=sortrows(sortm);
            fn=fieldnames(obj.loc);
            for k=1:length(fn)
                obj.loc.(fn{k})=obj.loc.(fn{k})(sortind);
            end
            
            %grouped
            sortm=zeros(length(obj.grouploc.xnm),length(varargin));
            for k=1:length(varargin)
                if isfield(obj.grouploc,varargin{k})
                    sortm(:,k)=(obj.grouploc.(varargin{k}));
                end
            end
            [~,sortind]=sortrows(sortm);
            fn=fieldnames(obj.grouploc);
            for k=1:length(fn)
                obj.grouploc.(fn{k})=obj.grouploc.(fn{k})(sortind);
            end
        end

        function locout= copy(obj,fields,indin,clearfilter)
            %creastes a deep copy of localization data.
            
            locout=interfaces.LocalizationData;
           locout.P=obj.P;
            locout.files=obj.files;
            locout.SE=obj.SE.copy;
             locout.SE.locData=locout;
             locout.history=obj.history;
             
            if nargin<2
                locout.loc=obj.loc;
                locout.grouploc=obj.grouploc;

    %             locout.iscopy=true;
                for k=1:length(obj.layer)
                    locout.layer(k).filter=obj.layer(k).filter;
                    locout.layer(k).groupfilter=obj.layer(k).groupfilter;
                end
            else %specify fields and inices
                if nargin<4
                    clearfilter=false;
                end
 
                fn=fieldnames(obj.loc);
                if isempty(fn)
                    return;
                end
                
                numlocs=length(obj.loc.(fn{1}));
                numglocs=length(obj.grouploc.(fn{1}));
                
                if nargin<3 || isempty(indin)
                    indu=true(numlocs,1);
                    indg=true(numglocs,1);
                else
                    if length(indin)==numlocs
                        indu=indin;
                        indg=indungrouped2grouped(indin,obj.loc.groupindex);
                    else
                        indg=indin;
                        indu=indgrouped2ungrouped(indin,obj.loc.groupindex);
                    end
                end
                if isempty(fields)
                    fields=fn;
                else
                    fields=intersect(fields,fn);
                end
                for k=1:length(fields)
                    locout.loc.(fields{k})=obj.loc.(fields{k})(indu);
                end
                for k=1:length(fields)
                    locout.grouploc.(fields{k})=obj.grouploc.(fields{k})(indg);
                end
%                 fields=fieldnames(obj.grouploc);
                if clearfilter
                    for l=1:length(obj.layer)
                    locout.layer(l).filter=[];
                    locout.layer(l).groupfilter=[];
                    end
                else
                    

    %                 l=length(obj.layer);
                    for l=1:length(obj.layer)
                        fields=fieldnames(obj.layer(l).filter);
                        for k=1:length(fields)
                            locout.layer(l).filter.(fields{k})=obj.layer(l).filter.(fields{k})(indu);
                        end        
                        fields=fieldnames(obj.layer(l).groupfilter);
                        for k=1:length(fields)
                            locout.layer(l).groupfilter.(fields{k})=obj.layer(l).groupfilter.(fields{k})(indg);
                        end                
                    end
                end
                
            end
            
        end
        
        function regroup(obj,dx,dt,mode)  
            %groups localization data.  
            %regroup(dx,dt): dx: maximum distance to be grouped. dt:
            %maximum number of frames localzaiton can be dark. if not
            %specified: taken from global parameters 'group_dx' and
            %'group_dt'
            %mode: 1. fix, 2. use locprec. Then dx is minmax vector
            if nargin<2
                dx=obj.getPar('group_dx');
            end
            if nargin<3
                dt=obj.getPar('group_dt');
            end
            if nargin<4
                mode=obj.getPar('group_mode');
                if mode==2 %locprec
                    dx=obj.getPar('group_dxminmax');
                end
            end
            
            if mode==3 %no grouping
                obj.grouploc=obj.loc;
                obj.filter
                return
            end
%             obj.setPar('status','grouping')
            grouper=Grouper;
            grouper.attachLocData(obj);
            groupfields={'frame','xnm','ynm','locprecnm','filenumber','channel'};
            if isfield(obj.loc,'class')
                groupfields{end+1}='class';
            end
            grouper.mode=mode;
            grouper.connect(dx,dt,groupfields{:})
            grouper.combine;  
            obj.filter;
%             locfields=fieldnames(obj.loc);
%             if ~obj.iscopy
%                 obj.setPar('locFields',locfields);
%             end
%             obj.setPar('status','grouping done')
        end  
        function [out,ind]=grouped2ungrouped(obj,indgrouped,values)
            % grouped2ungrouped writes values corresponding to grouped localizaitons in
            % indgrouped to corresponding values in ungrouped data.
            [indexgsort, indg]=sort(obj.grouploc.groupindex(indgrouped));
            [indexusort, indu]=sort(obj.loc.groupindex);
            vu=zeros(length(indexusort),1,'like',values);
            ku=1;
            for k=1:length(indexgsort)
                while ku<=length(indexusort)&&(indexusort(ku)<indexgsort(k))
                    ku=ku+1;
                end
%                 if (indexusort(ku)>indexgsort(k))
%                     continue;
%                 end
                ku2=ku;
                while ku2<=length(indexusort)&&(indexusort(ku2)==indexgsort(k))
                    vu(indu(ku2))=values(indg(k));
                    ku2=ku2+1;
                end
                ku=ku2; 
            end
            
            out=vu;
        end
    end    
end