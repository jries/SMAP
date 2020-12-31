classdef Grouper< interfaces.LocDataInterface
%     Combines (merges, groups) localizations which persit over several
%     frames into a single localization. Accessible in the File GUI.
    properties
        combinemodes
        indsortlist
        dllistsort
        mode=1;
    end
    methods
        function obj=Grouper(varargin)            
            if nargin>0
                obj.attachLocData(varargin{1});
            end
%             if nargin>1
%                 obj.attachPar(varargin{2});
%             end
%             obj.inputParameters={'group_dx','group_dt'};
          %define: meanx, meany, meanz to go with xerr, yerr, zerr when
          %weighting, same for locprecx, locprecy,locprecz.
          % or not: as averaged locs will have same x,y,z only relative
          % locprecs count, then it should not matter much which we take or
          % if we just take 1/sqrt(N). For performance: keep only one
          % weighting. E.g. locprecnm or N
            obj.combinemodes.xnm='mean';
            obj.combinemodes.ynm='mean';
            obj.combinemodes.znm='mean';
            obj.combinemodes.xpix='mean';
            obj.combinemodes.ypix='mean';
            obj.combinemodes.zpix='mean';
            obj.combinemodes.znm_SALM='mean';
            obj.combinemodes.znm_asSALM='mean';
            obj.combinemodes.znm_a='mean';
            obj.combinemodes.bg='sum';
            obj.combinemodes.bg1='sum';
            obj.combinemodes.bg2='sum';
            obj.combinemodes.phot='sum';
            obj.combinemodes.phot1='sum';
            obj.combinemodes.phot2='sum';
            obj.combinemodes.PSFxnm='mean';
            obj.combinemodes.PSFynm='mean';
            obj.combinemodes.PSFxpix='mean';
            obj.combinemodes.PSFypix='mean';
            obj.combinemodes.locprecnm='locp';
            obj.combinemodes.locprecxnm='locp';
            obj.combinemodes.locprecynm='locp';
            obj.combinemodes.locprecznm='locp';            
            
            obj.combinemodes.CRLBphot='min';
            obj.combinemodes.frame='min';
            obj.combinemodes.channel='first';
            obj.combinemodes.logLikelihood='max';
            obj.combinemodes.loglikelihood='max';
            obj.combinemodes.LLrel='max';
            obj.combinemodes.numberInGroup='first';
            obj.combinemodes.groupindex='first';
            obj.combinemodes.filenumber='first';
            obj.combinemodes.clusterdensity='mean';
            obj.combinemodes.cellnumbers='max';
            obj.combinemodes.sitenumbers='max';
        end
        function connect(obj,dx,dt,framef,xf,yf,locpf,varargin)
            fn=fieldnames(obj.locData.loc);
            lm=0;
            for k=1:length(fn)
                lm=max(lm,length(obj.locData.loc.(fn{k})));
            end
            for k=1:length(fn)
                if length(obj.locData.loc.(fn{k}))<lm
                    disp([fn{k} ' too short, padded with zeros']);
                    obj.locData.loc.(fn{k})(lm)=obj.locData.loc.(fn{k})(end);
                end
            end           
            x=double(obj.locData.getloc(xf).(xf));
            y=double(obj.locData.getloc(yf).(yf)); 
            frame=double(obj.locData.getloc(framef).(framef));
            numfields=length(varargin);
            %our algorithm requires localizations sorted by frame, then by
            %x:
            sortmatrix=zeros(length(frame),numfields+2);
            sortmatrix(:,end-1)=frame;
            sortmatrix(:,end)=x;
            for k=1:length(varargin)
                sortmatrix(:,k)=double(obj.locData.getloc(varargin{k}).(varargin{k}));
            end
            sm=size(sortmatrix);
            [sortmatrixsort,indsort]=sortrows(sortmatrix,1:sm(2));
            
            newgroup=false(length(frame)-1,1);
            for k=1:length(varargin)
                if ~any(isnan(diff(sortmatrixsort(:,k))))
                newgroup=newgroup | diff(sortmatrixsort(:,k));
                end
            end

            %later: remove multiple consequtive 1. For sites of length 1.
            %Always remove the first ones.
            fng=find(newgroup);
            if ~isempty(fng)
                if fng(1)==1
                fng(1)=[];
                end
                newgroup(fng-1)=0;
                sortmatrixsort(newgroup,end-1)=0;
            end       
            clear sortmatrix 
            
            % combine either with fixed search radius, or depending on
            % loalization precision
            switch obj.mode
                case 2 %locprec
                    sigmafactor=2.5;
                    locp=max(double(obj.locData.getloc(locpf).(locpf))*sigmafactor,dx(1)); %less than pixelsiz/10 is not resolvable also for HD
                    list=connectsinglesigma(sortmatrixsort(:,end),double(y(indsort)),sortmatrixsort(:,end-1),double(dx(end)),int32(dt),double(locp));
                case 1 %fix
                   maxactive=10000;
                   list=connectsingle2c(sortmatrixsort(:,end),double(y(indsort)),sortmatrixsort(:,end-1),double(dx),int32(dt),int32(maxactive));
            end
            clear  sortmatrixsort
            clear frame

            if list(end)==0
                list(end)=max(list)+1; %FIX connectsingle doesnt assign last loc. Fix later!
            end  
            if list(1)==0
                list(1)=1; %FIX connectsingle doesnt assign last loc. Fix later!
            end   
            
            numbers=1:sm(1);
            indold=numbers(indsort);
            clear numbers
            [~,indback]=sort(indold);
            listback=list(indback);
            clear indback;         
            obj.locData.setloc('groupindex',listback)
            
             %number of locs
            [listsort,indsort2]=sort(list);
            clear list

            if ~isempty(listsort)
                numbergroup=countlocs(double(listsort));
                obj.dllistsort=[diff(listsort);1];
                clear listsort          

                indold2=indold(indsort2);
                clear indsort2 indold
                [~,indback2]=sort(indold2);
                clear indold2
                obj.locData.setloc('numberInGroup',single(numbergroup(indback2)));
                [~,obj.indsortlist]=sort(listback);           
            end 
        end
        function combine(obj,field,combinemode,weights,gweights) %no field etc: group everythign for which we have combinemodes
            if nargin==1 %do all
                fn2=fieldnames(obj.locData.loc);
                if isempty(obj.locData.loc.(fn2{1}))  
                   return
                end

                %XXXXX weights should be 1/s2. Wikipedia.
                % also:if xerr, yerr, zerr: use proper ones for weighing.  
                if isfield(obj.locData.loc,'locprecnm')
                     weights=double(1./(obj.locData.getloc('locprecnm').locprecnm).^2); %w=1/s^2
                elseif isfield(obj.locData.loc,'xnmerr')
                     weights=double(1./(obj.locData.getloc('xnmerr').xnmerr).^2); %w=1/s^2
                else                   
                     weights=double((obj.locData.getloc('phot').phot)); %w=N
                end
               
                if isempty(weights)
                    weights=ones(size(obj.locData.loc.(fn2{1})));  
                else
                    weights(isinf(weights))=1;
                end
                
                list=obj.locData.getloc('groupindex').groupindex;
                gweights=sumcombineind(double(weights),double(list),double(obj.indsortlist));
                for k=1:length(fn2)
                    if isfield(obj.combinemodes,fn2{k})
                        combinemode=obj.combinemodes.(fn2{k});
                    elseif length(fn2{k})>=3 && strcmp(fn2{k}(end-2:end),'err')
                        combinemode='locp';
                    else
                        combinemode='mean';
                    end
                    obj.combine(fn2{k},combinemode,weights,gweights);
                end
            else %group only single field            
                if nargin>2
                    obj.combinemodes.(field)=combinemode;
                elseif  nargin==2
                    if isfield(obj.combinemodes,field)
                        combinemode=obj.combinemodes.(field);
                    elseif length(field)>=3 && strcmp(field(end-2:end),'err')
                        combinemode='locp';
                    else
                        combinemode='mean';
                    end
                end

                vtype=obj.locData.getloc(field).(field)(1);
                if iscell(vtype)
                    return
                end
                v=double(obj.locData.getloc(field).(field));

                list=obj.locData.getloc('groupindex').groupindex;

                if nargin <4
                    weights=ones(size(v),'like',v);              
                end   
                if nargin <5
                          gweights=sumcombineind(double(weights),double(list),double(obj.indsortlist));
                end

                cmode=0;
                switch combinemode
                    case 'sum'
                        v2=v;  
                    case 'mean'
                        v2=v.*weights;                 
                    case 'square' %currently not used
                        v2=v.^2.*weights.^2;  %really weights^2?
                        if nargin <6 
                            gweights2=sumcombineind(double(weights.^2),double(list),double(obj.indsortlist));
                        end
                    case 'min'
                        v2=-v;
                        cmode=1;
                    case 'max'
                        v2=v;
                        cmode=1;
                    case 'first'
                        cmode=2;
                        v2=v;   
                    case 'locp'
                        v2=1./v.^2;  
                end

            %combine modes: <x>, sqrt(<x^2>), first, sum
                indsort=obj.indsortlist;
                if cmode==2 %first
                    vwsort=v2(indsort);
                    vwout2=vwsort(obj.dllistsort>0);
                else 
                    if cmode==0 %sum up: for sum and average
                        vwout=sumcombineind(double(v2),double(list),(indsort));
                    switch combinemode
                        case 'sum'
                            vwout2=vwout;  
                        case 'mean'
                            vwout2=vwout./gweights;                 
                        case 'square'
                            vwout2=sqrt(vwout./gweights2); 
                        case 'locp'
                            vwout2=1./sqrt(vwout);
                    end   
                    elseif cmode==1 %min/max
                        vwout=maxcombine(double(v2),double(list),(indsort));
                        if strcmp(combinemode,'min')
                            vwout2=-vwout;
                        else
                            vwout2=vwout;
                        end
                    end
                end
                obj.locData.grouploc.(field)=cast(vwout2,'like',vtype);              
            end
        end
    end
end
