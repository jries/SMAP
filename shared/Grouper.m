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
            obj.combinemodes.xnm='mean';
            obj.combinemodes.ynm='mean';
            obj.combinemodes.znm='mean';
            obj.combinemodes.bg='sum';
            obj.combinemodes.bg1='sum';
            obj.combinemodes.bg2='sum';
            obj.combinemodes.phot='sum';
            obj.combinemodes.phot1='sum';
            obj.combinemodes.phot2='sum';
            obj.combinemodes.PSFxnm='square';
            obj.combinemodes.PSFynm='square';
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
            obj.combinemodes.neighbours='sum';
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
%             obj.status('group localizations')
            %here not clear. like this or exchanged?
            
            x=double(obj.locData.getloc(xf).(xf));
            y=double(obj.locData.getloc(yf).(yf));
            
             frame=double(obj.locData.getloc(framef).(framef));
            
            numfields=length(varargin);
            
%             sortmatrix=([frame,x,y]);  %modify:unconnected
%        sorting of y not necessary, as all y are compared.
            sortmatrix=zeros(length(frame),numfields+2);
            sortmatrix(:,end-1)=frame;
            sortmatrix(:,end)=x;
%             sortmatrix(:,end)=y;
            for k=1:length(varargin)
                sortmatrix(:,k)=double(obj.locData.getloc(varargin{k}).(varargin{k}));
%                 sortmatrix=[double(obj.locData.getloc(varargin{k}).(varargin{k})),sortmatrix];
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
            
%             frame(indsort(newgroup))=0;
            
            
            clear sortmatrix 
%             
%             list=connectsingle2c(sortmatrixsort(:,end),double(y(indsort)),sortmatrixsort(:,end-1),double(dx),int32(dt),int32(maxactive));
           
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
%              listm=connectsingle2mat(double(x(indsort)),double(y(indsort)),double(frame(indsort)),double(dx),int32(dt),int32(maxactive));
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
%             numbergroup=zeros(size(listsort));
%             whos listsort
            if ~isempty(listsort)
            numbergroup=countlocs(double(listsort));
            obj.dllistsort=[diff(listsort);1];
            clear listsort
           
%             groupc=listsort(1);
%             ng=0;
%             inds=1;
%             for k=1:length(numbergroup)
%                 if listsort(k)~=groupc
%                     groupc=listsort(k);
%                     numbergroup(inds:k-1)=ng;
%                     ng=0;
%                     inds=k;
%                 end
%                 ng=ng+1;
% %                 numbergroup(k)=ng;
%             end
            
            
            indold2=indold(indsort2);
            clear indsort2 indold
            [~,indback2]=sort(indold2);
            clear indold2
            obj.locData.setloc('numberInGroup',single(numbergroup(indback2)));
            
            [~,obj.indsortlist]=sort(listback);
            
            
            
  
            end
%             obj.status('group localizations done')

            
        end
        function combine(obj,field,combinemode,weights,gweights) %no field etc: group everythign for which we have combinemodes
            
%             obj.status('group: combine fields')
            if nargin==1 %do all
%                 fn=fieldnames(obj.combinemodes);
                fn2=fieldnames(obj.locData.loc);
                if isempty(obj.locData.loc.(fn2{1}))  
                   return
                end
%                 fnall=intersect(fn,fn2);
                %XXXXX weights should be 1/s2. Wikipedia. not really
                % also:if xerr, yerr, zerr: use proper ones for weighing.
                %intuitive...
                weights=1./(obj.locData.getloc('locprecnm').locprecnm);
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
                    else
                        combinemode='mean';
                    end
                    obj.combine(fn2{k},combinemode,weights,gweights);
                end
                
            else
                
            if nargin>2
                obj.combinemodes.(field)=combinemode;
            elseif  nargin==2
                if isfield(obj.combinemodes,'field')
                combinemode=obj.combinemodes.(field);
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
            mode=0;
            switch combinemode
                case 'sum'
                    v2=v;  
                case 'mean'
                    v2=v.*weights;                 
                case 'square'
                    v2=v.^2.*weights.^2;
                    gweights=sumcombineind(double(weights.^2),double(list),double(obj.indsortlist));
%                     weights=weights.^2;
                case 'min'
                    v2=-v;
                    mode=1;
                case 'max'
                    v2=v;
                    mode=1;
                case 'first'
                    mode=2;
                    v2=v;   
                case 'locp'
                    v2=1./v.^2;  
            end
            
            

            %combine modes
            %<x>
            %sqrt(<x^2>)
            %first
            %sum
%             [listsort,indsort]=sort(list);
            indsort=obj.indsortlist;
%             listsort=list(indsort);
%             vwsort=v2(indsort);

            if mode==2
%                 listsort=list(indsort);
                vwsort=v2(indsort);
                
%                 dl=[diff(listsort);1];
                vwout2=vwsort(obj.dllistsort>0);
%                 listsort(end)
%                 sum(dl==1)
            else 
                if mode==0
                    vwout=sumcombineind(double(v2),double(list),(indsort));
%                 vwout=sumcombine(double(vwsort),double(listsort));
%                 gweights=sumcombine(double(weights(indsort)),double(listsort));
                switch combinemode
                    case 'sum'
                        vwout2=vwout;  
                    case 'mean'
%                         gweights=sumcombineind(double(weights),double(list),(indsort));
                        vwout2=vwout./gweights;                 
                    case 'square'
%                         gweights=sumcombineind(double(weights),double(list),(indsort));
                        vwout2=sqrt(vwout./gweights); 
                    case 'locp'
                        vwout2=1./sqrt(vwout);
                end   
                elseif mode==1
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
%             obj.status('group: combine fields done')
        end
    end
end
