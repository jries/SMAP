function v=getFieldAsVectorInd(p,varargin)
% ls=length(p);
isarray=false;
ind=varargin{end};
if ~isnumeric(ind)
    ind=[];
    fieldnames=varargin(1:end);
else
    fieldnames=varargin(1:end-1);
end

fieldnames=takeapart(fieldnames);
for k=length(p):-1:1
%     try
    if iscell(p)
        ph=p{k};
    else
        ph=p(k);
    end
    if isfield(ph,fieldnames{1}) || isprop(ph,fieldnames{1})
    vh=ph.(fieldnames{1});
        for f=2:length(fieldnames)
           if isempty(vh) || ~(isfield(vh,fieldnames{f}) || isprop(vh,fieldnames{f}))
               vh=NaN;
               break 
           end
            nv=vh.(fieldnames{f});
            vh=nv;
        end
        if ~isempty(ind) %index passed
            vh=vh(ind);
        end
    else
       vh=NaN;
    end
    
%     if isarray||(numel(vh)==1 && (isnumeric(vh)||islogical(vh)))
%         v(k)=vh;
%         isarray=true;
    if isarray||((isnumeric(vh)||islogical(vh)))
        v(k,1:numel(vh))=vh(:);
        isarray=true;
    else
        v{k}=vh;
    end
%     catch err
%     if isarray
%         v(k)=NaN;
%     else
%         v{k}=NaN;
%     end
%     end

end


function fno=takeapart(fnin)
if ~iscell(fnin)
    fnin={fnin};
end
ind=1;
for k=1:length(fnin)
    stop=strfind(fnin{k},'.');
    if isempty(stop)
        fno{ind}=fnin{k};
        ind =ind +1;
    else
        stop=[0 stop length(fnin{k})+1];
        for s=1:length(stop)-1
            fhere=fnin{k}(stop(s)+1:stop(s+1)-1);
            fno{ind}=fhere;
            ind=ind+1;
        end
    end
end