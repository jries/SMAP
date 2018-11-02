function allmd=gethashtable(tab,exclude)
if nargin<2||isempty(exclude)
    exclude={};
end
                
               
k=tab.keys;
ind=1;
allmd=cell(0,2);
while k.hasNext
    kh=k.next;
    if ~any(strncmp(exclude,kh,10))
        try
            v=tab.get(kh);
            if isnumeric(v)
                v=num2str(v);
            elseif ~ischar(v)
                v=[];
            end
%             vs=v.toString;
            if ~isempty(v)
                allmd(ind,1:2)={kh ,v};
                ind=ind+1;
            end
        catch
        end

    end
end