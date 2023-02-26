function [svalfit, istepfit]=fitstepind(x,istepfit,mfun)
    for k=1:10
        svalfit=stepvalue(x,istepfit,mfun);
        istepfitold=istepfit;
        istepfit=moveind(x,svalfit,istepfit);   
        if all(istepfit==istepfitold)
%             k
            break
        end
    end
end

function istep=moveind(x,sval,istep)
istep=[istep; length(x)];
for k=2:length(istep)-1
    dmin=0;
    dplus=0;
    smin=inf;splus=inf;
    ih=istep(k);

    if abs(x(ih)-sval(k))<abs(x(ih)-sval(k-1)) %ih>1 && abs(x(ih-1)-sval(k))<abs(x(ih-1)-sval(k-1))
        smin=abs(x(ih)-sval(k));
        dmin=-1;
    elseif ih>1 && abs(x(ih-1)-sval(k))<abs(x(ih-1)-sval(k-1)) %ih>1 && abs(x(ih-1)-sval(k))<abs(x(ih-1)-sval(k-1))
        dmin=-2;
        smin=abs(x(ih-1)-sval(k));
    end
    if ih<length(x) && abs(x(ih+1)-sval(k-1))<abs(x(ih+1)-sval(k))% ih<length(x) && abs(x(ih)-sval(k))<abs(x(ih)-sval(k-1))
        dplus=1;
        splus=abs(x(ih+1)-sval(k-1));
        
    elseif ih<length(x)-1 && abs(x(ih+2)-sval(k-1))<abs(x(ih+2)-sval(k))
        splus=abs(x(ih+2)-sval(k-1));
        dplus=2;
    end
     if smin<splus
         istep(k)=max(1,max(istep(k)+dmin, istep(k-1)+1));
     else 
         istep(k)=max(1,min(istep(k)+dplus,istep(k+1)-1));
     end
end
istep=istep(1:end-1);
end