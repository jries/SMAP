function [sval,istep]=insertstep(sval,istep,insertind)
%     if insertind>=length(istep)
%         return
%     end
    istep=[istep(1:insertind) ;istep(insertind) ;istep(insertind+1:end)];
    sval=[sval(1:insertind-1) ;(sval(insertind-1)+sval(insertind))/2 ;sval(insertind:end)];  
    istep(insertind)=istep(insertind)-1;
    istep(insertind+1)=istep(insertind+1)+1;
end