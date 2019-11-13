l=g.locData.getloc({'numberInGroup','filenumber'},'grouping','grouped'); mean(l.numberInGroup)
lt=[];
for k=1:max(l.filenumber)
    indh=l.filenumber==k;
    lt(k)=mean(l.numberInGroup(indh));
end

times=[0.1 0.2 0.3 0.5 0.7 1 2 3 5 7 10];
figure(88);plot(times,lt,'o')

%mean(numberinGroup)=lifetime+1