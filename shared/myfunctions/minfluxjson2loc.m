function loc=minfluxjson2loc(jt)
valid=[jt.vld];
jtv=jt(valid);
for k=length(jtv):-1:1
    locs=single(jtv(k).itr(end).loc*1e9);
    loc.xnm(k,1)=locs(1);
    loc.ynm(k,1)=locs(2);
    if length(locs)>2
        loc.znm(k,1)=locs(3);
    end
    loc.time(k,1)=single(jtv(k).tim)*1e3;  %from seconds to milliseconds
    loc.iter(k,1)=single(jtv(k).itr(end).itr);
    loc.cfr(k,1)=single(jtv(k).itr(end).cfr);
    loc.dcr(k,1)=single(jtv(k).itr(end).dcr);
    loc.frame(k,1)=k;
end

loc.xnm=loc.xnm-min(loc.xnm);
loc.ynm=loc.ynm-min(loc.ynm);
end