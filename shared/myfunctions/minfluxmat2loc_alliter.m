function loc=minfluxmat2loc_alliter(jt, onlyvalid,loadall)
numiter=size(jt.cfr,2);
locs=single(jt.loc*1e9);
% if onlyvalid
%     valid=find(jt.vld);
%     indx=sub2ind(size(locs),valid,numiter*ones(size(valid)),ones(size(valid)));
%     indy=sub2ind(size(locs),valid,numiter*ones(size(valid)),2*ones(size(valid)));
%     indz=sub2ind(size(locs),valid,numiter*ones(size(valid)),3*ones(size(valid)));
%     ind2=sub2ind(size(jt.cfr),valid,numiter*ones(size(valid)));
%     ind1=jt.vld; 
% else  
    lx=jt.loc(:,:,1);
    lxg=~isnan(lx);
    cols=1:numiter;
    lxgc=cols.*lxg;
    lm=max(lxgc,[],2); %last iteration
    goodind=lm>0;
    g1=find(goodind);
    g2=lm(goodind); %take only last iteration
    indx=sub2ind(size(locs),g1,g2,ones(size(g1)));
    indy=sub2ind(size(locs),g1,g2,2*ones(size(g1)));
    indz=sub2ind(size(locs),g1,g2,3*ones(size(g1)));
    ind2=sub2ind(size(jt.cfr),g1,g2);
    ind1=goodind;
    
    loc.iterations(:,1)=g2;
% end

loc.xnm(:,1)=locs(indx);
loc.ynm(:,1)=locs(indy);
znm=locs(indz);
if any(znm>0)
    loc.znm(:,1)=locs(indz);
end
loc.time(:,1)=double(jt.tim(ind1))*1e3;  %from seconds to milliseconds
loc.frame(:,1)=double(1:length(loc.xnm));
loc.dcr(:,1)=single(jt.dcr(ind2));
loc.cfr(:,1)=single(jt.cfr(ind2));
loc.vld(:,1)=jt.vld(ind1);
loc.phot(:,1)=single(jt.eco(ind2));
loc.eco(:,1)=single(jt.eco(ind2));
loc.efo(:,1)=single(jt.efo(ind2));
loc.ecc(:,1)=single(jt.ecc(ind2));
loc.efc(:,1)=single(jt.efc(ind2));
if loadall
    
    loc.sta(:,1)=single(jt.sta(ind2));
    if isfield(jt,'fbg')
        loc.fbg(:,1)=single(jt.fbg(ind2));
    end
    loc.tid(:,1)=single(jt.tid(ind1));
    loc.act(:,1)=single(jt.act(ind1));
    loc.sky(:,1)=single(jt.sky(ind1));
    if isfield(jt,'lnc')
        locsnc=single(jt.lnc*1e9);
        loc.xncnm(:,1)=locsnc(indx);
        loc.yncnm(:,1)=locsnc(indy);
        zncnm=locsnc(indz);
        if any(zncnm>0)
            loc.zncnm(:,1)=locsnc(indz);
        end
    end
end

%determine localization precision
extfinal=jt.ext(:,end,:);
L=median(extfinal(:))*1e9;
loc.locprecnm(:,1)=L./2./sqrt(2*loc.phot);
loc.LLrel=loc.cfr;
end