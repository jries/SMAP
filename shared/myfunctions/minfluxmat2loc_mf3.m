function loc=minfluxmat2loc_mf3(jt)
numlocs=size(jt.cfr,2);
locs=single(jt.loc*1e9);
locx=locs(:,1);
locy=locs(:,2);
locz=locs(:,3);

indgood=~isnan(locx);

% numlocs=size(locs,1);
locnums=(1:numlocs)';

    % cols=1:numiter;
    % lxgc=cols.*lxg;
    % lm=max(lxgc,[],2); %last iteration
    % goodind=lm>0;
    % g1=find(goodind);
    % g2=lm(goodind); %take only last iteration
    % indx=sub2ind(size(locs),g1,g2,ones(size(g1)));
    % indy=sub2ind(size(locs),g1,g2,2*ones(size(g1)));
    % indz=sub2ind(size(locs),g1,g2,3*ones(size(g1)));
    % ind2=sub2ind(size(jt.cfr),g1,g2);
    % ind1=goodind;
    % loc.iterations(:,1)=g2;
% end

loc.xnm(:,1)=locx(indgood);
loc.ynm(:,1)=locy(indgood);
znm=locz(indgood);
if any(znm~=0)
    loc.znm(:,1)=znm;
end

loc.frame(:,1)=double(1:length(loc.xnm));
loc.locnum(:,1)=locnums(indgood);
loc.time(:,1)=double(jt.tim(loc.locnum))*1e3;  %from seconds to milliseconds
previous=loc.locnum-1; previous(previous<1)=1;
loc.dt(:,1)=single(jt.tim(loc.locnum)-jt.tim(previous))*1e3;
loc.dcr(:,1)=single(jt.dcr(indgood));
loc.cfr(:,1)=single(jt.cfr(indgood));
loc.vld(:,1)=jt.vld(loc.locnum);
loc.phot(:,1)=single(jt.eco(indgood));
loc.eco(:,1)=single(jt.eco(indgood));
loc.efo(:,1)=single(jt.efo(indgood));
loc.ecc(:,1)=single(jt.ecc(indgood));
loc.efc(:,1)=single(jt.efc(indgood));
loc.itr(:,1)=single(jt.itr(indgood));
loc.sta(:,1)=single(jt.sta(indgood));
loc.tid(:,1)=single(jt.tid(indgood));
loc.gri(:,1)=single(jt.gri(indgood));
loc.thi(:,1)=single(jt.thi(indgood));
loc.sqi(:,1)=single(jt.sqi(indgood));
loc.itr(:,1)=single(jt.itr(indgood));
loc.fbg(:,1)=single(jt.fbg(indgood));
loc.tid(:,1)=single(jt.tid(loc.locnum));
loc.fnl(:,1)=single(jt.fnl(loc.locnum));
loc.bot(:,1)=single(jt.bot(loc.locnum));
if isfield(jt,'lnc')
    locsnc=single(jt.lnc*1e9);
    locsxnc=locsnc(:,1);
    locsync=locsnc(:,2);
    locsznc=locsnc(:,3);
    loc.xncnm(:,1)=locsxnc(indgood);
    loc.yncnm(:,1)=locsync(indgood);
    zncnm=locsznc(indgood);
    if any(zncnm~=0)
        loc.zncnm(:,1)=zncnm;
    end
end
    %determine localization precision
   % ext=single(jt.ext*1e9);
   %      extx=ext(:,:,1);
   %      exty=ext(:,:,2);
   %      extz=ext(:,:,3);
        % loc.extx(:,1)=extx(indgood);
        % loc.exty(:,1)=exty(indgood);
        % extzg=extz(indgood);
        % if any(extzg~=0)
        %     loc.extz(:,1)=extzg;
        % end
        
        % L=median(extx(indgood),2);
        L=100; %xxxx nots sure where to find this
loc.locprecnm(:,1)=L./2./sqrt(2*loc.phot);
loc.LLrel=loc.cfr;
end