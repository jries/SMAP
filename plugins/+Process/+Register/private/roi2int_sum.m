function out=roi2int_sum(roi1,roi2,bg1,bg2,roisize)
sim=size(roi1);
mp=round(sim+1)/2;
dn=min(min(mp)+1,round((roisize+1)/2));


n1=squeeze(sum(sum(roi1(mp(1)-dn:mp(1)+dn,mp(2)-dn:mp(2)+dn,:),1),2));
n2=squeeze(sum(sum(roi2(mp(1)-dn:mp(1)+dn,mp(2)-dn:mp(2)+dn,:),1),2));
b1=squeeze(sum(sum(bg1(mp(1)-dn:mp(1)+dn,mp(2)-dn:mp(2)+dn,:),1),2));
b2=squeeze(sum(sum(bg2(mp(1)-dn:mp(1)+dn,mp(2)-dn:mp(2)+dn,:),1),2));
out.sum_n1=n1-b1;
out.sum_n2=n2-b2;

out.sum_bg1=b1/(2*dn+1)^2;
out.sum_bg2=b2/(2*dn+1)^2;