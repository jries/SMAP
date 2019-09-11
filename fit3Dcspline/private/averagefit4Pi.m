function out=averagefit4Pi(vf,varargin)
zwindow=5; %+/- zwin

%varargin: name, vector
%pass on 'shared' variable to then calculate the proper averages
%also allow one variable to be passed on 

 P=double(vf.fit.P);CRLB1=double(vf.fit.CRLB);
%collect fitted parameters
sim=vf.sim;
for k=size(CRLB1,2):-1:1
    Pr(:,k,:)=reshape(P(:,k),[],sim(4));
    Cr(:,k,:)=reshape(CRLB1(:,k),[],sim(4));
end
Pr(:,k+1,:)=reshape(P(:,k+1),[],sim(4)); %iterations, not in crlb

fields={'x','y','N','BG','z','ph'};
shared=vf.fit.shared;
ind=1;
for k=1:length(shared)
    if shared(k)
       l.(fields{k})=(Pr(:,ind,:));  dl.(fields{k})=(Cr(:,ind,:)); 
       ind=ind+1;
    else
        l.(fields{k})=(Pr(:,ind:ind+3,:));  dl.(fields{k})=(Cr(:,ind:ind+3,:)); 
        ind=ind+4;
    end
end
l.ph=mod(l.ph,2*pi);
% xfit=squeeze(Pr(:,1,:));dx=squeeze(Cr(:,1,:));
% yfit=squeeze(Pr(:,2,:));dy=squeeze(Cr(:,2,:));
% Nfit=squeeze(Pr(:,3,:));
% Bg=squeeze(Pr(:,4,:));
% phase=mod(squeeze(Pr(:,6,:)),2*pi);
% zastigf=squeeze(Pr(:,5,:));
%determine average positions in small window


mpz=ceil(size(vf.fit.PSF.Ispline,3)/2+0.5);
droi=(size(vf.imstacksq,1)-1)/2;
zrange=mpz-zwindow:mpz+zwindow;
numbeads=sim(4);


xn=1:sim(1);yn=1:sim(2);zn=1:sim(3);
[Xq,Yq,Zq]=meshgrid(yn,xn,zn);

for k=numbeads:-1:1
    %linear fit to phase, remove jumps
    for c=1:size(l.ph,2)
        phaseh=squeeze(l.ph(zrange,c,k));
        dp=find(diff(phaseh)<-1);
        if ~isempty(dp)
            phaseh(dp+1:end)=phaseh(dp+1:end)+2*pi;
        end
        Xmat=horzcat(zrange(:)-mpz, zrange(:)*0+1);
        linfit=Xmat\(phaseh);   
        out.phase(k,c)=mod(linfit(2),2*pi);
    end
    for c=1:size(l.z,2)
        zastigh=squeeze(l.z(zrange,c,k))-mpz;
        Xmat=horzcat(zastigh(:), zastigh(:)*0+1);
        linfit=Xmat\(zrange(:)-mpz);
        out.z0(k,c)=linfit(2);
    end
%     for c=1:size(l.x,2)
    out.x0(k,:)=squeeze(sum(l.x(zrange,:,k)./dl.x(zrange,:,k),1)./sum(1./dl.x(zrange,:,k),1))-droi+1;
    out.y0(k,:)=squeeze(sum(l.y(zrange,:,k)./dl.y(zrange,:,k),1)./sum(1./dl.y(zrange,:,k),1))-droi+1;
    out.N(k,:)=squeeze(sum(l.N(zrange,:,k)./dl.N(zrange,:,k),1)./sum(1./dl.N(zrange,:,k),1));
    if nargin>1
        for f=1:2:length(varargin)
            field=varargin{f};
            addv=varargin{f+1};
            out.(field)(k,:)=squeeze(nanmean(addv(zrange,k,:),1));
        end
    end
%         out.x0(k,:)=squeeze(sum(xfit(zrange,k)./dx(zrange,k),1)./sum(1./dx(zrange,k),1))-droi+1;
%     end
%     out.y0(k,:)=squeeze(sum(yfit(zrange,k)./dy(zrange,k),1)./sum(1./dy(zrange,k),1))-droi+1;
%     
%     out.N(k,:)=squeeze(mean(Nfit(zrange,k),1));
%     for c=1:size(vf.imstack,5)
%         imh=squeeze(imstack(:,:,:,k,c));
%         xshift=-y0(k,c); %works empirically
%         yshift=-x0(k,c);
%         zshift=0; %shift IAB in z only
% %         zshift=-z0(k);
%         shiftedh=interp3(imh(:,:,:),Xq-xshift,Yq-yshift,Zq-zshift,'cubic',0);
%         imstackaligned(:,:,:,k,c)=shiftedh;
%     end
%     [I,A,B]=make4Pimodel(squeeze(imstackaligned(:,:,:,k,:)),phaseshifts+phasem(k),ph.frequency,1./Nhere(k,:));
%        [I,A,B]=make4Pimodel(squeeze(imstackaligned(:,:,:,k,:)),phaseshifts+phasem(k),ph.frequency,PSF.normf);
%     Is=interp3(I,Xq,Yq,Zq+z0(k),'cubic',0);
%     As=interp3(A,Xq,Yq,Zq+z0(k),'cubic',0);
%     Bs=interp3(B,Xq,Yq,Zq+z0(k),'cubic',0);
%     
%     Aa(:,:,:,k,:)=As;
%     Ba(:,:,:,k,:)=Bs;
%     Ia(:,:,:,k,:)=Is;
end
