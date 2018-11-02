function [fitpos,outim,outimnorm,ci]=my2Dgaussfit(image,startp,cas)

%fit par=(x,y,a,bg,sx,sy,r)
if nargin<3
    cas=2;
    show=false;
end

if nargout>1
    show=true;
end
% if ishandle(show)
%     hfig=show;
%     show=true;
% elseif show
%    hf=figure(2);
%    hfig=axes('Parent',hf);
% end
% %normalize by max pixel
% maximage=max(image(:));
% image=image/maximage;

weighted=0;
    
  
    switch cas
        case 1 %only x, y
            fitinit.fitf=@gaussPALMVfix; fitinit.nfitp=4;
        case 2% also sigma
            fitinit.fitf=@gaussPALMVsigma; fitinit.nfitp=5;
        case 3 %sx sy rho
            fitinit.fitf=@gaussPALMVfree; fitinit.nfitp=7;
    end

    %determine startparameters
    s=size(image);
    xi=1:s(2);
    yi=1:s(1);
    gxi=sum(image,1);
    gyi=sum(image,2)';
    
    
    xm=sum(xi.*gxi)/sum(gxi);
    ym=sum(yi.*gyi)/sum(gyi);
    
    gxi(gxi<0)=0;gyi(gyi<0)=0;
    sxm=sqrt(sum((xi-xm).^2.*gxi)/sum(gxi));
    sym=sqrt(sum((yi-ym).^2.*gyi)/sum(gyi));
    if ~isreal(sxm)
        sxm=1;
    end
    if ~isreal(sym)
        sym=1;
    end
    
[Xi,Yi]=meshgrid(xi,yi);
  rsxsy= sum( sum((Xi-xm).*(Yi-ym).*image))/sum(image(:));
  rho=rsxsy/sxm/sym;
    
    bg=min(image(:));
    a=max(image(:))-bg;
%     rho=0;
    
    %startparameter
%     sx=startp(5);
%     sy=startp(6); rho=startp(7);
   

    mIter=100;
    
    %construct matrix and grid
%     fitp.V=[sx^2,rho*sx*sy;rho*sx*sy,sy^2];
    fitp.V=[sxm^2,rho*sxm*sym;rho*sxm*sym,sym^2];
    fitp.Vinv=inv(fitp.V);
    s=size(image);
    x=1:s(1);
    y=1:s(2);
    [fitp.X,fitp.Y]=meshgrid(x,y);
    fitinit.fitp=fitp;
    %handle to function with differences
    fitferr=str2func([func2str(fitinit.fitf) 'err']);
%     %initialize fi



%     oldopts=optimset('lsqnonneg');
       oldopts=optimset('lsqnonlin');
%     options=optimset(oldopts,'TolX',10^-6,'TolFun',10^-6,'Jacobian','on','Algorithm',{'levenberg-marquardt',0.005},'Display','on','Diagnostics','off','MaxIter',mIter);
options=optimset(oldopts,'TolX',10^-6,'TolFun',10^-6,'Jacobian','on','Display','off','Diagnostics','off','MaxIter',mIter);


   
    %initialize mYlsq
    %x,y, A,BG,V1, V2, V4
%     fitinit.startp=double([startp(1),startp(2),100,min(image(:)),fitp.Vinv([1 2 4])]); 
%     fitinit.lsqstruc =mYlsqnonlininit(fitferr,fitinit.startp(1:fitinit.nfitp),[],[], options,fitp,fitp.X); %fitp.X dummy insize of smallframe
fitinit.startp=double([xm,ym,a,bg,fitp.Vinv([1 2 4])]); 
if ~isempty(find(isnan(fitinit.startp)))
    fitinit.startp=zeros(8,1);
end

       startph=fitinit.startp(1:fitinit.nfitp);
%        startph(4)=startp(4);
%        startph(3)=startp(3);
      if weighted
           image=sqrt(image);
       end
       [fitout,resnorm,residual,exitflag,output,lambda,jacobian] =lsqnonlin(fitferr,startph,[],[], options,fitp,image);
          fitpos=fitout2fitpos(fitout,0,0,fitinit.nfitp,fitinit.startp);
           

            if weighted

            chi2=(resnorm)/((s(1)*s(2))-fitinit.nfitp-1); 
            else

                chi2=sum(residual(:).^2./image(:))/((s(1)*s(2))-fitinit.nfitp-1); 
            end
    
           fitpos(9)=chi2;

            fitpos(8)=output.funcCount; 
%             fitpos(3:4)=fitpos(3:4)*maximage;
            ci=nlparci(fitout,residual,'jacobian',jacobian);
            
      

        if show            
            
            im2=fitinit.fitf(startph,fitinit.fitp);
            im3=fitinit.fitf(fitout,fitinit.fitp);
            im4=fitinit.fitf(fitout,fitinit.fitp)-image;
            outim=[image,im2;im3,im4];
            
            outimnorm=[imnorm(image),imnorm(im2);imnorm(im3),imnorm(im4)];
%             imagesc([image,im2;im3,im4],'Parent',hfig)
%             waitforbuttonpress
%          drawnow
        end


function fitpos=fitout2fitpos(fitout,mx,my,nfitp,startp)
fitpos(1)=fitout(2)+mx;
fitpos(2)=fitout(1)+my; %%%%%%%%%%%%%%%%%%%%%%%%%%%% careful
fitpos(3:4)=fitout(3:4);
if nfitp==4
    fitpos(5)=Vinv2sigma(startp(5));
    fitpos(6)=fitpos(5); %sx=xy
    fitpos(7)=0; %rho=0
elseif nfitp==5
   fitpos(5)=Vinv2sigma(fitout(5));
   fitpos(6)=fitpos(5); %sx=xy
   fitpos(7)=0;
elseif nfitp==7
   fitpos(5:7)=Vinv2sigma(fitout(5),fitout(6),fitout(7));
else
    disp('error')
end

function [out,J]=gaussPALMVfixerr(par,fitp,frame)
if nargout==1
out=gaussPALMVfix(par,fitp)-(frame);
else
    [fitframe,J]=gaussPALMVfix(par,fitp);
    out=fitframe-frame;
end
out=out(:);

function [out,J]=gaussPALMVfix(par,fitp)
%matrix fix. free parameters: x0,y0,b,A
x0=par(1);y0=par(2);A=par(3);b=par(4);
xpon=-0.5*(fitp.Vinv(1,1)*(fitp.X-x0).^2+2*fitp.Vinv(1,2)*(fitp.X-x0).*(fitp.Y-y0)+fitp.Vinv(2,2)*(fitp.Y-y0).^2);
FA=exp(xpon);
out=A*FA+b;
if nargout>1

   x1=A*FA.* (fitp.Vinv(1,1)*(fitp.X - x0) +  fitp.Vinv(1,2)*(fitp.Y - y0));
    x2=A* FA.*(fitp.Vinv(1,2)*(fitp.X - x0) + fitp.Vinv(2,2)*(fitp.Y - y0));
   x3= FA;
   x4=ones(length(FA));
   J=[x1(:),x2(:),x3(:),x4(:)];
%    size(J)

end

function [out,J]=gaussPALMVfreeerr(par,fitp,frame)
if nargout==1
out=(frame)-gaussPALMVfree(par,fitp);
else
    [fitframe,J]=gaussPALMVfree(par,fitp);
    out=fitframe-frame;
end

function [out,J]=gaussPALMVfree(par,fitp)
%matrix fix. free parameters: x0,y0,b,A,V11,V12,V22
x0=par(1);y0=par(2);A=par(3);b=par(4);V11=par(5);V12=par(6);V22=par(7);
xpon=-0.5*(V11*(fitp.X-x0).^2+2*V12*(fitp.X-x0).*(fitp.Y-y0)+V22*(fitp.Y-y0).^2);
FA=exp(xpon);
out=A*FA+b;
if nargout>1

   x1=A*FA.* (V11*(fitp.X - x0) +  V12*(fitp.Y - y0));
    x2=A* FA.*(V12*(fitp.X - x0) + V22*(fitp.Y - y0));
   x3= FA;
   x4=ones(length(FA));

   x5=-A*FA.* ((fitp.X - x0).^2)/2;
   x6=-A*FA.* (fitp.X - x0).*(fitp.Y - y0);
   x7=-A*FA.* ((fitp.Y - y0).^2)/2;
   J=[x1(:),x2(:),x3(:),x4(:),x5(:),x6(:),x7(:)];
%    size(J)

end

function [out,J]=gaussPALMVsigmaerr(par,fitp,frame)
if nargout==1
out=gaussPALMVsigma(par,fitp)-(frame);
else
    [fitframe,J]=gaussPALMVsigma(par,fitp);
    out=fitframe-frame;
end
out=out(:);

function [out,J]=gaussPALMVsigma(par,fitp)
%matrix fix. free parameters: x0,y0,b,A
x0=par(1);y0=par(2);A=par(3);b=par(4);Vi1=par(5);
xpon=-0.5*(Vi1*(fitp.X-x0).^2+Vi1*(fitp.Y-y0).^2);
FA=exp(xpon);%no norm
out=A*FA+b;
if nargout>1

   x1=A*FA.* (Vi1*(fitp.X - x0));
    x2=A* FA.*( Vi1*(fitp.Y - y0));
   x3= FA;
   x4=ones(length(FA));
   x5=-A*FA/2.*((fitp.X - x0).^2+(fitp.Y - y0).^2);
   J=[x1(:),x2(:),x3(:),x4(:),x5(:)];
%    size(J)

end

function s= Vinv2sigma(V11,V12,V22)
if nargin==3
    s=real([sqrt(V22./(-V12.^2 + V11.*V22)),  sqrt(V11./(-V12.^2 + V11.*V22)), V12./sqrt(V11.*V22)]);
elseif nargin==1;
    s=real(1/sqrt(V11));
elseif nargin==0;
    s=[];
end

function V= sigma2Vinv(sx,sy,r)
V=[(sx.^2 - r.^2.*sx.^2).^(-1), -(r./(sx.*sy - r.^2.*sx.*sy)), (sy.^2 - r.^2.*sy.^2).^(-1)];


function [out,J]=gaussPALMVfixwerr(par,fitp,frame)
if nargout==1
out=gaussPALMVfixw(par,fitp)-(frame);
else
    [fitframe,J]=gaussPALMVfixw(par,fitp);
    out=fitframe-frame;
end
out=out(:);

function [out,J]=gaussPALMVfixw(par,fitp)
%matrix fix. free parameters: x0,y0,b,A
x0=par(1);y0=par(2);A=par(3);b=par(4);
xpon=-0.5*(fitp.Vinv(1,1)*(fitp.X-x0).^2+2*fitp.Vinv(1,2)*(fitp.X-x0).*(fitp.Y-y0)+fitp.Vinv(2,2)*(fitp.Y-y0).^2);
FA=exp(xpon);
out=sqrt(A*FA+b);
if nargout>1
    x4=0.5./out; %ok
    x3= FA.*x4; %ok
    
   x1=A*x3.* (fitp.Vinv(1,1)*(fitp.X - x0) +  fitp.Vinv(1,2)*(fitp.Y - y0));
    x2=A* x3.*(fitp.Vinv(1,2)*(fitp.X - x0) + fitp.Vinv(2,2)*(fitp.Y - y0));
  
   
   J=[x1(:),x2(:),x3(:),x4(:)];
%    size(J)

end

function [out,J]=gaussPALMVsigmawerr(par,fitp,frame)
if nargout==1
out=gaussPALMVsigmaw(par,fitp)-(frame);
else
    [fitframe,J]=gaussPALMVsigmaw(par,fitp);
    out=fitframe-frame;
end
out=out(:);

function [out,J]=gaussPALMVsigmaw(par,fitp)
%matrix fix. free parameters: x0,y0,b,A
x0=par(1);y0=par(2);A=par(3);b=par(4);Vi1=par(5);
xpon=-0.5*(Vi1*(fitp.X-x0).^2+Vi1*(fitp.Y-y0).^2);
FA=exp(xpon);%no norm
out=sqrt(A*FA+b);
if nargout>1
    x4=0.5./out; %ok
    x3= FA.*x4; %ok
    
    
   x1=A*x3.* (Vi1*(fitp.X - x0));
    x2=A* x3.*( Vi1*(fitp.Y - y0));

   x5=-A*x3/2.*((fitp.X - x0).^2+(fitp.Y - y0).^2);
   J=[x1(:),x2(:),x3(:),x4(:),x5(:)];
%    size(J)

end

function [out,J]=gaussPALMVfreewerr(par,fitp,frame)
if nargout==1
out=(frame)-gaussPALMVfreew(par,fitp);
else
    [fitframe,J]=gaussPALMVfreew(par,fitp);
    out=fitframe-frame;
end

function [out,J]=gaussPALMVfreew(par,fitp)
%matrix fix. free parameters: x0,y0,b,A,V11,V12,V22
x0=par(1);y0=par(2);A=par(3);b=par(4);V11=par(5);V12=par(6);V22=par(7);
xpon=-0.5*(V11*(fitp.X-x0).^2+2*V12*(fitp.X-x0).*(fitp.Y-y0)+V22*(fitp.Y-y0).^2);
FA=exp(xpon);
out=sqrt(A*FA+b);
if nargout>1
    x4=0.5./out; %ok
    x3= FA.*x4; %ok
   x1=A*x3.* (V11*(fitp.X - x0) +  V12*(fitp.Y - y0));
    x2=A* x3.*(V12*(fitp.X - x0) + V22*(fitp.Y - y0));


   x5=-A*x3.* ((fitp.X - x0).^2)/2;
   x6=-A*x3.* (fitp.X - x0).*(fitp.Y - y0);
   x7=-A*x3.* ((fitp.Y - y0).^2)/2;
   J=[x1(:),x2(:),x3(:),x4(:),x5(:),x6(:),x7(:)];
%    size(J)

end

function imout=imnorm(imin)
minim=min(imin(:));maxim=max(imin(:));
imout=(imin-minim)/(maxim-minim);
