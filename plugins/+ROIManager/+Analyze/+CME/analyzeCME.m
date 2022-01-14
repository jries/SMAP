global sitepar
% sitelist=sitepar.sitelist.list;
 direct='D:\Data\CME\rbCHC_autompicked_classified_roi2\';
 mkdir(direct)
     cm=hot(256);
x=0;
y=0;
z=0;
dz=0;
dx=0;
dy=0;
dtheta=0;
dtheta10=0;
ori=0;
use=0;
curv=0;
rfit=0;
class=0;
keep=0;
keeptemp=0;
coverage_area_raw=0;
coverage_area_main=0;
coverage_fraction_raw=0;
coverage_fraction_main=0;
numberoflocs_mask=0;
footprint=0;
newthetafit=0;
newrfit=0;

analysis='theta plots new'
% analysis='save'


switch analysis
    case 'save'
%      for k=1:length(sitelist)


% for k=1:186
for k=1:length(sitelist)
    si=sitelist(k).siteinfo.evaluation.clathrin;
    siten=sitelist(k).number;
    celln=sitelist(k).cell.number;
    filen=sitelist(k).file.number;
    v=sitelist(k).siteinfo.parlist.values;
    lis=['L' num2str(v{1}-1) num2str(v{2}-1) num2str(v{3}-1) num2str(v{4}-1)];
    
    fn=['F' num2str(filen) 'C' num2str(celln) 'S' num2str(siten) lis '.tif'];
    x(k)=si.eval.qz(3)-si.eval.qz(1);
    y(k)=si.eval.qr;
    imallrgb=ind2rgb(uint8(si.all*255+1),cm);
     imwrite(imallrgb, [direct fn]);
%     figure(33)   
%     image(imallrgb)
%     colormap hot
%     title(fn)
%     waitforbuttonpress
     end


case 'x z plots'
%     for k=1:length(sitelist)
    for k=1:179
        si=sitelist(k).siteinfo.evaluation.clathrin;
        saveclathrin(k)=si;
        dx(k)=si.eval.qx(3)-si.eval.qx(1);
        dy(k)=si.eval.qy(3)-si.eval.qy(1);
%         dz(k)=si.eval.qz(3)-si.eval.qz(1);
        dz(k)=si.eval.allqz(end-1)-si.eval.allqz(1);
        
        z(k)=si.eval.qz(3);
        zc(k)=si.eval.qz(2);
        posx(k)=sitelist(k).pos(1);
        posy(k)=sitelist(k).pos(2);
    end
    si.eval.allrange([end-3 3])
    
    figure(90)
    subplot(2,2,1)
    plot(dx,dz,'.')
    xlabel('dx');ylabel('dz');
        subplot(2,2,2)
        
    plot(z,dz,'.')
        xlabel('z');ylabel('dz');
        subplot(2,2,3)
    plot(z,dx,'.')
        xlabel('z');ylabel('dx');
        
        ind_glass=z>0.1;
        
    subplot(2,2,4)
    plot(dx(ind_glass),dz(ind_glass),'.')
    xlabel('dx');ylabel('dz');
        
    zint=zc-min(zc);
    zint=round(zint/max(zint)*254+1);
    max(zint)
    col=jet(256);
    coloval=col(zint,:);
    figure(3)
    scatter(posx,posy,30,coloval)
    colormap jet
    colorbar

case 'theta plots'
%     for k=991:length(sitelist)
    for k=1:length(sitelist)
        si=sitelist(k).siteinfo.evaluation.clathrin;
        saveclathrin(k)=si;
        
        dx(k)=si.eval.allqx(end)-si.eval.allqx(1);
        dy(k)=si.eval.qy(3)-si.eval.qy(1);
%         dz(k)=si.eval.qz(3)-si.eval.qz(1);
        dz(k)=si.eval.allqz(end)-si.eval.allqz(1);
%         allqz(1) is 5%, allqz(end) is 95%f
        dtheta(k)=si.eval.allqtheta(end-1);
        dtheta10(k)=si.eval.allqtheta(2);
        diff(k)=abs(si.eval.allqtheta(end-3))-abs(si.eval.allqtheta(4));
        if ( abs(si.eval.allqtheta(end-3)) < abs(si.eval.allqtheta(4)))
            orient(k)=1;
%             right way around
        else
            orient(k)=0;
%             upside down
        end        
        rfit(k)=si.eval.spherefit(1);
        z(k)=si.eval.allqz(end);
        zc(k)=si.eval.qz(2);
        posx(k)=sitelist(k).pos(1);
        posy(k)=sitelist(k).pos(2);
        class(k)=sitelist(k).siteinfo.parlist.values{3};
        if sitelist(k).siteinfo.parlist.values{1}==1
            keep(k)=1;
        else
            keep(k)=0;
        end
        sphere_segm(k)=2*pi*rfit(k)*(rfit(k)*(1-sin(-dtheta(k))));
        % M =2pi*rh, h=r*(1-sin(-theta))
        % M =2pi*rh,    for dtheta<0: r= dx/2*cos(theta)
        %               for dtheta>0: r=dx/2
        %               

        if dtheta < 0
            sphere_segm_dx_dth(k)=2*pi*(1-sin(-1*dtheta(k)))*((0.5*dx(k)/cos(-1*dtheta(k))).^2);
        else
            sphere_segm_dx_dth(k)=2*pi*(1-sin(-1*dtheta(k)))*((0.5*dx(k)).^2);
        end
            
%         class(k)=sitelist(k).siteinfo.parlist.values{1};
%         if class(k) == 1
%             col(k)='red'
%         elseif class(k) == 2
%             col(k) = 'blue'
%         else
%             col(k) = 'green'
%         end
        
    end
    si.eval.allrange([end-3 3])

    
    
    
    

%     
%     ind_bg_fl=class==1; %red, big flat
%     ind_sm_fl=class==2; 
%     ind_sm_rnd=class==3;
%     ind_sm=class==2 | class==3; %green, both small's
%     ind_late=class==4; %blue, late

    
    
    % e.g. for 2_SKmel2_ClaRFP_2rbCHC-CLC_647_640i100_TIRF_3DA_new_annotation_2apr15_sites
%     ind_bg_fl=class==4; %red, big flat
%     ind_sm_fl=class==2; 
%     ind_sm=class>1&class<4; %green, both small's
%     ind_sm_rnd=class==3;
%     ind_late=class==1; %blue, late
%     
        
% for 1_SKmel2_ClaRFP_2rbCHC-CLC_647_640i100_TIRF_3DA_x15
%     ind_bg_fl=class==1; %not vesicles, red
%     ind_sm=class==2; %green, vesicles
%     ind_late=class==4;
%     
    
% 3 color code

% for new classification: 1-small flat, 2-big flat, 3-slightly curved,
% 4-curvedC, 5-curvedU, 6-omega, 7-closed
    ind_sm_fl=(class==1 & keep);  %green, small    
    ind_bg_fl=(class==2 & keep); %red, big flat
    ind_sm=(class==1 & keep); %green
    ind_slightly=(class==3 & keep);%magenta, slightly curved
    ind_C=(class==4 & keep);%cyan, curved C
    ind_U=(class==5 & keep); %yellow, curved U
    ind_omega=(class==6 & keep); %blue, omega
    ind_closed=(class==7 & keep); %gray, closed
    
    ind_late=((class==4 | class==5 | class==6 | class==7) & keep); %blue, late

    ind_keep=keep==1;
    
    
    
    % GLASS LEVEL - CHANGE HERE
    glass=0.0;
    topcov=-.8;
    ind_glass=(z>glass & keep); %all sites that are so far down, so that they start at 0 or lower
    ind_10=(dtheta10<topcov & keep); %all sites that are sufficiently covered on top
    
    
    
    % many plots of dz, dx, rfit ... uzsing 3 color code
if 0
        
    figure(91)
    clf

    subplot(3,4,1)
    plot(dx(ind_late),dz(ind_late),'.',dx(ind_sm),dz(ind_sm),'.g',dx(ind_bg_fl),dz(ind_bg_fl),'.r')
    xlabel('width (dx)');ylabel('height (dz)');
    
    subplot(3,4,2) 
    plot(z(ind_late),dz(ind_late),'.',z(ind_sm),dz(ind_sm),'.g',z(ind_bg_fl),dz(ind_bg_fl),'.r')
    xlabel('z start');ylabel('height (dz)');
   
    subplot(3,4,3)
    plot(z(ind_late),dx(ind_late),'.',z(ind_sm),dx(ind_sm),'.g',z(ind_bg_fl),dx(ind_bg_fl),'.r')
    xlabel('z start');ylabel('width (dx)');
        
    subplot(3,4,5)
    plot(dx(ind_late),dtheta(ind_late),'.',dx(ind_sm),dtheta(ind_sm),'.g',dx(ind_bg_fl),dtheta(ind_bg_fl),'.r')
    xlabel('width (dx)');ylabel('dtheta');
         
    subplot(3,4,6)
    plot(dx(ind_late & ind_glass),dtheta(ind_late & ind_glass),'.',dx(ind_sm & ind_glass),dtheta(ind_sm & ind_glass),'.g',dx(ind_bg_fl & ind_glass),dtheta(ind_bg_fl & ind_glass),'.r')
    xlabel(sprintf('dx with z>%.2f', glass));ylabel('dtheta');
 
    subplot(3,4,4)
    plot(rfit(ind_late),dx(ind_late),'.',rfit(ind_sm),dx(ind_sm),'.g',rfit(ind_bg_fl),dx(ind_bg_fl),'.r')
    xlabel('rfit');ylabel('dx');
    
    subplot(3,4,7)
    plot(rfit(ind_late),dtheta(ind_late),'.',rfit(ind_sm),dtheta(ind_sm),'.g',rfit(ind_bg_fl),dtheta(ind_bg_fl),'.r')
    xlabel('rfit');ylabel('dtheta');
    
    subplot(3,4,8)
    plot(rfit(ind_late & ind_glass),dtheta(ind_late & ind_glass),'.',rfit(ind_sm & ind_glass),dtheta(ind_sm & ind_glass),'.g',rfit(ind_bg_fl & ind_glass),dtheta(ind_bg_fl & ind_glass),'.r')
    xlabel(sprintf('rfit with z>%.2f', glass));ylabel('dtheta');

    subplot(3,4,9)
    plot(rfit(ind_late),dtheta10(ind_late),'.',rfit(ind_sm),dtheta10(ind_sm),'.g',rfit(ind_bg_fl),dtheta10(ind_bg_fl),'.r')
    xlabel('rfit');ylabel('dtheta10');
    legend('late','small','big_ flat');
    
    subplot(3,4,10)
    plot(dx(ind_late & ind_glass & ind_10),dtheta(ind_late & ind_glass & ind_10),'.',dx(ind_sm & ind_glass & ind_10),dtheta(ind_sm & ind_glass & ind_10),'.g',dx(ind_bg_fl & ind_glass & ind_10),dtheta(ind_bg_fl & ind_glass & ind_10),'.r')
    xlabel(strcat(sprintf('dx with z>%.2f and d10<', glass),sprintf('%.2f',topcov)));ylabel('dtheta');    
    
    subplot(3,4,12)
    plot(rfit(ind_late & ind_glass & ind_10),dtheta(ind_late & ind_glass & ind_10),'.',rfit(ind_sm & ind_glass & ind_10),dtheta(ind_sm & ind_glass & ind_10),'.g',rfit(ind_bg_fl & ind_glass & ind_10),dtheta(ind_bg_fl & ind_glass & ind_10),'.r')
    xlabel(strcat(sprintf('rfit with z>%.2f and d10<', glass),sprintf('%.2f',topcov)));ylabel('dtheta');    
    
end
    
   

    % many plots of dz, dx, rfit ..... using complex 7 color code
if 0
    figure(93)
    clf;
    
    subplot(3,4,1)
    plot(dx(ind_sm_fl),dz(ind_sm_fl),'.g');
    hold on
    plot(dx(ind_bg_fl),dz(ind_bg_fl),'.r');
    plot(dx(ind_slightly),dz(ind_slightly),'.m');
    plot(dx(ind_C),dz(ind_C),'.c');
    plot(dx(ind_U),dz(ind_U),'.y');
    plot(dx(ind_omega),dz(ind_omega),'.');
    plot(dx(ind_closed),dz(ind_closed),'.k');
    hold off
    xlabel('width (dx)');ylabel('height (dz)');
    
    subplot(3,4,2) 
    plot(z(ind_late),dz(ind_late),'.',z(ind_sm),dz(ind_sm),'.g',z(ind_bg_fl),dz(ind_bg_fl),'.r')
    xlabel('z start');ylabel('height (dz)');
   
    subplot(3,4,3)
    plot(z(ind_late),dx(ind_late),'.',z(ind_sm),dx(ind_sm),'.g',z(ind_bg_fl),dx(ind_bg_fl),'.r')
    xlabel('z start');ylabel('width (dx)');
    
    subplot(3,4,4)
    plot(rfit(ind_late),dx(ind_late),'.',rfit(ind_sm),dx(ind_sm),'.g',rfit(ind_bg_fl),dx(ind_bg_fl),'.r')
    xlabel('rfit');ylabel('dx');
         
%     subplot(3,4,5)
%     plot(dx(ind_sm_fl),dtheta(ind_sm_fl),'.g');
%     hold on
%     plot(dx(ind_bg_fl),dtheta(ind_bg_fl),'.r');
%     plot(dx(ind_slightly),dtheta(ind_slightly),'.m');
%     plot(dx(ind_C),dtheta(ind_C),'.c');
%     plot(dx(ind_U),dtheta(ind_U),'.y');
%     plot(dx(ind_omega),dtheta(ind_omega),'.');
%     plot(dx(ind_closed),dtheta(ind_closed),'.k');
%     hold off
%     xlabel('width (dx)');ylabel('dtheta');
             
    subplot(3,4,5)
    plot(dx(ind_sm_fl),dtheta(ind_sm_fl),'.g');
    hold on
    plot(dx(ind_bg_fl | ind_slightly),dtheta(ind_bg_fl | ind_slightly),'.r');
%     plot(dx(ind_slightly),dtheta(ind_slightly),'.m');
    plot(dx(ind_C),dtheta(ind_C),'.c');
    plot(dx(ind_U),dtheta(ind_U),'.y');
    plot(dx(ind_omega),dtheta(ind_omega),'.');
    plot(dx(ind_closed),dtheta(ind_closed),'.','MarkerEdgeColor',[.7 .7 .7]);
    hold off
    xlabel('width (dx)');ylabel('dtheta');
    
    subplot(3,4,6)
    plot(dx(ind_late & ind_glass),dtheta(ind_late & ind_glass),'.',dx(ind_sm & ind_glass),dtheta(ind_sm & ind_glass),'.g',dx(ind_bg_fl & ind_glass),dtheta(ind_bg_fl & ind_glass),'.r')
    xlabel(sprintf('dx with z>%.2f', glass));ylabel('dtheta');
 
    subplot(3,4,7)
    plot(rfit(ind_late),dtheta(ind_late),'.',rfit(ind_sm),dtheta(ind_sm),'.g',rfit(ind_bg_fl),dtheta(ind_bg_fl),'.r')
    xlabel('rfit');ylabel('dtheta');
    
    subplot(3,4,8)
    plot(rfit(ind_late & ind_glass),dtheta(ind_late & ind_glass),'.',rfit(ind_sm & ind_glass),dtheta(ind_sm & ind_glass),'.g',rfit(ind_bg_fl & ind_glass),dtheta(ind_bg_fl & ind_glass),'.r')
    xlabel(sprintf('rfit with z>%.2f', glass));ylabel('dtheta');

    subplot(3,4,9)
    plot(rfit(ind_late),dtheta10(ind_late),'.',rfit(ind_sm),dtheta10(ind_sm),'.g',rfit(ind_bg_fl),dtheta10(ind_bg_fl),'.r')
    xlabel('rfit');ylabel('dtheta10');
    legend('late','small','big_ flat');
    
    subplot(3,4,10)
    plot(dx(ind_late & ind_glass & ind_10),dtheta(ind_late & ind_glass & ind_10),'.',dx(ind_sm & ind_glass & ind_10),dtheta(ind_sm & ind_glass & ind_10),'.g',dx(ind_bg_fl & ind_glass & ind_10),dtheta(ind_bg_fl & ind_glass & ind_10),'.r')
    xlabel(strcat(sprintf('dx with z>%.2f and d10<', glass),sprintf('%.2f',topcov)));ylabel('dtheta');    
    
    subplot(3,4,12)
    plot(rfit(ind_late & ind_glass & ind_10),dtheta(ind_late & ind_glass & ind_10),'.',rfit(ind_sm & ind_glass & ind_10),dtheta(ind_sm & ind_glass & ind_10),'.g',rfit(ind_bg_fl & ind_glass & ind_10),dtheta(ind_bg_fl & ind_glass & ind_10),'.r')
    xlabel(strcat(sprintf('rfit with z>%.2f and d10<', glass),sprintf('%.2f',topcov)));ylabel('dtheta');    
end    
    


% magnified plots of dtheta vs rfit and vs dx, using color code for 7 time
% windows

if 0    
    figure(94)
    clf
    
     subplot(1,2,1)
    plot(dx(ind_sm_fl),dtheta(ind_sm_fl),'.g','MarkerSize',20);
    hold on
    plot(dx(ind_bg_fl | ind_slightly),dtheta(ind_bg_fl | ind_slightly),'.r','MarkerSize',20);
%     plot(dx(ind_slightly),dtheta(ind_slightly),'.m');
    plot(dx(ind_C),dtheta(ind_C),'.c','MarkerSize',20);
    plot(dx(ind_U),dtheta(ind_U),'.y','MarkerSize',20);
    plot(dx(ind_omega),dtheta(ind_omega),'.','MarkerSize',20);
    plot(dx(ind_closed),dtheta(ind_closed),'.','MarkerEdgeColor',[.7 .7 .7],'MarkerSize',20);
    hold off
    xlabel('width (dx)');ylabel('dtheta');
    legend('small flat','big flat','curved C','curved U','curved Omega','closed');
    
     subplot(1,2,2)
    plot(rfit(ind_sm_fl),dtheta(ind_sm_fl),'.g','MarkerSize',20);
    hold on
    plot(rfit(ind_bg_fl | ind_slightly),dtheta(ind_bg_fl | ind_slightly),'.r','MarkerSize',20);
%     plot(rfit(ind_slightly),dtheta(ind_slightly),'.m');
    plot(rfit(ind_C),dtheta(ind_C),'.c','MarkerSize',20);
    plot(rfit(ind_U),dtheta(ind_U),'.y','MarkerSize',20);
    plot(rfit(ind_omega),dtheta(ind_omega),'.','MarkerSize',20);
    plot(rfit(ind_closed),dtheta(ind_closed),'.','MarkerEdgeColor',[.7 .7 .7],'MarkerSize',20);
    hold off
    xlabel('rfit');ylabel('dtheta');
    
end

% calculate means, stds, quantiles for dx,dtheta,rfit in time windows

if 1
    
    dx1=dx(ind_sm_fl);
    dx2=dx(ind_bg_fl | ind_slightly);
    dx3=dx(ind_C);
    dx4=dx(ind_U);
    dx5=dx(ind_omega);
    dx6=dx(ind_closed);
    
    dth1=dtheta(ind_sm_fl);
    dth2=dtheta(ind_bg_fl | ind_slightly);
    dth3=dtheta(ind_C);
    dth4=dtheta(ind_U);
    dth5=dtheta(ind_omega);
    dth6=dtheta(ind_closed);
    
    rfit1=rfit(ind_sm_fl);
    rfit2=rfit(ind_bg_fl | ind_slightly);
    rfit3=rfit(ind_C);
    rfit4=rfit(ind_U);
    rfit5=rfit(ind_omega);
    rfit6=rfit(ind_closed);
    
    dx1m=mean(dx1); dth1m=mean(dth1);
    dx2m=mean(dx2); dth2m=mean(dth2);
    dx3m=mean(dx3); dth3m=mean(dth3);
    dx4m=mean(dx4); dth4m=mean(dth4);
    dx5m=mean(dx5); dth5m=mean(dth5);
    dx6m=mean(dx6); dth6m=mean(dth6);
    
    dx1s=std(dx1); dth1s=std(dth1);
    dx2s=std(dx2); dth2s=std(dth2);
    dx3s=std(dx3); dth3s=std(dth3);
    dx4s=std(dx4); dth4s=std(dth4);
    dx5s=std(dx5); dth5s=std(dth5);
    dx6s=std(dx6); dth6s=std(dth6);
    
    dr1m=mean(rfit1);
    dr2m=mean(rfit2);
    dr3m=mean(rfit3);
    dr4m=mean(rfit4);
    dr5m=mean(rfit5);
    dr6m=mean(rfit6);
    
    dr1s=std(rfit1);
    dr2s=std(rfit2);
    dr2sem=std(rfit2)/sqrt(length(rfit2));
    dr3s=std(rfit3);
    dr4s=std(rfit4);
    dr5s=std(rfit5);
    dr6s=std(rfit6);
        
    dx1q=myquantile(dx1,[0.05,0.5,0.95]);
    dx2q=myquantile(dx2,[0.05,0.5,0.95]);
    dx3q=myquantile(dx3,[0.05,0.5,0.95]);
    dx4q=myquantile(dx4,[0.05,0.5,0.95]);
    dx5q=myquantile(dx5,[0.05,0.5,0.95]);
    dx6q=myquantile(dx6,[0.05,0.5,0.95]);
    
    dth1q=myquantile(dth1,[0.05,0.5,0.95]);
    dth2q=myquantile(dth2,[0.05,0.5,0.95]);
    dth3q=myquantile(dth3,[0.05,0.5,0.95]);
    dth4q=myquantile(dth4,[0.05,0.5,0.95]);
    dth5q=myquantile(dth5,[0.05,0.5,0.95]);
    dth6q=myquantile(dth6,[0.05,0.5,0.95]);
    
    dr1q=myquantile(rfit1,[0.05,0.5,0.95]);
    dr2q=myquantile(rfit2,[0.05,0.5,0.95]);
    dr3q=myquantile(rfit3,[0.05,0.5,0.95]);
    dr4q=myquantile(rfit4,[0.05,0.5,0.95]);
    dr5q=myquantile(rfit5,[0.05,0.5,0.95]);
    dr6q=myquantile(rfit6,[0.05,0.5,0.95]);
end




% magnified plots of dtheta vs rfit and vs dx MEANS+-STD. using color code
% for 7 time windows

if 0    
    figure(99)
    clf
        
    subplot(1,2,1)
    H=errorbarxy(dx1m,dth1m,dx1s,dth1s,{'gs','g','g'});
    hold on
    H=errorbarxy(dx2m,dth2m,dx2s,dth2s,{'rs','r','r'});
    H=errorbarxy(dx3m,dth3m,dx3s,dth3s,{'cs','c','c'});
    H=errorbarxy(dx4m,dth4m,dx4s,dth4s,{'ys','y','y'});
    H=errorbarxy(dx5m,dth5m,dx5s,dth5s,{'bs','b','b'});
    H=errorbarxy(dx6m,dth6m,dx6s,dth6s,{'ks','k','k'});
    hold off
    xlabel('mean(dx)+-sd');ylabel('mean(dtheta)+-sd');
%     xlabel('median(dx)+-sd');ylabel('median(dtheta)+-sd');
    title('dtheta vs dx - mean+-sd');
    subplot(1,2,2)
    H=errorbarxy(dr1m,dth1m,dr1s,dth1s,{'gs','g','g'});
    hold on
    H=errorbarxy(dr2m,dth2m,dr2s,dth2s,{'rs','r','r'});
    H=errorbarxy(dr3m,dth3m,dr3s,dth3s,{'cs','c','c'});
    H=errorbarxy(dr4m,dth4m,dr4s,dth4s,{'ys','y','y'});
    H=errorbarxy(dr5m,dth5m,dr5s,dth5s,{'bs','b','b'});
    H=errorbarxy(dr6m,dth6m,dr6s,dth6s,{'ks','k','k'});
    hold off
    xlabel('mean(rfit)+-sd');ylabel('mean(dtheta)+-sd');
%     xlabel('median(rfit)+-sd');ylabel('median(dtheta)+-sd');
    title('dtheta vs rfit - mean+-sd');
    
%     
%      subplot(1,2,2)
%     plot(mean(rfit(ind_sm_fl)),mean(dth1),'.g','MarkerSize',20);
%     hold on
%     plot(mean(rfit(ind_bg_fl | ind_slightly)),mean(dtheta(ind_bg_fl | ind_slightly)),'.r','MarkerSize',20);
% %     plot(rfit(ind_slightly),dtheta(ind_slightly),'.m');
%     plot(mean(rfit(ind_C)),mean(dtheta(ind_C)),'.c','MarkerSize',20);
%     plot(mean(rfit(ind_U)),mean(dtheta(ind_U)),'.y','MarkerSize',20);
%     plot(mean(rfit(ind_omega)),mean(dtheta(ind_omega)),'.','MarkerSize',20);
%     plot(mean(rfit(ind_closed)),mean(dtheta(ind_closed)),'.','MarkerEdgeColor',[.7 .7 .7],'MarkerSize',20);
%     hold off
%     xlabel('rfit');ylabel('dtheta');
%     
end

if 1

    figure(100)
    clf
    subplot(1,2,1)
    H=ploterr(dx1m,dth1m,{dx1q(1),dx1q(end)},{dth1q(1),dth1q(end)},'g.'); set(H(1:3),'LineWidth',1);
%     rectangle('Position',[0,1,0.1,0.5],'Curvature',[1,1],'EdgeColor','g');
%     rectangle('Position',[0.2,-0.5,0.1,0.5],'Curvature',[1,1],'EdgeColor','g');
    hold on
    H=ploterr(dx2m,dth2m,{dx2q(1),dx2q(end)},{dth2q(1),dth2q(end)},'r.'); set(H(1:3),'LineWidth',1);
    H=ploterr(dx3m,dth3m,{dx3q(1),dx3q(end)},{dth3q(1),dth3q(end)},'c.'); set(H(1:3),'LineWidth',1);
    H=ploterr(dx4m,dth4m,{dx4q(1),dx4q(end)},{dth4q(1),dth4q(end)},'y.'); set(H(1:3),'LineWidth',1);
    H=ploterr(dx5m,dth5m,{dx5q(1),dx5q(end)},{dth5q(1),dth5q(end)},'b.'); set(H(1:3),'LineWidth',1);
    H=ploterr(dx6m,dth6m,{dx6q(1),dx6q(end)},{dth6q(1),dth6q(end)},'k.'); set(H(1:3),'LineWidth',1);
    hold off
    xlabel('mean(dx) +- 5/95 percentile');ylabel('mean(dtheta) +- 5/95 percentile');
    title('dtheta vs dx - mean +- 5/95 percentile');
    
    subplot(1,2,2)
%     H=ploterr(dr1m,dth1m,{dr1q(1),dr1q(end)},{dth1q(1),dth1q(end)},'g.'); set(H(1:3),'LineWidth',1);
    rectangle('Position',[0,1,0.1,0.5],'Curvature',[1,1],'EdgeColor','g');
    rectangle('Position',[0.2,-0.5,0.1,0.5],'Curvature',[1,1],'EdgeColor','g');
    hold on
    H=ploterr(dr2m,dth2m,{dr2q(1),dr2q(end)},{dth2q(1),dth2q(end)},'r.'); set(H(1:3),'LineWidth',1);
    H=ploterr(dr3m,dth3m,{dr3q(1),dr3q(end)},{dth3q(1),dth3q(end)},'c.'); set(H(1:3),'LineWidth',1);
    H=ploterr(dr4m,dth4m,{dr4q(1),dr4q(end)},{dth4q(1),dth4q(end)},'y.'); set(H(1:3),'LineWidth',1);
    H=ploterr(dr5m,dth5m,{dr5q(1),dr5q(end)},{dth5q(1),dth5q(end)},'b.'); set(H(1:3),'LineWidth',1);
    H=ploterr(dr6m,dth6m,{dr6q(1),dr6q(end)},{dth6q(1),dth6q(end)},'k.'); set(H(1:3),'LineWidth',1);
    hold off
    xlabel('mean(rfit) +- 5/95 percentile');ylabel('mean(dtheta) +- 5/95 percentile');
    title('dtheta vs rfit - mean +- 5/95 percentile');
    
    

end

% not color coded graph of dtheta/height vs dx/rfit
if 0
    figure (98)
    clf
    subplot(1,3,1)
    plot(dx(ind_keep),dz(ind_keep),'.','MarkerSize',10);
    xlabel('width');ylabel('height');
    subplot(1,3,2)
    plot(dx(ind_keep),dtheta(ind_keep),'.','MarkerSize',10);
    xlabel('width');ylabel('dtheta');
    subplot(1,3,3)
    plot(rfit(ind_keep),dtheta(ind_keep),'.','MarkerSize',10);
    xlabel('radius');ylabel('dtheta');
end
    
% box plots of rfit, dtheta and surface calculated from fit

if 0
    figure(92)
    clf
    stages={'small (flat)','big flat','slightly curved','curved C','curved U','curved Omega','closed'};
    subplot(3,2,1)
    notBoxPlot(sphere_segm(ind_sm_fl),1);
    notBoxPlot(sphere_segm(ind_bg_fl),2);
    notBoxPlot(sphere_segm(ind_slightly),3);
    notBoxPlot(sphere_segm(ind_C),4);
    notBoxPlot(sphere_segm(ind_U),5);
    notBoxPlot(sphere_segm(ind_omega),6);
    notBoxPlot(sphere_segm(ind_closed),7);
    ylabel('calculated surface /um2');
    title('all sites');
    set(gca,'XLim',[0.1 7.9],'XTick',[1:7],'XTickLabel',stages);
    subplot(3,2,2)
    notBoxPlot(sphere_segm(ind_sm_fl & ind_glass & ind_10),1);
    notBoxPlot(sphere_segm(ind_bg_fl & ind_glass & ind_10),2);
    notBoxPlot(sphere_segm(ind_slightly & ind_glass & ind_10),3);
    notBoxPlot(sphere_segm(ind_C & ind_glass & ind_10),4);
    notBoxPlot(sphere_segm(ind_U & ind_glass & ind_10),5);
    notBoxPlot(sphere_segm(ind_omega & ind_glass & ind_10),6);
    notBoxPlot(sphere_segm(ind_closed & ind_glass & ind_10),7);
    ylabel('calculated surface /um2');
    title(strcat(sprintf('only sites close to glass and covered on top (z>%.2f and d10<', glass),sprintf('%.2f',topcov)));
    set(gca,'XLim',[0.1 7.9],'XTick',[1:7],'XTickLabel',stages);
    subplot(3,2,3)
    notBoxPlot(dx(ind_sm_fl),1);
    notBoxPlot(dx(ind_bg_fl),2);
    notBoxPlot(dx(ind_slightly),3);
    notBoxPlot(dx(ind_C),4);
    notBoxPlot(dx(ind_U),5);
    notBoxPlot(dx(ind_omega),6);
    notBoxPlot(dx(ind_closed),7);
    ylabel('width /um (=dx)');
    title('all sites');
    set(gca,'XLim',[0.1 7.9],'XTick',[1:7],'XTickLabel',stages);
    subplot(3,2,4)
%     notBoxPlot(dx((ind_sm_fl | ind_bg_fl) & ind_glass & ind_10),1);
    notBoxPlot(dx(ind_sm_fl & ind_glass & ind_10),1);
    notBoxPlot(dx(ind_bg_fl & ind_glass & ind_10),2);
    notBoxPlot(dx(ind_slightly & ind_glass & ind_10),3);
    notBoxPlot(dx(ind_C & ind_glass & ind_10),4);
    notBoxPlot(dx(ind_U & ind_glass & ind_10),5);
    notBoxPlot(dx(ind_omega & ind_glass & ind_10),6);
    notBoxPlot(dx(ind_closed & ind_glass & ind_10),7);
    ylabel('width /um (=dx)');
    title(strcat(sprintf('only sites close to glass and covered on top (z>%.2f and d10<', glass),sprintf('%.2f',topcov)));
    set(gca,'XLim',[0.1 7.9],'XTick',[1:7],'XTickLabel',stages);
    subplot(3,2,5)
    notBoxPlot(dtheta(ind_sm_fl),1);
    notBoxPlot(dtheta(ind_bg_fl),2);
    notBoxPlot(dtheta(ind_slightly),3);
    notBoxPlot(dtheta(ind_C),4);
    notBoxPlot(dtheta(ind_U),5);
    notBoxPlot(dtheta(ind_omega),6);
    notBoxPlot(dtheta(ind_closed),7);
    ylabel('how far down (=dtheta)');
    title('all sites');
    set(gca,'XLim',[0.1 7.9],'XTick',[1:7],'XTickLabel',stages);
    subplot(3,2,6)
%     notBoxPlot(dtheta((ind_sm_fl | ind_bg_fl) & ind_glass & ind_10),1);
    notBoxPlot(dtheta(ind_sm_fl & ind_glass & ind_10),1);
    notBoxPlot(dtheta(ind_bg_fl & ind_glass & ind_10),2);
    notBoxPlot(dtheta(ind_slightly & ind_glass & ind_10),3);
    notBoxPlot(dtheta(ind_C & ind_glass & ind_10),4);
    notBoxPlot(dtheta(ind_U & ind_glass & ind_10),5);
    notBoxPlot(dtheta(ind_omega & ind_glass & ind_10),6);
    notBoxPlot(dtheta(ind_closed & ind_glass & ind_10),7);
    ylabel('how far down (=dtheta)');
    title(strcat(sprintf('only sites close to glass and covered on top (z>%.2f and d10<', glass),sprintf('%.2f',topcov)));
    set(gca,'XLim',[0.1 7.9],'XTick',[1:7],'XTickLabel',stages);
end    
    
    

% box plots of surface calculated from rfit vs calculated from dz and dx

if 0
    figure(95)
    clf
    stages={'small (flat)','big flat','slightly curved','curved C','curved U','curved Omega','closed'};
    subplot(2,2,1)
    notBoxPlot(sphere_segm(ind_sm_fl),1);
    notBoxPlot(sphere_segm(ind_bg_fl),2);
    notBoxPlot(sphere_segm(ind_slightly),3);
    notBoxPlot(sphere_segm(ind_C),4);
    notBoxPlot(sphere_segm(ind_U),5);
    notBoxPlot(sphere_segm(ind_omega),6);
    notBoxPlot(sphere_segm(ind_closed),7);
    ylabel('calculated surface /um2');
    title('all sites');
    set(gca,'XLim',[0.1 7.9],'XTick',[1:7],'XTickLabel',stages);
    
   
    subplot(2,2,2)
    notBoxPlot(sphere_segm(ind_sm_fl & ind_glass & ind_10),1);
    notBoxPlot(sphere_segm(ind_bg_fl & ind_glass & ind_10),2);
    notBoxPlot(sphere_segm(ind_slightly & ind_glass & ind_10),3);
    notBoxPlot(sphere_segm(ind_C & ind_glass & ind_10),4);
    notBoxPlot(sphere_segm(ind_U & ind_glass & ind_10),5);
    notBoxPlot(sphere_segm(ind_omega & ind_glass & ind_10),6);
    notBoxPlot(sphere_segm(ind_closed & ind_glass & ind_10),7);
    ylabel('calculated surface /um2');
    title(strcat(sprintf('only sites close to glass and covered on top (z>%.2f and d10<', glass),sprintf('%.2f',topcov)));
    set(gca,'XLim',[0.1 7.9],'XTick',[1:7],'XTickLabel',stages);
    
    subplot(2,2,3)
    notBoxPlot(sphere_segm_dx_dth(ind_sm_fl),1);
    notBoxPlot(sphere_segm_dx_dth(ind_bg_fl),2);
    notBoxPlot(sphere_segm_dx_dth(ind_slightly),3);
    notBoxPlot(sphere_segm_dx_dth(ind_C),4);
    notBoxPlot(sphere_segm_dx_dth(ind_U),5);
    notBoxPlot(sphere_segm_dx_dth(ind_omega),6);
    notBoxPlot(sphere_segm_dx_dth(ind_closed),7);
    ylabel('calculated surface /um2');
    title('calculated from dx and dtheta - all sites');
    set(gca,'XLim',[0.1 7.9],'XTick',[1:7],'XTickLabel',stages);
    
    subplot(2,2,4)
    notBoxPlot(sphere_segm_dx_dth(ind_sm_fl & ind_glass & ind_10),1);
    notBoxPlot(sphere_segm_dx_dth(ind_bg_fl & ind_glass & ind_10),2);
    notBoxPlot(sphere_segm_dx_dth(ind_slightly & ind_glass & ind_10),3);
    notBoxPlot(sphere_segm_dx_dth(ind_C & ind_glass & ind_10),4);
    notBoxPlot(sphere_segm_dx_dth(ind_U & ind_glass & ind_10),5);
    notBoxPlot(sphere_segm_dx_dth(ind_omega & ind_glass & ind_10),6);
    notBoxPlot(sphere_segm_dx_dth(ind_closed & ind_glass & ind_10),7);
    ylabel('calculated surface /um2');
    title(strcat(sprintf('calculated from dx and dtheta - only sites close to glass and covered on top (z>%.2f and d10<', glass),sprintf('%.2f',topcov)));
    set(gca,'XLim',[0.1 7.9],'XTick',[1:7],'XTickLabel',stages);
    
end

% correlation between different calculated surfaces

if 0
    figure(96)
    clf
    stages={'small (flat)','big flat','slightly curved','curved C','curved U','curved Omega','closed'};
    subplot(1,2,1)
    plot(sphere_segm(ind_keep),sphere_segm_dx_dth(ind_keep),'.')
end

















case 'theta plots new'

    
    global se;
         for k=1:length(se.sites)

        si=se.sites(k).evaluation.CME3DDSpherefit;

        dx(k)=si.allqx(end)-si.allqx(1);
        dy(k)=si.qy(3)-si.qy(1);
%         dz(k)=si.eval.qz(3)-si.eval.qz(1);
        dz(k)=si.allqz(end)-si.allqz(1);
%         allqz(1) is 5%, allqz(end) is 95%
        dtheta(k)=si.allqtheta(end-1);
        dtheta10(k)=si.allqtheta(2);
        diff(k)=abs(si.allqtheta(end-3))-abs(si.allqtheta(4));
        if ( abs(si.allqtheta(end-3)) < abs(si.allqtheta(4)))
            orient(k)=1;
%             right way around
        else
            orient(k)=0;
%             upside down
        end        
        rfit(k)=si.spherefit(1);
        curv(k)=1/(si.spherefit(1));
%         ori(k)=(pi-(-1*(dtheta(k))+(pi/2)))*(180/pi);
        ori(k)=(dtheta(k)+(pi/2))*(180/pi);
        z(k)=si.allqz(end);
        zc(k)=si.qz(2);

        posx(k)=se.sites(k).pos(1);
        posy(k)=se.sites(k).pos(2);
        class(k)=se.sites(k).annotation.list3.value;
        fitbad(k)=se.sites(k).annotation.list4.value==2;
        
%         for new analysis: Nov 24th, new surface and area analysis
        coverage_area_raw(k)=si.map3D.coverageArea;
        coverage_area_main(k)=si.map3D.mainArea;
        coverage_fraction_raw(k)=si.map3D.coverageFraction;
        coverage_fraction_main(k)=si.map3D.mainFraction;
        numberoflocs_mask(k)=si.numberOfLocsin3DMask;
        footprint(k)=si.map2D.areanm;
        newthetafit(k)=(si.fitcoverage_thetabottom+pi/2)*(180/pi);
        newrfit(k)=si.map3D.rSphere;
        newcurv(k)=1/newrfit(k);
        

        

    
    

        
% if (se.sites(k).annotation.list1.value < 3) && (se.sites(k).annotation.list2.value==1) && (se.sites(k).annotation.list4.value==1) 
%             keeptemp(k)=1;
%         else
%             keeptemp(k)=0;
% end
%         keep=logical(keeptemp);
        
        
                
%         sphere_segm(k)=2*pi*rfit(k)*(rfit(k)*(1-sin(-dtheta(k))));
        sphere_segm(k)=2*pi*rfit(k)*(rfit(k)*(1-cos((pi/2)+dtheta(k))));
        % M =2pi*rh, h=r*(1-sin(-theta))
        % M =2pi*rh,    for dtheta<0: r= dx/2*cos(theta)
        %               for dtheta>0: r=dx/2
        %               

        if dtheta < 0
            sphere_segm_dx_dth(k)=2*pi*(1-sin(-1*dtheta(k)))*((0.5*dx(k)/cos(-1*dtheta(k))).^2);
        else
            sphere_segm_dx_dth(k)=2*pi*(1-sin(-1*dtheta(k)))*((0.5*dx(k)).^2);
        end
            
%         class(k)=sitelist(k).siteinfo.parlist.values{1};
%         if class(k) == 1
%             col(k)='red'
%         elseif class(k) == 2
%             col(k) = 'blue'
%         else
%             col(k) = 'green'
%         end
        
     end

     %only sites with "use" in ROImanager Apr13 2017
        keep=getFieldAsVector(se.sites,'annotation','use');
      %  keep=logical(ones(1,length(se.sites)));
        
     % for new classification: 1-small flat, 2-big flat, 3-slightly curved,
% 4-curvedC, 5-curvedU, 6-omega, 7-closed
    ind_sm_fl=(class==1 & keep);  %green, small    
    ind_bg_fl=(class==2 & keep); %red, big flat
%     ind_sm=(class==1 & keep); %green
%     ind_slightly=(class==3 & keep);%magenta, slightly curved
    ind_C=(class==3 & keep);%cyan, curved C
    ind_U=(class==4 & keep); %yellow, curved U
    ind_omega=(class==5 & keep); %blue, omega
    ind_closed=(class==6 & keep); %gray, closed
    
    ind_late=((class==3 | class==4 | class==5 | class==6) & keep); %blue, late

    ind_keep=keep==1;
    
    
    
    % GLASS LEVEL - CHANGE HERE
    glass=0.0;
    topcov=-.8;
    ind_glass=(z>glass & keep); %all sites that are so far down, so that they start at 0 or lower
    ind_10=(dtheta10<topcov & keep); %all sites that are sufficiently covered on top
    
    
    
    % many plots of dz, dx, rfit ... uzsing 3 color code
if 0
        
    figure(91)
    clf

    subplot(3,4,1)
    plot(dx(ind_late),dz(ind_late),'.',dx(ind_sm_fl),dz(ind_sm_fl),'.g',dx(ind_bg_fl),dz(ind_bg_fl),'.r')
    xlabel('width (dx)');ylabel('height (dz)');
    
    subplot(3,4,2) 
    plot(z(ind_late),dz(ind_late),'.',z(ind_sm_fl),dz(ind_sm_fl),'.g',z(ind_bg_fl),dz(ind_bg_fl),'.r')
    xlabel('z start');ylabel('height (dz)');
   
    subplot(3,4,3)
    plot(z(ind_late),dx(ind_late),'.',z(ind_sm_fl),dx(ind_sm_fl),'.g',z(ind_bg_fl),dx(ind_bg_fl),'.r')
    xlabel('z start');ylabel('width (dx)');
        
    subplot(3,4,5)
    plot(dx(ind_late),dtheta(ind_late),'.',dx(ind_sm_fl),dtheta(ind_sm_fl),'.g',dx(ind_bg_fl),dtheta(ind_bg_fl),'.r')
    xlabel('width (dx)');ylabel('dtheta');
         
    subplot(3,4,6)
    plot(dx(ind_late & ind_glass),dtheta(ind_late & ind_glass),'.',dx(ind_sm_fl & ind_glass),dtheta(ind_sm_fl & ind_glass),'.g',dx(ind_bg_fl & ind_glass),dtheta(ind_bg_fl & ind_glass),'.r')
    xlabel(sprintf('dx with z>%.2f', glass));ylabel('dtheta');
 
    subplot(3,4,4)
    plot(rfit(ind_late),dx(ind_late),'.',rfit(ind_sm_fl),dx(ind_sm_fl),'.g',rfit(ind_bg_fl),dx(ind_bg_fl),'.r')
    xlabel('rfit');ylabel('dx');
    
    subplot(3,4,7)
    plot(rfit(ind_late),dtheta(ind_late),'.',rfit(ind_sm_fl),dtheta(ind_sm_fl),'.g',rfit(ind_bg_fl),dtheta(ind_bg_fl),'.r')
    xlabel('rfit');ylabel('dtheta');
    
    subplot(3,4,8)
    plot(rfit(ind_late & ind_glass),dtheta(ind_late & ind_glass),'.',rfit(ind_sm_fl & ind_glass),dtheta(ind_sm_fl & ind_glass),'.g',rfit(ind_bg_fl & ind_glass),dtheta(ind_bg_fl & ind_glass),'.r')
    xlabel(sprintf('rfit with z>%.2f', glass));ylabel('dtheta');

    subplot(3,4,9)
    plot(rfit(ind_late),dtheta10(ind_late),'.',rfit(ind_sm_fl),dtheta10(ind_sm_fl),'.g',rfit(ind_bg_fl),dtheta10(ind_bg_fl),'.r')
    xlabel('rfit');ylabel('dtheta10');
    legend('late','small','big_ flat');
    
    subplot(3,4,10)
    plot(dx(ind_late & ind_glass & ind_10),dtheta(ind_late & ind_glass & ind_10),'.',dx(ind_sm_fl & ind_glass & ind_10),dtheta(ind_sm_fl & ind_glass & ind_10),'.g',dx(ind_bg_fl & ind_glass & ind_10),dtheta(ind_bg_fl & ind_glass & ind_10),'.r')
    xlabel(strcat(sprintf('dx with z>%.2f and d10<', glass),sprintf('%.2f',topcov)));ylabel('dtheta');    
    
    subplot(3,4,12)
    plot(rfit(ind_late & ind_glass & ind_10),dtheta(ind_late & ind_glass & ind_10),'.',rfit(ind_sm_fl & ind_glass & ind_10),dtheta(ind_sm_fl & ind_glass & ind_10),'.g',rfit(ind_bg_fl & ind_glass & ind_10),dtheta(ind_bg_fl & ind_glass & ind_10),'.r')
    xlabel(strcat(sprintf('rfit with z>%.2f and d10<', glass),sprintf('%.2f',topcov)));ylabel('dtheta');    
    
end
    
   

    % many plots of dz, dx, rfit ..... using complex 7 color code
if 0
    figure(93)
    clf;
    
    subplot(3,4,1)
    plot(dx(ind_sm_fl),dz(ind_sm_fl),'.g');
    hold on
    plot(dx(ind_bg_fl),dz(ind_bg_fl),'.r');
%     plot(dx(ind_slightly),dz(ind_slightly),'.m');
    plot(dx(ind_C),dz(ind_C),'.c');
    plot(dx(ind_U),dz(ind_U),'.y');
    plot(dx(ind_omega),dz(ind_omega),'.');
    plot(dx(ind_closed),dz(ind_closed),'.k');
    hold off
    xlabel('width (dx)');ylabel('height (dz)');
    
    subplot(3,4,2) 
    plot(z(ind_late),dz(ind_late),'.',z(ind_sm_fl),dz(ind_sm_fl),'.g',z(ind_bg_fl),dz(ind_bg_fl),'.r')
    xlabel('z start');ylabel('height (dz)');
   
    subplot(3,4,3)
    plot(z(ind_late),dx(ind_late),'.',z(ind_sm_fl),dx(ind_sm_fl),'.g',z(ind_bg_fl),dx(ind_bg_fl),'.r')
    xlabel('z start');ylabel('width (dx)');
    
    subplot(3,4,4)
    plot(rfit(ind_late),dx(ind_late),'.',rfit(ind_sm_fl),dx(ind_sm_fl),'.g',rfit(ind_bg_fl),dx(ind_bg_fl),'.r')
    xlabel('rfit');ylabel('dx');
         
%     subplot(3,4,5)
%     plot(dx(ind_sm_fl),dtheta(ind_sm_fl),'.g');
%     hold on
%     plot(dx(ind_bg_fl),dtheta(ind_bg_fl),'.r');
%     plot(dx(ind_slightly),dtheta(ind_slightly),'.m');
%     plot(dx(ind_C),dtheta(ind_C),'.c');
%     plot(dx(ind_U),dtheta(ind_U),'.y');
%     plot(dx(ind_omega),dtheta(ind_omega),'.');
%     plot(dx(ind_closed),dtheta(ind_closed),'.k');
%     hold off
%     xlabel('width (dx)');ylabel('dtheta');
             
    subplot(3,4,5)
    plot(dx(ind_sm_fl),dtheta(ind_sm_fl),'.g');
    hold on
    plot(dx(ind_bg_fl),dtheta(ind_bg_fl),'.r');
%     plot(dx(ind_slightly),dtheta(ind_slightly),'.m');
    plot(dx(ind_C),dtheta(ind_C),'.c');
    plot(dx(ind_U),dtheta(ind_U),'.y');
    plot(dx(ind_omega),dtheta(ind_omega),'.');
    plot(dx(ind_closed),dtheta(ind_closed),'.','MarkerEdgeColor',[.7 .7 .7]);
    hold off
    xlabel('width (dx)');ylabel('dtheta');
    
    subplot(3,4,6)
    plot(dx(ind_late & ind_glass),dtheta(ind_late & ind_glass),'.',dx(ind_sm_fl & ind_glass),dtheta(ind_sm_fl & ind_glass),'.g',dx(ind_bg_fl & ind_glass),dtheta(ind_bg_fl & ind_glass),'.r')
    xlabel(sprintf('dx with z>%.2f', glass));ylabel('dtheta');
 
    subplot(3,4,7)
    plot(rfit(ind_late),dtheta(ind_late),'.',rfit(ind_sm_fl),dtheta(ind_sm_fl),'.g',rfit(ind_bg_fl),dtheta(ind_bg_fl),'.r')
    xlabel('rfit');ylabel('dtheta');
    
    subplot(3,4,8)
    plot(rfit(ind_late & ind_glass),dtheta(ind_late & ind_glass),'.',rfit(ind_sm_fl & ind_glass),dtheta(ind_sm_fl & ind_glass),'.g',rfit(ind_bg_fl & ind_glass),dtheta(ind_bg_fl & ind_glass),'.r')
    xlabel(sprintf('rfit with z>%.2f', glass));ylabel('dtheta');

    subplot(3,4,9)
    plot(rfit(ind_late),dtheta10(ind_late),'.',rfit(ind_sm_fl),dtheta10(ind_sm_fl),'.g',rfit(ind_bg_fl),dtheta10(ind_bg_fl),'.r')
    xlabel('rfit');ylabel('dtheta10');
    legend('late','small','big_ flat');
    
    subplot(3,4,10)
    plot(dx(ind_late & ind_glass & ind_10),dtheta(ind_late & ind_glass & ind_10),'.',dx(ind_sm_fl & ind_glass & ind_10),dtheta(ind_sm_fl & ind_glass & ind_10),'.g',dx(ind_bg_fl & ind_glass & ind_10),dtheta(ind_bg_fl & ind_glass & ind_10),'.r')
    xlabel(strcat(sprintf('dx with z>%.2f and d10<', glass),sprintf('%.2f',topcov)));ylabel('dtheta');    
    
    subplot(3,4,11)
    plot(ori(ind_sm_fl),dz(ind_sm_fl),'.g');
    hold on
    plot(ori(ind_bg_fl),dz(ind_bg_fl),'.r');
%     plot(ori(ind_slightly),dz(ind_slightly),'.m');
    plot(ori(ind_C),dz(ind_C),'.c');
    plot(ori(ind_U),dz(ind_U),'.y');
    plot(ori(ind_omega),dz(ind_omega),'.');
    plot(ori(ind_closed),dz(ind_closed),'.k');
    hold off
    xlabel('ori-dtheta');ylabel('height (dz)');
    
    
    
    subplot(3,4,12)
    plot(rfit(ind_late & ind_glass & ind_10),dtheta(ind_late & ind_glass & ind_10),'.',rfit(ind_sm_fl & ind_glass & ind_10),dtheta(ind_sm_fl & ind_glass & ind_10),'.g',rfit(ind_bg_fl & ind_glass & ind_10),dtheta(ind_bg_fl & ind_glass & ind_10),'.r')
    xlabel(strcat(sprintf('rfit with z>%.2f and d10<', glass),sprintf('%.2f',topcov)));ylabel('dtheta');    
end    
    


% magnified plots of dtheta vs rfit and vs dx, using color code for 7 time
% windows
% PLUS theta vs r-fit as Ori Avinaom defined it
% PLUS curvature vs coverage

if 1    
%     figure(94)
%     clf
%     
%      subplot(1,2,1)
%     plot(dx(ind_sm_fl),dtheta(ind_sm_fl),'.g','MarkerSize',20);
%     hold on
%     plot(dx(ind_bg_fl),dtheta(ind_bg_fl),'.r','MarkerSize',20);
% %     plot(dx(ind_slightly),dtheta(ind_slightly),'.m');
%     plot(dx(ind_C),dtheta(ind_C),'.c','MarkerSize',20);
%     plot(dx(ind_U),dtheta(ind_U),'.y','MarkerSize',20);
%     plot(dx(ind_omega),dtheta(ind_omega),'.','MarkerSize',20);
%     plot(dx(ind_closed),dtheta(ind_closed),'.','MarkerEdgeColor',[.7 .7 .7],'MarkerSize',20);
%     hold off
%     xlabel('width (dx)');ylabel('dtheta');
%     legend('small flat','big flat','curved C','curved U','curved Omega','closed');
%     
%      subplot(1,2,2)
%     plot(rfit(ind_sm_fl),dtheta(ind_sm_fl),'.g','MarkerSize',20);
%     hold on
%     plot(rfit(ind_bg_fl),dtheta(ind_bg_fl),'.r','MarkerSize',20);
% %     plot(rfit(ind_slightly),dtheta(ind_slightly),'.m');
%     plot(rfit(ind_C),dtheta(ind_C),'.c','MarkerSize',20);
%     plot(rfit(ind_U),dtheta(ind_U),'.y','MarkerSize',20);
%     plot(rfit(ind_omega),dtheta(ind_omega),'.','MarkerSize',20);
%     plot(rfit(ind_closed),dtheta(ind_closed),'.','MarkerEdgeColor',[.7 .7 .7],'MarkerSize',20);
%     hold off
%     xlabel('rfit');ylabel('dtheta');
% no    
    figure (999)
    clf
%     subplot (2,2,1)
%     plot(curv(keep),dtheta(keep),'.b','MarkerSize',20);
%     xlabel('curvature');ylabel('coverage');
    
    subplot (2,2,1)
    plot(coverage_fraction_raw(keep),newcurv(keep),'.b','MarkerSize',10);
    xlabel('cov_frct_raw');ylabel('curvature_new');
    
    subplot(2,2,2)
    plot(newthetafit(keep),newrfit(keep),'.b','MarkerSize',10);
    xlabel('theta new fit');ylabel('radius new fit');
    
    subplot(2,2,3)
    histogram(newthetafit(keep),50)
    xlabel('theta new fit');
    
    
    
    
end

    
% surface calculated from fit or calculated from dx and dz

if 0
 figure (996)
    clf
%     subplot(1,2,1)
    plot(ori(keep),sphere_segm(keep),'.b','MarkerSize',20);
    xlabel('theta ori');ylabel('surface from fit');
%     subplot(1,2,2)
%     plot(ori(ind_sm_fl),sphere_segm_dx_dth(ind_sm_fl),'.b','MarkerSize',20);
%     xlabel('theta ori');ylabel('surface from dx, dz');

end

% new analysis nov24

if 1
    figure(997)
    clf
    subplot(4,3,1)
    plot(ori(keep),coverage_area_raw(keep),'.');
    xlabel('theta ori');ylabel('surface raw new');
    subplot(4,3,2)
    plot(ori(keep),coverage_area_main(keep),'.');
    xlabel('theta ori');ylabel('surface main, only bottom hole');
    subplot(4,3,3)
    plot(sphere_segm(keep),coverage_area_raw(keep),'.');
    xlabel('sphere segm');ylabel('surface raw');
    subplot(4,3,4)
    plot(ori(keep),numberoflocs_mask(keep),'.');
    xlabel('theta ori');ylabel('number of locs');
    subplot(4,3,5)
    plot(numberoflocs_mask(keep),coverage_area_raw(keep),'.');
    xlabel('number of locs');ylabel('surface raw');
    subplot(4,3,6)
%     plot(coverage_area_raw(keep),coverage_area_main(keep),'.');
%     xlabel('surface raw');ylabel('surface main');
    plot(ori(keep),coverage_fraction_main(keep),'.');
    xlabel('theta ori');ylabel('coverage fraction main');
    subplot(4,3,7)
    hist(footprint(keep));
    xlabel('footprint');
    subplot(4,3,8)
    plot(coverage_fraction_main(keep),footprint(keep),'.');
    xlabel('coverage fraction main new');ylabel('footprint');
    subplot(4,3,9)
    plot(footprint(keep),coverage_area_raw(keep),'.');
    xlabel('footprint');ylabel('surface raw');
    subplot(4,3,10)
    plot(ori(keep),newthetafit(keep),'.');
    xlabel('theta ori');ylabel('theta new fit');
    subplot(4,3,11)
    plot(coverage_fraction_raw(keep),coverage_area_raw(keep),'.');
    xlabel('coverage fraction raw');ylabel('surface area raw');
    
    
    
end

% coverage_area_raw(k)=si.map3D.coverageArea;
% coverage_area_main(k)=si.map3D.mainArea;
% coverage_fraction_raw(k)=si.map3D.coverageFraction;
% coverage_fraction_main(k)=si.map3D.mainFraction;
% numberoflocs_mask(k)=si.numberOfLocsin3DMask;
% footprint(k)=si.map2D.areanm;
% newthetafit=0;
% newrfit=0;


 if 0
    figure(195)
    clf
    stages={'small (flat)','big flat','slightly curved','curved C','curved U','curved Omega','closed'};
    subplot(1,2,1)
    notBoxPlot(sphere_segm(ind_sm_fl),1);
    notBoxPlot(sphere_segm(ind_bg_fl),2);
    notBoxPlot(sphere_segm(ind_slightly),3);
    notBoxPlot(sphere_segm(ind_C),4);
    notBoxPlot(sphere_segm(ind_U),5);
    notBoxPlot(sphere_segm(ind_omega),6);
    notBoxPlot(sphere_segm(ind_closed),7);
    ylabel('calculated surface /um2');
    title('all sites');
    set(gca,'XLim',[0.1 7.9],'XTick',[1:7],'XTickLabel',stages);
    
       
    subplot(1,2,2)
    notBoxPlot(sphere_segm_dx_dth(ind_sm_fl),1);
    notBoxPlot(sphere_segm_dx_dth(ind_bg_fl),2);
    notBoxPlot(sphere_segm_dx_dth(ind_slightly),3);
    notBoxPlot(sphere_segm_dx_dth(ind_C),4);
    notBoxPlot(sphere_segm_dx_dth(ind_U),5);
    notBoxPlot(sphere_segm_dx_dth(ind_omega),6);
    notBoxPlot(sphere_segm_dx_dth(ind_closed),7);
    ylabel('calculated surface /um2');
    title('calculated from dx and dtheta - all sites');
    set(gca,'XLim',[0.1 7.9],'XTick',[1:7],'XTickLabel',stages);
    
 end
    
case 'theta plots new saved'

    
    global se;
         for k=1:length(infoCME)

        si=infoCME(k).clathrin;

        dx(k)=si.allqx(end)-si.allqx(1);
        dy(k)=si.qy(3)-si.qy(1);
%         dz(k)=si.eval.qz(3)-si.eval.qz(1);
        dz(k)=si.allqz(end)-si.allqz(1);
%         allqz(1) is 5%, allqz(end) is 95%
        dtheta(k)=si.allqtheta(end-1);
        dtheta10(k)=si.allqtheta(2);
        diff(k)=abs(si.allqtheta(end-3))-abs(si.allqtheta(4));
        if ( abs(si.allqtheta(end-3)) < abs(si.allqtheta(4)))
            orient(k)=1;
%             right way around
        else
            orient(k)=0;
%             upside down
        end        
        rfit(k)=si.spherefit(1)*1000;
        curv(k)=1/(si.spherefit(1));
%         ori(k)=(pi-(-1*(dtheta(k))+(pi/2)))*(180/pi);
        ori(k)=(dtheta(k)+(pi/2))*(180/pi);
        z(k)=si.allqz(end);
        zc(k)=si.qz(2);

%         posx(k)=se.sites(k).pos(1);
%         posy(k)=se.sites(k).pos(2);
        class(k)=infoCME(k).parlist.values{3};

if infoCME(k).parlist.values{1}==1
            keep(k)=1;
        else
            keep(k)=0;
        end

%         sphere_segm(k)=2*pi*rfit(k)*(rfit(k)*(1-sin(-dtheta(k))));
        sphere_segm(k)=2*pi*rfit(k)*(rfit(k)*(1-cos((pi/2)+dtheta(k))));
        % M =2pi*rh, h=r*(1-sin(-theta))
        % M =2pi*rh,    for dtheta<0: r= dx/2*cos(theta)
        %               for dtheta>0: r=dx/2
        %               

        if dtheta < 0
            sphere_segm_dx_dth(k)=2*pi*(1-sin(-1*dtheta(k)))*((0.5*dx(k)/cos(-1*dtheta(k))).^2);
        else
            sphere_segm_dx_dth(k)=2*pi*(1-sin(-1*dtheta(k)))*((0.5*dx(k)).^2);
        end
            
%         class(k)=sitelist(k).siteinfo.parlist.values{1};
%         if class(k) == 1
%             col(k)='red'
%         elseif class(k) == 2
%             col(k) = 'blue'
%         else
%             col(k) = 'green'
%         end
        
     end

     
     % for new classification: 1-small flat, 2-big flat, 3-slightly curved,
% 4-curvedC, 5-curvedU, 6-omega, 7-closed
    ind_sm_fl=(class==1 & keep);  %green, small    
    ind_bg_fl=(class==2 & keep); %red, big flat
    ind_sm=(class==1 & keep); %green
    ind_slightly=(class==3 & keep);%magenta, slightly curved
    ind_C=(class==4 & keep);%cyan, curved C
    ind_U=(class==5 & keep); %yellow, curved U
    ind_omega=(class==6 & keep); %blue, omega
    ind_closed=(class==7 & keep); %gray, closed
    
    ind_late=((class==4 | class==5 | class==6 | class==7) & keep); %blue, late

    ind_keep=keep==1;
    
    
    
    % GLASS LEVEL - CHANGE HERE
    glass=0.0;
    topcov=-.8;
    ind_glass=(z>glass & keep); %all sites that are so far down, so that they start at 0 or lower
    ind_10=(dtheta10<topcov & keep); %all sites that are sufficiently covered on top
    
    
    
    % many plots of dz, dx, rfit ... uzsing 3 color code
if 0
        
    figure(91)
    clf

    subplot(3,4,1)
    plot(dx(ind_late),dz(ind_late),'.',dx(ind_sm),dz(ind_sm),'.g',dx(ind_bg_fl),dz(ind_bg_fl),'.r')
    xlabel('width (dx)');ylabel('height (dz)');
    
    subplot(3,4,2) 
    plot(z(ind_late),dz(ind_late),'.',z(ind_sm),dz(ind_sm),'.g',z(ind_bg_fl),dz(ind_bg_fl),'.r')
    xlabel('z start');ylabel('height (dz)');
   
    subplot(3,4,3)
    plot(z(ind_late),dx(ind_late),'.',z(ind_sm),dx(ind_sm),'.g',z(ind_bg_fl),dx(ind_bg_fl),'.r')
    xlabel('z start');ylabel('width (dx)');
        
    subplot(3,4,5)
    plot(dx(ind_late),dtheta(ind_late),'.',dx(ind_sm),dtheta(ind_sm),'.g',dx(ind_bg_fl),dtheta(ind_bg_fl),'.r')
    xlabel('width (dx)');ylabel('dtheta');
         
    subplot(3,4,6)
    plot(dx(ind_late & ind_glass),dtheta(ind_late & ind_glass),'.',dx(ind_sm & ind_glass),dtheta(ind_sm & ind_glass),'.g',dx(ind_bg_fl & ind_glass),dtheta(ind_bg_fl & ind_glass),'.r')
    xlabel(sprintf('dx with z>%.2f', glass));ylabel('dtheta');
 
    subplot(3,4,4)
    plot(rfit(ind_late),dx(ind_late),'.',rfit(ind_sm),dx(ind_sm),'.g',rfit(ind_bg_fl),dx(ind_bg_fl),'.r')
    xlabel('rfit');ylabel('dx');
    
    subplot(3,4,7)
    plot(rfit(ind_late),dtheta(ind_late),'.',rfit(ind_sm),dtheta(ind_sm),'.g',rfit(ind_bg_fl),dtheta(ind_bg_fl),'.r')
    xlabel('rfit');ylabel('dtheta');
    
    subplot(3,4,8)
    plot(rfit(ind_late & ind_glass),dtheta(ind_late & ind_glass),'.',rfit(ind_sm & ind_glass),dtheta(ind_sm & ind_glass),'.g',rfit(ind_bg_fl & ind_glass),dtheta(ind_bg_fl & ind_glass),'.r')
    xlabel(sprintf('rfit with z>%.2f', glass));ylabel('dtheta');

    subplot(3,4,9)
    plot(rfit(ind_late),dtheta10(ind_late),'.',rfit(ind_sm),dtheta10(ind_sm),'.g',rfit(ind_bg_fl),dtheta10(ind_bg_fl),'.r')
    xlabel('rfit');ylabel('dtheta10');
    legend('late','small','big_ flat');
    
    subplot(3,4,10)
    plot(dx(ind_late & ind_glass & ind_10),dtheta(ind_late & ind_glass & ind_10),'.',dx(ind_sm & ind_glass & ind_10),dtheta(ind_sm & ind_glass & ind_10),'.g',dx(ind_bg_fl & ind_glass & ind_10),dtheta(ind_bg_fl & ind_glass & ind_10),'.r')
    xlabel(strcat(sprintf('dx with z>%.2f and d10<', glass),sprintf('%.2f',topcov)));ylabel('dtheta');    
    
        
    
    subplot(3,4,12)
    plot(rfit(ind_late & ind_glass & ind_10),dtheta(ind_late & ind_glass & ind_10),'.',rfit(ind_sm & ind_glass & ind_10),dtheta(ind_sm & ind_glass & ind_10),'.g',rfit(ind_bg_fl & ind_glass & ind_10),dtheta(ind_bg_fl & ind_glass & ind_10),'.r')
    xlabel(strcat(sprintf('rfit with z>%.2f and d10<', glass),sprintf('%.2f',topcov)));ylabel('dtheta');    
    
end
    
   

    % many plots of dz, dx, rfit ..... using complex 7 color code
if 1
    figure(1193)
    clf;
    
    subplot(3,4,1)
    plot(dx(ind_sm_fl),dz(ind_sm_fl),'.g');
    hold on
    plot(dx(ind_bg_fl | ind_slightly),dz(ind_bg_fl | ind_slightly),'.r');
%     plot(dx(ind_slightly),dz(ind_slightly),'.m');
    plot(dx(ind_C),dz(ind_C),'.c');
    plot(dx(ind_U),dz(ind_U),'.y');
    plot(dx(ind_omega),dz(ind_omega),'.');
    plot(dx(ind_closed),dz(ind_closed),'.k');
    hold off
    xlabel('width (dx)');ylabel('height (dz)');
    
    subplot(3,4,2) 
    plot(z(ind_late),dz(ind_late),'.',z(ind_sm),dz(ind_sm),'.g',z(ind_bg_fl),dz(ind_bg_fl),'.r')
    xlabel('z start');ylabel('height (dz)');
   
    subplot(3,4,3)
    plot(z(ind_late),dx(ind_late),'.',z(ind_sm),dx(ind_sm),'.g',z(ind_bg_fl),dx(ind_bg_fl),'.r')
    xlabel('z start');ylabel('width (dx)');
    
    subplot(3,4,4)
    plot(rfit(ind_late),dx(ind_late),'.',rfit(ind_sm),dx(ind_sm),'.g',rfit(ind_bg_fl),dx(ind_bg_fl),'.r')
    xlabel('rfit');ylabel('dx');
         
%     subplot(3,4,5)
%     plot(dx(ind_sm_fl),dtheta(ind_sm_fl),'.g');
%     hold on
%     plot(dx(ind_bg_fl),dtheta(ind_bg_fl),'.r');
%     plot(dx(ind_slightly),dtheta(ind_slightly),'.m');
%     plot(dx(ind_C),dtheta(ind_C),'.c');
%     plot(dx(ind_U),dtheta(ind_U),'.y');
%     plot(dx(ind_omega),dtheta(ind_omega),'.');
%     plot(dx(ind_closed),dtheta(ind_closed),'.k');
%     hold off
%     xlabel('width (dx)');ylabel('dtheta');
             
    subplot(3,4,5)
    plot(dx(ind_sm_fl),dtheta(ind_sm_fl),'.g');
    hold on
    plot(dx(ind_bg_fl | ind_slightly),dtheta(ind_bg_fl | ind_slightly),'.r');
%     plot(dx(ind_slightly),dtheta(ind_slightly),'.m');
    plot(dx(ind_C),dtheta(ind_C),'.c');
    plot(dx(ind_U),dtheta(ind_U),'.y');
    plot(dx(ind_omega),dtheta(ind_omega),'.');
    plot(dx(ind_closed),dtheta(ind_closed),'.','MarkerEdgeColor',[.7 .7 .7]);
    hold off
    xlabel('width (dx)');ylabel('dtheta');
    
    subplot(3,4,6)
    plot(dx(ind_late & ind_glass),dtheta(ind_late & ind_glass),'.',dx(ind_sm & ind_glass),dtheta(ind_sm & ind_glass),'.g',dx(ind_bg_fl & ind_glass),dtheta(ind_bg_fl & ind_glass),'.r')
    xlabel(sprintf('dx with z>%.2f', glass));ylabel('dtheta');
 
    subplot(3,4,7)
    plot(rfit(ind_late),dtheta(ind_late),'.',rfit(ind_sm),dtheta(ind_sm),'.g',rfit(ind_bg_fl),dtheta(ind_bg_fl),'.r')
    xlabel('rfit');ylabel('dtheta');
    
    subplot(3,4,8)
    plot(rfit(ind_late & ind_glass),dtheta(ind_late & ind_glass),'.',rfit(ind_sm & ind_glass),dtheta(ind_sm & ind_glass),'.g',rfit(ind_bg_fl & ind_glass),dtheta(ind_bg_fl & ind_glass),'.r')
    xlabel(sprintf('rfit with z>%.2f', glass));ylabel('dtheta');

    subplot(3,4,9)
    plot(rfit(ind_late),dtheta10(ind_late),'.',rfit(ind_sm),dtheta10(ind_sm),'.g',rfit(ind_bg_fl),dtheta10(ind_bg_fl),'.r')
    xlabel('rfit');ylabel('dtheta10');
    legend('late','small','big_ flat');
    
    subplot(3,4,10)
    plot(dx(ind_late & ind_glass & ind_10),dtheta(ind_late & ind_glass & ind_10),'.',dx(ind_sm & ind_glass & ind_10),dtheta(ind_sm & ind_glass & ind_10),'.g',dx(ind_bg_fl & ind_glass & ind_10),dtheta(ind_bg_fl & ind_glass & ind_10),'.r')
    xlabel(strcat(sprintf('dx with z>%.2f and d10<', glass),sprintf('%.2f',topcov)));ylabel('dtheta');    
    
    subplot(3,4,11)
    plot(ori(ind_sm_fl),dz(ind_sm_fl),'.g');
    hold on
    plot(ori(ind_bg_fl | ind_slightly),dz(ind_bg_fl | ind_slightly),'.r');
%     plot(ori(ind_slightly),dz(ind_slightly),'.m');
    plot(ori(ind_C),dz(ind_C),'.c');
    plot(ori(ind_U),dz(ind_U),'.y');
    plot(ori(ind_omega),dz(ind_omega),'.');
    plot(ori(ind_closed),dz(ind_closed),'.k');
    hold off
    xlabel('ori-dtheta');ylabel('height (dz)');
    
    
    
    subplot(3,4,12)
    plot(rfit(ind_late & ind_glass & ind_10),dtheta(ind_late & ind_glass & ind_10),'.',rfit(ind_sm & ind_glass & ind_10),dtheta(ind_sm & ind_glass & ind_10),'.g',rfit(ind_bg_fl & ind_glass & ind_10),dtheta(ind_bg_fl & ind_glass & ind_10),'.r')
    xlabel(strcat(sprintf('rfit with z>%.2f and d10<', glass),sprintf('%.2f',topcov)));ylabel('dtheta');    
end    
    


% magnified plots of dtheta vs rfit and vs dx, using color code for 7 time
% windows
% PLUS theta vs r-fit as Ori Avinaom defined it
% PLUS curvature vs coverage

if 1    
    figure(1194)
    clf
    
     subplot(1,2,1)
    plot(dx(ind_sm_fl),dtheta(ind_sm_fl),'.g','MarkerSize',20);
    hold on
    plot(dx(ind_bg_fl | ind_slightly),dtheta(ind_bg_fl | ind_slightly),'.r','MarkerSize',20);
%     plot(dx(ind_slightly),dtheta(ind_slightly),'.m');
    plot(dx(ind_C),dtheta(ind_C),'.c','MarkerSize',20);
    plot(dx(ind_U),dtheta(ind_U),'.y','MarkerSize',20);
    plot(dx(ind_omega),dtheta(ind_omega),'.','MarkerSize',20);
    plot(dx(ind_closed),dtheta(ind_closed),'.','MarkerEdgeColor',[.7 .7 .7],'MarkerSize',20);
    hold off
    xlabel('width (dx)');ylabel('dtheta');
    legend('small flat','big flat','curved C','curved U','curved Omega','closed');
    
     subplot(1,2,2)
    plot(rfit(ind_sm_fl),dtheta(ind_sm_fl),'.g','MarkerSize',20);
    hold on
    plot(rfit(ind_bg_fl | ind_slightly),dtheta(ind_bg_fl | ind_slightly),'.r','MarkerSize',20);
%     plot(rfit(ind_slightly),dtheta(ind_slightly),'.m');
    plot(rfit(ind_C),dtheta(ind_C),'.c','MarkerSize',20);
    plot(rfit(ind_U),dtheta(ind_U),'.y','MarkerSize',20);
    plot(rfit(ind_omega),dtheta(ind_omega),'.','MarkerSize',20);
    plot(rfit(ind_closed),dtheta(ind_closed),'.','MarkerEdgeColor',[.7 .7 .7],'MarkerSize',20);
    hold off
    xlabel('rfit');ylabel('dtheta');
    
    figure (1999)
    clf
    plot(curv(ind_sm_fl),dtheta(ind_sm_fl),'.g','MarkerSize',20);
    hold on
    plot(curv(ind_bg_fl | ind_slightly),dtheta(ind_bg_fl | ind_slightly),'.r','MarkerSize',20);
%     plot(rfit(ind_slightly),dtheta(ind_slightly),'.m');
    plot(curv(ind_C),dtheta(ind_C),'.c','MarkerSize',20);
    plot(curv(ind_U),dtheta(ind_U),'.y','MarkerSize',20);
    plot(curv(ind_omega),dtheta(ind_omega),'.','MarkerSize',20);
    plot(curv(ind_closed),dtheta(ind_closed),'.','MarkerEdgeColor',[.7 .7 .7],'MarkerSize',20);
    hold off
    xlabel('curvature');ylabel('coverage');
    
    figure (1995)
    clf
    plot(ori(ind_sm_fl),rfit(ind_sm_fl),'.g','MarkerSize',20);
    hold on
    plot(ori(ind_bg_fl | ind_slightly),rfit(ind_bg_fl | ind_slightly),'.r','MarkerSize',20);
%     plot(rfit(ind_slightly),dtheta(ind_slightly),'.m');
    plot(ori(ind_C),rfit(ind_C),'.c','MarkerSize',20);
    plot(ori(ind_U),rfit(ind_U),'.y','MarkerSize',20);
    plot(ori(ind_omega),rfit(ind_omega),'.','MarkerSize',20);
    plot(ori(ind_closed),rfit(ind_closed),'.','MarkerEdgeColor',[.7 .7 .7],'MarkerSize',20);
    hold off
    xlabel('theta ori');ylabel('radius');
    
    
    
    
end

% comparison between control and released dtheta-ori histogram
if 0
    figure (101)
    clf
    ctr=histogram(ori_control,'Normalization','probability');
    hold on
    rel=histogram(ori_released,'Normalization','probability');
    hold off
    xlabel('Theta');ylabel('frequency');
    legend('control','released 1 min','Location','northwest');
    

end


% surface calculated from fit or calculated from dx and dz

if 1
 figure (1996)
    clf
    plot(ori(ind_sm_fl),sphere_segm(ind_sm_fl),'.g','MarkerSize',20);
    hold on
    plot(ori(ind_bg_fl | ind_slightly),sphere_segm(ind_bg_fl | ind_slightly),'.r','MarkerSize',20);
%     plot(rfit(ind_slightly),dtheta(ind_slightly),'.m');
    plot(ori(ind_C),sphere_segm(ind_C),'.c','MarkerSize',20);
    plot(ori(ind_U),sphere_segm(ind_U),'.y','MarkerSize',20);
    plot(ori(ind_omega),sphere_segm(ind_omega),'.','MarkerSize',20);
    plot(ori(ind_closed),sphere_segm(ind_closed),'.','MarkerEdgeColor',[.7 .7 .7],'MarkerSize',20);
    hold off
    xlabel('theta ori');ylabel('surface from fit');
    
end

  figure(1997)
    clf
    hist(ori(ind_sm_fl | ind_bg_fl | ind_slightly | ind_C | ind_U | ind_omega | ind_closed))
    xlabel('theta ori');
    
    


 if 0
    figure(195)
    clf
    stages={'small (flat)','big flat','slightly curved','curved C','curved U','curved Omega','closed'};
    subplot(1,2,1)
    notBoxPlot(sphere_segm(ind_sm_fl),1);
    notBoxPlot(sphere_segm(ind_bg_fl),2);
    notBoxPlot(sphere_segm(ind_slightly),3);
    notBoxPlot(sphere_segm(ind_C),4);
    notBoxPlot(sphere_segm(ind_U),5);
    notBoxPlot(sphere_segm(ind_omega),6);
    notBoxPlot(sphere_segm(ind_closed),7);
    ylabel('calculated surface /um2');
    title('all sites');
    set(gca,'XLim',[0.1 7.9],'XTick',[1:7],'XTickLabel',stages);
    
       
    subplot(1,2,2)
    notBoxPlot(sphere_segm_dx_dth(ind_sm_fl),1);
    notBoxPlot(sphere_segm_dx_dth(ind_bg_fl),2);
    notBoxPlot(sphere_segm_dx_dth(ind_slightly),3);
    notBoxPlot(sphere_segm_dx_dth(ind_C),4);
    notBoxPlot(sphere_segm_dx_dth(ind_U),5);
    notBoxPlot(sphere_segm_dx_dth(ind_omega),6);
    notBoxPlot(sphere_segm_dx_dth(ind_closed),7);
    ylabel('calculated surface /um2');
    title('calculated from dx and dtheta - all sites');
    set(gca,'XLim',[0.1 7.9],'XTick',[1:7],'XTickLabel',stages);
    
 end
    


   
    
 
end
   
    
