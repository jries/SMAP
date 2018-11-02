global se


saveon=1;

%segment with free roi
% if saveon==1
figure(24)
file=se.files(se.indexFromID(se.files,se.currentfile.ID));
imagesc(file.image.image)
colormap  hot
h = imfreehand;
bw=createMask(h);
bw(1,:)=false;bw(:,1)=false;
% end

%%

cutoff=80;
sigma=2;

pf=file.image.parameters.pixrec;
sf=size(file.image.image);

cells=se.cells;
for cn=1:length(cells)
    cell=cells(cn);
%     srim=sum(double(cell.image.image),3);
    srim=cell.image.layers(1).images.rawimage.image;
    srims=size(srim);


    
    h=fspecial('gaussian',2*sigma,sigma);
    srimf=filter2(h,srim);

    maximaout=NMS2DBlockCcall(srimf,sigma);
    pc=cell.image.parameters.pixrec;
    sc=size(cell.image.image)/2;
    
    
    
    ysite=maximaout(:,1)*pc+cell.image.rangey(1)*1000;
    xsite=maximaout(:,2)*pc+cell.image.rangex(1)*1000;
    
    %initial mask
    ym=ysite-file.image.rangey(1)*1000;
    xm=xsite-file.image.rangex(1)*1000;
    xr=round(xm/pf);
    yr=round(ym/pf);
    xr(xr<1)=1;yr(yr<1)=1;xr(xr>sf(1))=1;yr(yr>sf(1))=1;
    
    linind=sub2ind(sf(1:2),yr,xr);
    indinroi=bw(linind);
    % cutoff=max(maximaout(:,3));
    indgood=maximaout(:,3)>=cutoff;
    indgood=indgood&indinroi;
%     newsites=maximaout(indgood,:);

xsite=xsite(indgood);
ysite=ysite(indgood);

    figure(223)
    subplot(1,2,1)
    imagesc(srim)
    hold on
    plot(maximaout(indgood,2),maximaout(indgood,1),'wo')
    hold off
        subplot(1,2,2)
    imagesc(srimf)
    hold on
    plot(maximaout(indgood,2),maximaout(indgood,1),'wo')
    hold off
    colorbar
    
    if saveon
    %convert to micrometers
%     poscell=cell.pos;
%     newsitesum=(newsites-srims(1)/2)*cell.globpar.pixrec/1000;
%     possites=newsitesum(:,[2 1]);
%     possites(:,1)=possites(:,1)+poscell(1);
%     possites(:,2)=possites(:,2)+poscell(2);

    for k=1:length(xsite)
        thissite=interfaces.SEsites;
        thissite.pos=[xsite(k) ysite(k)];
        thissite.info.cell=cell.ID;
        thissite.info.filenumber=cell.info.filenumber;
        % thissite.cellnumber=sitepar.currentcell.number;
%         thissite.number=sitepar.sitelist.cellnumber+1;
        se.addSite(thissite)

    end
    else
        waitforbuttonpress
    end
end
se.processors.preview.updateSitelist
