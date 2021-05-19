%% load
if ~exist('pfad','var')
    pfad='';
end
[file,pfad]=uigetfile([pfad '*.tif']);
img=imread([pfad file]);

%% segment
figure(88);imagesc(img)
axis('equal')
roinumber=1;

%% add a line
roi(roinumber)=drawpolyline;
roinumber=roinumber+1;

%% assemble output table
x=[];y=[];fibril=[]; interp=[];
tab=table(x,y,fibril,interp);
for k=1:length(roi)
    roih=roi(k);
    pos=roih.Position;
    

    
    newind=length(tab.x)+1:length(tab.x)+size(pos,1);
    tab.x(newind)=pos(:,1);
    tab.y(newind)=pos(:,2);
    tab.fibril(newind)=k;
    tab.interp(newind)=0;
    
    %interpolate position
    mask=roih.createMask;
    linind=find(mask);
    [yi,xi]=ind2sub(size(mask),linind);
    newind=length(tab.x)+1:length(tab.x)+length(xi);
    tab.x(newind)=xi;
    tab.y(newind)=yi;
    tab.fibril(newind)=k;
    tab.interp(newind)=1;   
end

%% save outputtable
outname=strrep(file,'.tif','.csv');
% outname='newname.csv';
writetable(tab,[pfad outname]);