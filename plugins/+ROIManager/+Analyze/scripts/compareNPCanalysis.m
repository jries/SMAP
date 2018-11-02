%analyze npc corners
global path
[f,path]=uigetfile([path '*.csv']);
t=readtable([path f],'TreatAsEmpty',{'.','NA','N/A'} );

fields={'truth','assignedcorners','numfoundint','numfoundrat','prediction_icp_svm_','prediction_trr_'};

p.corners=8;p.rings=4;p.fitrange=[3 8];
% labeling efficiencies
plabel=unique(t.plabel);
figure(99)
for field=1:length(fields)
    subplot(3,3,field)
    hold off
    if field==1
        pc=[0 0 0];
    elseif field==2
        pc=pf;
    end
    for k=1:length(plabel)
        indhere=t.plabel==plabel(k);
        cornershere=t.(fields{field})(indhere);
        pf(k)=fitNPClabeling(cornershere,p);
    end

    title([fields{field} ': ' num2str(pf-pc,'%1.3f, ')]);
end

figure(100)
for field=1:length(fields)
    subplot(3,3,field)
    hold off
%     for k=1:length(plabel)
%         indhere=t.plabel==plabel(k);
%         cornershere=t.(fields{field})(indhere);
%         pf(k)=fitNPClabeling(cornershere,p);
%     end
    dc=t.(fields{field})-t.truth;
    histogram(dc,-8:8)
    title([fields{field} ', mean: ' num2str(nanmean(dc), '%1.2f')  ', std: ' num2str(nanstd(dc), '%1.2f')]);
end
