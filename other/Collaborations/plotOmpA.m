tt=readtable('/Users/ries/Data/Collaborations/OmpA/OmpA_mito_cc_output.xlsx');

omm=tt.ompaMito./tt.locsmito;
omb=tt.ompaBg./tt.area_um2;
omA=tt.locsOmpa./tt.area_um2;

ombR=tt.ompaBgRand./tt.area_um2;
ommR=tt.ompaMitoRand./tt.locsmito;

c=tt.Exp;

mito=strcmp(tt.Comments,'Mito');
bg=contains(tt.Comments,'BG');
normal=strcmp(tt.Comments,'');
what(bg,1)=0;
what(normal,1)=2;
what(mito,1)=1;

figure(88)
subplot(2,3,2)
hold off
subplot(2,3,1)
hold off
subplot(2,3,3)
hold off
subplot(2,3,4)
hold off
subplot(2,3,5)
hold off
subplot(2,3,6)
hold off

whatlegend={'bg','mito','bg+mito'};
explegend={'WT','infected','expressed'};
whatstr={'o','x','+'};
expstr={'r','k','b'};

whatstr={'o','*','x'};
legendtext={};
for whatind=0:2
    for expind=0:2
        subplot(2,3,2)
        indh=(c==expind & what ==whatind);
        plot(omA(indh),omb(indh),[whatstr{whatind+1} expstr{expind+1}])
        hold on
        subplot(2,3,1)
        plot(omA(indh),omm(indh),[whatstr{whatind+1} expstr{expind+1}])
        hold on
        legendtext{end+1}=[whatlegend{whatind+1} ': ' explegend{expind+1}];
    end
end

legendtext2={};
for expind=0:2
    subplot(2,3,5)
    indh=(c==expind & what ==2);
    plot(omA(indh),omb(indh),[whatstr{whatind+1} expstr{expind+1}])
    hold on
    subplot(2,3,4)
    plot(omA(indh),omm(indh),[whatstr{whatind+1} expstr{expind+1}])
    hold on
    legendtext2{end+1}=[ explegend{expind+1}];
    
    subplot(2,3,6)
     plot(omb(indh),omm(indh),[whatstr{whatind+1} expstr{expind+1}])
    hold on
end
    


% plot(omA(c==0),omb(c==0),'bo')
% hold on
% plot(omA(c==1),omb(c==1),'r+')
% plot(omA(c==2),omb(c==2),'gx')

legend(legendtext)

ylabel('ompa at bg / area')
xlabel('ompa / area')

subplot(2,3,1)

legend(legendtext)

ylabel('ompa at mito / mito locs')
xlabel('ompa / area')

subplot(2,3,4)
legend(legendtext2)
ylabel('ompa at mito / mito locs')
xlabel('ompa / area')

subplot(2,3,5)
legend(legendtext2)
ylabel('ompa at bg / area')
xlabel('ompa / area')

subplot(2,3,6)
legend(legendtext2)
xlabel('ompa at bg / area')
ylabel('ompa at mito / mito locs')

figure(89)
legendtext3={};
for expind=0:2
    indh=(c==expind & what ==2);
     plot(omb(indh),omm(indh),[whatstr{expind+1} 'k'])
     legendtext3{end+1}=[explegend{expind+1}];
    hold on
    plot(ombR(indh),ommR(indh),[whatstr{expind+1} 'r'])
    legendtext3{end+1}=[explegend{expind+1} ' rand'];
end
legend(legendtext3)
xlabel('ompa at bg / area')
ylabel('ompa at mito / mito locs')