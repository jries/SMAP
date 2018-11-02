global se

savetiff=true;
savecoordinates=true;
savelist=true;
savesml=true;

sites=se.sites;
ls=length(sites);
s1=sites(1);
im1=s1.image.image;
sim=size(im1);

% nlocs=zeros(ls,1);
% corners=zeros(ls,1);
filenumber=zeros(ls,1);


[f,p]=uiputfile('NPCsimul.tif');
[~ ,f]=fileparts(f);
% f='test'
r=[];
% rlocs=[];
rall=[];
    indtab=1;
    indlocs=1;
for ff=1:length(se.files)
    indim=1;
    imall=zeros(sim(1),sim(2),ls,'uint8');
    for k=1:ls
        sh=sites(k);
        if sh.info.filenumber==ff
            imall(:,:,indim)=uint8(sum(sh.image.image,3)*255/3);
            r.nlocs(indtab)=sh.evaluation.generalStatistics.Nlayers;
            r.cornersfiltered(indtab)=sh.evaluation.NPCLabelingQuantify.numcornersfiltered;
            r.corners(indtab)=sh.evaluation.NPCLabelingQuantify.numcorners;
            
            r.assignedcorners(indtab)=sh.evaluation.NPCLabelingQuantify.numbercornerassined;
            r.numfoundint(indtab)=sh.evaluation.NPCLabelingQuantify.numfoundint;
            
            r.numfoundrat(indtab)=sh.evaluation.NPCLabelingQuantify.numfoundrat;
            r.plabel(indtab)=se.files(ff).info.simulationParameters.labeling_efficiency;
            r.reactivations(indtab)=se.files(ff).info.simulationParameters.blinks;
            r.photons(indtab)=se.files(ff).info.simulationParameters.photons;
            r.background(indtab)=se.files(ff).info.simulationParameters.background;
            
            r.filenumber(indtab)=ff;
            r.ind(indtab)=indtab;

            
            %full list of localizations
            if isfield(sh.evaluation.NPCLabelingQuantify,'coordinates')
            coordinates=sh.evaluation.NPCLabelingQuantify.coordinates;
            numl=length(coordinates.rho);
            rangeh=(indlocs:indlocs+numl-1)';
            rall.ind(rangeh,1)=indtab;
            rall.rho(rangeh,1)=coordinates.rho;
            rall.theta(rangeh,1)=coordinates.theta;
            rall.drho(rangeh,1)=coordinates.drho;
            rall.dtheta(rangeh,1)=coordinates.dtheta;   
            rall.x(rangeh,1)=coordinates.x;
            rall.y(rangeh,1)=coordinates.y;
            rall.xyerr(rangeh,1)=coordinates.drho;
            rall.cornersfiltered(rangeh,1)=r.cornersfiltered(indtab);
            
            rall.labels(rangeh,1)=r.plabel(indtab)*100;
            rall.blinks(rangeh,1)=r.reactivations(indtab);
            rall.photons(rangeh,1)=r.photons(indtab);
            
             indlocs=indlocs+numl;
            end
            
            indtab=indtab+1;
            indim=indim+1;
           
            
        end

    end
    imall(:,:,indim:end)=[];
    labels=num2str(r.plabel(indtab-1)*100,'%2.0f');
    phots=num2str(r.photons(indtab-1),'%3.0f');
    blinks=num2str(r.reactivations(indtab-1),'%3.0f');
    filename=[num2str(ff) '_' f 'L' labels(1:2) 'P' phots 'B' blinks '.tif']
    if savetiff
        saveastiff(imall,[p filesep filename]);
    end
end
fn=fieldnames(r);
for k=1:length(fn)
    r.(fn{k})=r.(fn{k})';
end
tout=struct2table(r);

% tout=table((1:ls)',r.nlocs,nlocs,corners,filenumber);


[~,ff]=fileparts(f);
if savelist
writetable(tout,[p filesep ff '.csv']);
end

if savecoordinates
    tout=struct2table(rall);
    writetable(tout,[p filesep ff '_c.csv']);
end

if savesml
    l=se.locData;
    l.savelocs([p filesep ff '_sml.mat'])
end
% 
