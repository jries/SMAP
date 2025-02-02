classdef Cluster_MINFLUX_Roi<interfaces.SEEvaluationProcessor
    % LINEPROFILE Calculates profiles along a linear ROI and fits it with a
    % model of choice. Flat: step function convolved with Gaussian
    % (=Erf). Disk: Projection of a homogeneously filled disk, convolved
    % with Gaussian. Ring: Projection of a ring, convolved with
    % Gaussian. Distance: Two Gaussians in a distance d.
    methods
        function obj=Cluster_MINFLUX_Roi(varargin)        
            obj@interfaces.SEEvaluationProcessor(varargin{:});
            % obj.inputParameters={'sr_layerson','linewidth_roi','znm_min','znm_max','sr_pixrec','layernames'};
            % obj.showresults=true;
        end
        
        function out=run(obj,p)      
            usefields={'xnm','ynm','znm','time','groupindex','numberInGroup','filenumber','efo','cfr','eco','ecc','efc','tid','fbg','phot'};
            [locsfind,indr]=obj.getLocs({'tid','groupindex','filenumber'},'layer',find(obj.getPar('sr_layerson')),'size',obj.getPar('se_siteroi')/2,'removeFilter',{'time'});
            locs=obj.locData.getloc(usefields,'layer',find(obj.getPar('sr_layerson')),'Position','all','removeFilter',{'filenumber','time'});
            switch p.link.selection
                case 'StepsMINFLUX'
                    ind=obj.site.evaluation.StepsMINFLUX.steps.allindices;
                    locs=obj.locData.loc;
                case 'id'
                    ind=locs.tid==mode(locsfind.tid) & locs.filenumber==mode(locsfind.filenumber);
                case 'group'
                    ind=locs.groupindex==mode(locsfind.groupindex);
                case 'all'
                    ind=locs.filenumber==mode(locsfind.filenumber) & (locs.xnm-obj.site.pos(1)).^2+(locs.ynm-obj.site.pos(2)).^2<(obj.getPar('se_siteroi')/2)^2;
                    % ind=true(size(locsfind.tid));
                    % locs=obj.locData.loc;
            end
            if p.skipfirst>0
                ind=find(ind);
                ind=ind(p.skipfirst+1:end);
            end
            
            filelist=obj.getPar('filelist_short');
            filename=filelist.String{mode(locs.filenumber(ind))};
            dt=diff(locs.time(ind));
            timedt=locs.time(ind);
            dtmin=min(dt(dt>0));
            dtmedian=median(dt);
            dtmean=mean(dt);
            efo=median(locs.efo(ind),'omitnan');
            cfr=median(locs.cfr(ind),'omitnan');
            eco=median(locs.eco(ind),'omitnan');
            ecc=median(locs.ecc(ind),'omitnan');
            efc=median(locs.efc(ind),'omitnan');
            fbg=median(locs.fbg(ind),'omitnan');
            nlocs=length(locs.time(ind));
            ontime=max(locs.time(ind))-min(locs.time(ind));

            

            sigmax=std(locs.xnm(ind));sigmay=std(locs.ynm(ind));
            sxdetrend=std(diff(locs.xnm(ind)))/sqrt(2);sydetrend=std(diff(locs.ynm(ind)))/sqrt(2);
            [~, sxrobust]=robustMean(locs.xnm(ind)); [~, syrobust]=robustMean(locs.ynm(ind));
            
            %graphs
                 ff='%2.1f';
                 mx=mean(locs.xnm(ind));
                 my=mean(locs.ynm(ind));

            dxplot=locs.xnm(ind)-mx;
            dyplot=locs.ynm(ind)-my;
            ltime=locs.time(ind)-min(locs.time(ind));
            
            
            if p.groupon && strcmp(p.link.selection,'all')
                ids=(locs.tid(ind));
                allids=unique(ids);
                dxph=zeros(length(allids),1); dyph=dxph; lth=dxph;sxh=dxph;syh=dxph;
                for k=1:length(allids)
                    indh=(ids==allids(k));
                    dxph(k)=mean(dxplot(indh));
                    dyph(k)=mean(dyplot(indh));
                    sxh(k)=std(dxplot(indh))/sqrt(sum(indh));
                    syh(k)=std(dyplot(indh))/sqrt(sum(indh));
                    lth(k)=mean(ltime(indh));
                end
                dxplot=dxph;
                dyplot=dyph;
                ltime=lth;
            else
                sxh=[]; syh=[];
            end

           
            tzero=ltime*0;

            axy=obj.setoutput('xy');
            plot(axy,dxplot,dyplot,'c')
            hold(axy,'on')
            plot(axy,dxplot,dyplot,'k.')
            hold(axy,'off')
            xlabel(axy,'x (nm)')
            ylabel(axy,'y (nm)')
            axis(axy,'equal')

            axx=obj.setoutput('x');
            
            plot(axx,ltime,dxplot,ltime,tzero,'k', ...
                ltime,tzero+sigmax,'c',ltime,tzero-sigmax,'c',...
                ltime,tzero+sxdetrend,'m',ltime,tzero-sxdetrend,'m')
            if ~isempty(sxh)
                hold(axx,'on')
                errorbar(axx,ltime, dxplot,sxh,'k.')
                 hold(axx,'off')
            end
            legend(axx,'data','','std','','robust std','','detrend std','')
            xlabel(axx,'time (ms)')
            ylabel(axx,'x (nm)')
               title(axx,['std(x) = ' num2str(sigmax,ff) ' nm, std(x) robust = ' num2str(sxrobust,ff) ' nm, std(x) detrend = ' num2str(sxdetrend,ff) ' nm.'])

             plotavtrace(p,axx, ltime,locs.xnm(ind)-mean(locs.xnm(ind)));


            axy=obj.setoutput('y');
            
            plot(axy,ltime,dyplot,ltime,tzero,'k', ...
                ltime,tzero+sigmay,'c',ltime,tzero-sigmay,'c',...
                ltime,tzero+sydetrend,'m',ltime,tzero-sydetrend,'m')

            if ~isempty(syh)
                hold(axy,'on')
                errorbar(axy,ltime, dyplot,syh,'k.')
                 hold(axy,'off')
            end

            legend(axy,'data','','std','','robust std','','detrend std','')
            xlabel(axy,'time (ms)')
            ylabel(axy,'y (nm)')
            title(axy,['std(y) = ' num2str(sigmay,ff) ' nm, std(y) robust = ' num2str(syrobust,ff) ' nm, std(y) detrend = ' num2str(sydetrend,ff) ' nm.'])
            plotavtrace(p,axy, ltime,locs.ynm(ind)-mean(locs.ynm(ind)));
            
            axbb=obj.setoutput('xbin');
            plotstdbin(p,locs.xnm(ind)-mx, locs.phot(ind), axbb)
            axbby=obj.setoutput('ybin');
            plotstdbin(p,locs.ynm(ind)-my, locs.phot(ind), axbby)

            if isfield(locs,'znm') && ~isempty(locs.znm) & any(locs.znm ~= 0)
                sigmaz=std(locs.znm(ind));
                szdetrend=std(diff(locs.znm(ind)))/sqrt(2);
                [~, szrobust]=robustMean(locs.znm(ind)); 
                axz=obj.setoutput('z');
                plot(axz,ltime,locs.znm(ind)-mean(locs.znm(ind)),ltime,0*locs.time(ind),'k', ...
                    ltime,0*locs.time(ind)+sigmaz,'c',ltime,0*locs.time(ind)-sigmaz,'r',...
                    ltime,0*locs.time(ind)+szrobust,'r',ltime,0*locs.time(ind)-szrobust,'c',...
                    ltime,0*locs.time(ind)+szdetrend,'m',ltime,0*locs.time(ind)-szdetrend,'m')
                plotavtrace(p,axz, ltime,locs.znm(ind)-mean(locs.znm(ind)));
                legend(axz,'data','','std','','robust std','','detrend std','')
                xlabel(axz,'time (ms)')
                ylabel(axz,'z (nm)')
                title(axz,['std(z) = ' num2str(sigmaz,ff) ' nm, std(z) robust = ' num2str(szrobust,ff) ' nm, std(z) detrend = ' num2str(szdetrend,ff) ' nm.'])
                tz='nlocs \t on-time \t dtmin \t dtmedian \t <dt> \t sigmax \t sigmay \t sigmaz \t sigmax robust \t sigmay robust \t sigmaz robust \t sigmax detrend \t sigmay detrend \t sigmaz detrend \t efo med \t cfr med  \t eco med  \t ecc med  \t efc med \t fbg med  \t filename';
                tsig=['\t' num2str(sigmax) '\t' num2str(sigmay) '\t' num2str(sigmaz)  '\t' num2str(sxrobust)  '\t' num2str(syrobust) '\t' num2str(szrobust)  '\t' num2str(sxdetrend)  '\t' num2str(sydetrend) '\t' num2str(szdetrend)];
                outsig.sigmax=sigmax;outsig.sigmay=sigmay;outsig.sigmaz=sigmaz;outsig.sxrobust=sxrobust;outsig.syrobust=syrobust;outsig.szrobust=szrobust;outsig.sxdetrend=sxdetrend;outsig.sydetrend=sydetrend;outsig.szdetrend=szdetrend;

                axbbz=obj.setoutput('zbin');
                plotstdbin(p,locs.znm(ind)-mz, locs.phot(ind), axbbz)
            else
                tz='nlocs \t on-time \t dtmin \t dtmedian \t <dt> \t sigmax \t sigmay \t sigmax robust \t sigmay robust \t sigmax detrend \t sigmay detrend \t efo med \t cfr med  \t eco med  \t ecc med  \t efc med \t fbg med \t filename' ;
                tsig=['\t' num2str(sigmax) '\t' num2str(sigmay)  '\t' num2str(sxrobust)  '\t' num2str(syrobust)  '\t' num2str(sxdetrend)  '\t' num2str(sydetrend)];
                outsig.sigmax=sigmax;outsig.sigmay=sigmay;outsig.sxrobust=sxrobust;outsig.syrobust=syrobust;outsig.sxdetrend=sxdetrend;outsig.sydetrend=sydetrend;
            end

            axt=obj.setoutput('time');
            hold(axt,'off')
            histogram(axt,dt,dtmin/2:dtmin:max(quantile(dt,0.995),dtmin/2+2*dtmin))
            hold(axt,'on')
            histogram(axt,dt,0:dtmin*0.1:quantile(dt,0.995))
            xlabel(axt,'dt (ms)')
            ylabel(axt,'frequency')
            title(axt,['dtmin = ' num2str(dtmin,ff) ' ms, dtmedian = ' num2str(dtmedian,ff) ' ms, dtmean = ' num2str(dtmean,ff) ' ms.'])
            
            axdt=obj.setoutput('dt');
            plot(axdt,timedt(2:end),dt)
            xlabel(axdt,'time (ms)')
            ylabel(axdt,'dt (ms)')
            plotavtrace(p,axdt, ltime(2:end),dt);
            
            if ~isempty(locs.efo)
                axe=obj.setoutput('efo');
                plot(axe,timedt,locs.efo(ind))
                xlabel(axe,'time (ms)')
                ylabel(axe,'efo')
                plotavtrace(p,axe, timedt,locs.efo(ind));
            end
            if ~isempty(locs.cfr)
                axc=obj.setoutput('cfr');
                plot(axc,timedt,locs.cfr(ind))
                xlabel(axc,'time (ms)')
                ylabel(axc,'cfr')
                plotavtrace(p,axc, timedt,locs.cfr(ind));
            end
            if ~isempty(locs.eco)
                axec=obj.setoutput('eco');
                plot(axec,timedt,locs.eco(ind))
                xlabel(axec,'time (ms)')
                ylabel(axec,'eco')
                plotavtrace(p,axec, timedt,locs.eco(ind));
            end
            if ~isempty(locs.ecc)
                axcc=obj.setoutput('ecc');
                plot(axcc,timedt,locs.ecc(ind))
                xlabel(axcc,'time (ms)')
                ylabel(axcc,'ecc')
                plotavtrace(p,axcc, timedt,locs.ecc(ind));
            end
            
            out.nocs=nlocs;out.ontime=ontime;out.dtmin=dtmin; out.dtmedian=dtmedian;out.dtmean=dtmean;
            out.sigmax=sigmax;out.sigmay=sigmay;out.sxrobust=sxrobust;out.syrobust=syrobust;out.sxdetrend=sxdetrend;out.sydetrend=sydetrend;
            out.efo=efo;out.cfr=cfr;out.eco=eco;out.ecc=ecc;out.efc=efc;out.fbg=fbg;
            % [~,filenamep]=fileparts(filename);
            out.filenames=filename;
            header=sprintf(tz);
            disp(header)
            results=sprintf([num2str(nlocs)  '\t' num2str(ontime)  '\t' num2str(dtmin)  '\t' num2str(dtmedian) '\t' num2str(dtmean) tsig '\t' num2str(efo) '\t' ...
                 num2str(cfr) '\t' num2str(eco) '\t' num2str(ecc) '\t' num2str(efc) '\t' num2str(fbg) '\t' filename]  );
            clipboard("copy",results)
            % out.clipboard=results;

        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function plotavtrace(p,axx, time, x)
    if p.avtraceon
        oldhold=ishold(axx);
        hold(axx, 'on')
        fun=@(x) median(x(:));
        % y=nlfilter(x,[p.avwin 1], fun);
        % switch p.avmode.selection
        %     case 'median'
                y=runningWindowAnalysisPoints(x,p.avwin,p.avmode.selection);
        %     case 'mean'
        %         y=runningaverage(x,p.avwin);
        % end
        plot(axx, time, y,'k','LineWidth',1.)
        if ~oldhold
            hold(axx,"off")
        end
    end
end

function plotstdbin(p,x, phot, axx)
xh=x;photh=phot;
nump=20;
xb=zeros(nump,1);photb=xb;
for k=1:nump
    xb(k)=std(xh);
    photb(k)=mean(photh);
    xh=sumbintrace2(xh)/2;
    if length(xh)<10
        break
    end
    photh=sumbintrace2(photh);
end
hold(axx,'off')
semilogx(axx,photb(1:k), xb(1:k),'ko-')
hold(axx,'on')
sigsmlm=120./sqrt(photb(1:k));

sigminflux=p.L./sqrt(8*photb(1:k));
semilogx(axx,photb(1:k), sigsmlm,'m--')
semilogx(axx,photb(1:k), sigminflux,'b-.')
ylabel(axx,'std pos (nm)')
xlabel(axx,'photons')
legend(axx,'std','SMLM','MINFLUX')
end

function xb=sumbintrace2(x)
len=floor(length(x)/2)*2;
xb=x(1:2:len-1)+x(2:2:len);
end

function pard=guidef(obj)

pard.linkt.object=struct('String','Selection','Style','text');
pard.linkt.position=[1,1];
pard.link.object=struct('String',{{'StepsMINFLUX','group','id','all'}},'Style','popupmenu','Value',3);
pard.link.position=[1,2];
pard.link.Width=1.5;


pard.skipt.object=struct('String','Skip first','Style','text');
pard.skipt.position=[2,1];
pard.skipfirst.object=struct('String',0,'Style','edit');
pard.skipfirst.position=[2,2];
pard.skipfirst.Width=1;

pard.avtraceon.object=struct('String','average trace','Style','checkbox');
pard.avtraceon.position=[3,1];
pard.avtraceon.Width=2;
pard.avmode.object=struct('String',{{'median','mean'}},'Style','popupmenu');
pard.avmode.position=[3,3];
pard.avmode.Width=2;
pard.avwint.object=struct('String','window size (locs)','Style','text');
pard.avwint.position=[4,2];
pard.avwint.Width=2;
pard.avwin.object=struct('String','20','Style','edit');
pard.avwin.position=[4,4];
pard.avwin.Width=1;

pard.Lt.object=struct('String','scan size L (nm)','Style','text');
pard.Lt.position=[5,1];
pard.Lt.Width=2;
pard.L.object=struct('String','75','Style','edit');
pard.L.position=[5,3];
pard.L.Width=1;


pard.groupon.object=struct('String','group individual tracks (selection all)','Style','checkbox','Value',1);
pard.groupon.position=[6,1];
pard.groupon.Width=4;


pard.plugininfo.description=sprintf('');
pard.plugininfo.type='ROI_Evaluate';
end