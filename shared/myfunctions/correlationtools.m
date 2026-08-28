classdef correlationtools<handle
    % ct=correlationtools(x,bin), then all is calculated, stored at
    % properties
    % ct.fft.f, ct.fft.mag ct.fft.peaks, ct.fft.mainpeak
    % ct.plot('fft',arguments)
    % cross-correlation: function taht takes two correlation objects (or input: additional one),
    % outputs cc structure
    %AC change to two-sided, then same for CC
    properties
        x
        bin
        fft
        profile
        autocorrelation
        crosscorrelation
        periodguess
        sigmafilter
        maxcorr
        
        fftav
        autocorrelationav
        crosscorrelationav
        numav
        
    end
    methods
        function obj = correlationtools(posin,bin,options)
            arguments
                posin
                bin
                options.periodguess=[];
                options.sigmafilter=2;
                options.maxcorr=[];
            end
            obj.periodguess=options.periodguess;
            obj.sigmafilter=options.sigmafilter;
            obj.maxcorr=options.maxcorr;
            if isa(posin, 'struct')
                if isfield(posin,'xnmline')
                    obj.x=posin.xnmline;
                else
                    obj.x=posin.xnm;
                end
            else
                obj.x=posin;
            end
            obj.bin=bin;
            obj.calculateprofile;
            obj.calculatefft;
            obj.calculateautocorrelation;
            
        end
        function calculateprofile(obj)
            obj.profile.counts = histcounts(obj.x, obj.nx);  

            proffilt=smoothdata(obj.profile.counts,'gaussian',6*obj.sigmafilter);
            [ypos,xPeaksIdx]=findpeaks(proffilt);
            [yRPeaks,xRPeaks] = refinepeaks(proffilt,xPeaksIdx,obj.nxplot);

            obj.profile.nx=obj.nxplot;
            obj.profile.n=obj.nx;
            obj.profile.filter.counts=proffilt;
            obj.profile.filter.nx=obj.nxplot;
            obj.profile.filter.peaks.x=xRPeaks;
            obj.profile.filter.peaks.y=yRPeaks;
        end
        function calculateautocorrelation(obj)
            ac=xcorr(obj.profile.counts,obj.profile.counts)/mean(obj.profile.counts)^2;
            % ach=ac(ceil(end/2):end)/max(ac);
            % nac=(0:length(ach)-1)'*obj.bin;
            lm=floor(length(ac)/2);
            nac=(-lm:lm)'*obj.bin;
            obj.autocorrelation.mag=ac;
            obj.autocorrelation.lag=nac;
            obj.autocorrelation.nx=nac;
            analyzeautocorrelation(obj,'autocorrelation',1)
        end
        function calculatecrosscorrelation(obj,ca2)
            cc=xcorr(obj.profile.counts,ca2.profile.counts)/mean(obj.profile.counts)/mean(ca2.profile.counts);
            lm=floor(length(cc)/2);
            nac=(-lm:lm)'*obj.bin;
            % ach=cc(ceil(end/2):end)/max(cc);
            % nac=(0:length(ach)-1)'*obj.bin;
            obj.crosscorrelation.mag=cc;
            obj.crosscorrelation.lag=nac;
            obj.crosscorrelation.nx=nac;
            analyzeautocorrelation(obj,'crosscorrelation',1)
        end

        function analyzeautocorrelation(obj,f,norm)
            if isempty(obj.(f))
                return
            end
            ach=obj.(f).mag/norm;
            nac=obj.(f).lag;
            [yRPeaks,xPeaksIdx]=findpeaks(detrend(ach,2));
            xPeaksIdx(xPeaksIdx<3)=[]; 
            [obj.(f).period,obj.(f).periodmag]=obj.selectpeak(nac(xPeaksIdx),ach(xPeaksIdx));


            obj.(f).peaks.x=nac(xPeaksIdx);
            obj.(f).peaks.y=ach(xPeaksIdx);

            acs=ach;
            for k=1:3
                acs=detrend(acs,2);
                acs=xcorr(acs-mean(acs),acs-mean(acs));
                acs=acs/max(acs);
                acs=acs(ceil(end/4):ceil(end/4*3));
            end

            [yRPeaks,xPeaksIdx]=findpeaks(acs);
            [yRPeaks,xRPeaks] = refinepeaks(acs,xPeaksIdx,nac);
            [obj.(f).filter.period,obj.(f).filter.periodmag]=obj.selectpeak(xRPeaks,yRPeaks);

            % obj.autocorrelation.filter.period=periodfilt;
            obj.(f).filter.peaks.x=xRPeaks;
            obj.(f).filter.peaks.y=yRPeaks;
            obj.(f).filter.mag=acs;
            obj.(f).filter.nx=nac;
            

            % acmulticorr=0.5*(1+acs(ceil(end/2):ceil(end/2)+length(nac)-1));
            % acmulticorr=acmulticorr-min(acmulticorr);
            % acmulticorr=acmulticorr/mean(acmulticorr)*mean(ach);
        end
        function calculatefft(obj)
            [f, mag] = fft_one_sided(obj.nx,obj.profile.counts);
            obj.fft.mag=mag;
            obj.fft.f=f;
            obj.fft.xplot=f;
            obj.analyzefft('fft',1)
        end
        function analyzefft(obj,fi,norm)
            mag=obj.(fi).mag/norm;
            f=obj.(fi).f;
            [yRPeaks,xPeaksIdx]=findpeaks(mag);
            [yRPeaks,xRPeaks] = refinepeaks(mag,xPeaksIdx,f);
            periods=1./xRPeaks;
            [obj.(fi).period,obj.(fi).periodmag]=obj.selectpeak(periods,yRPeaks);
            
            obj.(fi).peaks.x=xRPeaks;
            obj.(fi).peaks.y=yRPeaks;
        end
        % function crosscorrelationi(obj)
        % end
        function nx=nx(obj)
            nx=double(min(obj.x):obj.bin:max(obj.x))';
        end
        function nxp=nxplot(obj)
            nx=obj.nx;
            nxp=nx(1:end-1)+obj.bin/2;
        end
        function plot(obj,prop,options)
            arguments
                obj
                prop
                % options.filter=0;
                options.plotpeaks=true;
                options.axis=gca;
                options.color='b';
                options.average=false;
            end
            ff="%2.1f";
            switch prop
                case 'profile'
                    plot(options.axis,obj.profile.nx,obj.profile.counts,'Color',options.color,'DisplayName',"profile")
                    hold(options.axis,'on')
                    plot(options.axis,obj.profile.filter.nx,obj.profile.filter.counts,'-','Color',options.color,'DisplayName','filt','LineWidth',2)
                    if options.plotpeaks
                        xpos=obj.profile.filter.peaks.x;
                        dx=diff(xpos);
                        tpos=(xpos(1:end-1)+xpos(2:end))/2;
                        text(options.axis,tpos,max(obj.profile.filter.counts)+0*tpos,compose(ff,dx))
                        plot(options.axis,obj.profile.filter.peaks.x,obj.profile.filter.peaks.y,'o','Color',options.color)
                    end
                     xlabel(options.axis,"position (nm)")
                     ylabel(options.axis,"counts")
                case {'fft','fftav'}
                    ft=obj.(prop);
                    if options.average
                        norm=obj.numav;
                    else
                        norm=1;
                    end
                    plot(options.axis, ft.f, ft.mag/norm,'Color',options.color,'DisplayName',"P: "+num2str(ft.period,ff));
                    hold(options.axis,'on')
                    if options.plotpeaks
                        plot(options.axis,1/ft.period,ft.periodmag,'o','Color',options.color,'HandleVisibility','off')
                        plot(options.axis,ft.peaks.x,ft.peaks.y,'x','Color',options.color,'HandleVisibility','off')

                    end
                    legend(options.axis)
                    xlim(options.axis,[0 10/ft.period])
                    xlabel(options.axis,"frequency (1/nm)")
                    ylabel(options.axis,"magnitude")

                case {'autocorrelation','autocorrelationav','crosscorrelation','crosscorrelationav'}
                    ac=obj.(prop);
                    if options.average
                        norm=obj.numav;
                    else
                        norm=1;
                    end
                    h1=plot(options.axis,ac.nx,ac.mag/norm,'Color',options.color,'DisplayName',"correlation, P="+num2str(ac.period,ff));
                    hold(options.axis,'on')


                    acf=ac.filter.mag;
                    off=min(acf);
                    acf=acf-off;
                    fac=mean(acf)*mean(ac.mag/norm);
                    acf=acf*fac;
                    h2=plot(options.axis,ac.nx,acf,':','Color',options.color,'DisplayName',"multi-c: P="+num2str(ac.filter.period,ff));
                    legend(options.axis)
                    if ~isempty(obj.maxcorr)
                        xlim(options.axis,[-obj.maxcorr obj.maxcorr])
                    end

                    if options.plotpeaks
                        plot(options.axis,ac.peaks.x,(ac.peaks.y),'x','Color',options.color,'HandleVisibility','off')
                        plot(options.axis,ac.period,ac.periodmag,'o','Color',options.color,'HandleVisibility','off')
                        plot(options.axis,ac.filter.peaks.x,(ac.filter.peaks.y-off)*fac,'x','Color',options.color,'HandleVisibility','off')
                        plot(options.axis,ac.filter.period,(ac.filter.periodmag-off)*fac,'o','Color',options.color,'HandleVisibility','off')
                    end

                    xlabel(options.axis,"dx (nm)")
                    ylabel(options.axis,"g(r)")
                otherwise
                    disp('not implemented')
            end
           

        end

        function add(obj,co)
            if isempty(obj.autocorrelationav) %also true for fft, always together
                obj.autocorrelationav=obj.autocorrelation;
                obj.crosscorrelationav=obj.crosscorrelation;
                obj.fftav=obj.fft;
                obj.numav=1;

            else
                obj.autocorrelationav=addminlenstruc(obj.autocorrelationav,co.autocorrelation,{'mag'});
                obj.autocorrelationav=copyminlenstruc(obj.autocorrelationav,co.autocorrelation,{'lag','nx'});
                % if ~isempty(co.crosscorrelation)&& ~isempty(obj.crosscorrelationav)
                obj.crosscorrelationav=addminlenstruc(obj.crosscorrelationav,co.crosscorrelation,{'mag'});
                obj.crosscorrelationav=copyminlenstruc(obj.crosscorrelationav,co.crosscorrelation,{'lag','nx'});
                % end
                obj.fftav=addminlenstruc(obj.fftav,co.fft,{'mag'});
                obj.fftav=copyminlenstruc(obj.fftav,co.fft,{'f','xplot'});
                obj.numav=obj.numav+1;
                %add
            end
            analyzeautocorrelation(obj,'autocorrelationav',obj.numav)
            % if ~isempty(obj.crosscorrelationav)
            analyzeautocorrelation(obj,'crosscorrelationav',obj.numav)
            % end
            analyzefft(obj,'fftav',obj.numav)
            %recalculate peaks, averages etc

        end
        function [peak,peaky]=selectpeak(obj,xpeaks,ypeaks)
            if ~isempty(obj.periodguess) && obj.periodguess~=0
                [~,imxp]=min(abs(xpeaks-obj.periodguess));
            else
                ypeaks(abs(xpeaks)==0)=0; %exclude zero peak
                [~,imxp]=max(ypeaks);
            end
            peak=xpeaks(imxp);
            peaky=ypeaks(imxp);
        end
           
    end
end

function so=addminlenstruc(s1,s2,fields)
if isempty(s1) || isempty(s2)
    so=[];
    return
end
so=s1;
for k=1:length(fields)
    so.(fields{k})=addminlen(s1.(fields{k}),s2.(fields{k}));
end

end
function so=addminlen(s1,s2)
if length(s1)<=length(s2)
    so=s1+s2(1:length(s1));
else
    so=s1(1:length(s2))+s2;
end
end


function so=copyminlenstruc(s1,s2,fields)
if isempty(s1) || isempty(s2)
    so=[];
    return
end
so=s1;
for k=1:length(fields)
    if length(s1.(fields{k}))<length(s2.(fields{k}))
        so.(fields{k})=s1.(fields{k});
    else
        so.(fields{k})=s2.(fields{k});
    end
end
end



function [f, mag1] = fft_one_sided(t, x)
    % Returns one-sided frequency axis and magnitude |X(f)| for real x(t)
    t = t(:); x = x(:);
    dt = diff(t);
    if max(abs(dt - mean(dt))) > 1e-6*mean(dt)
        error('Time vector is not uniformly sampled.');
    end
    Fs = 1/mean(dt);
    N = numel(x);
    X = fft(x);
    mag1 = abs(X(1:floor(N/2)+1))/N;
    if N > 2
        mag1(2:end-1) = 2*mag1(2:end-1);
    end
    f = (0:floor(N/2))*(Fs/N);
end