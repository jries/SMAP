classdef dyeblinking<interfaces.DialogProcessor
    % PLUGIN_TEMPLATE Summary of this plugin goes here
    % put a description of your plugin here.
        %replace Plugin_Template by filename   
    properties
        %define class properties if needed
    end
    methods
        function obj=dyeblinking(varargin)   %replace by filename        
            obj@interfaces.DialogProcessor(varargin{:}) 
            obj.showresults=true; %set true, if results are shown by default
            obj.history=false; %if set true, every time the plugin is called, its parameters are added to the history. Makes sense only if plugin changes the data
            % obj.guiselector.show=true; %if true, the selector for simple vs complex GUI is shown.
            % obj.excludeFromSave={'field'}; %dont save or load this GUI component
            % obj.propertiesToSave={'property'}; %save these properties along with the parameters, e.g. in the history
        end      
        function out=run(obj,p)
            out=[];
            % FP blinking, later make plugin out of it
            % maxdTs= 5; %s
            ax1=obj.initaxis('ontime');
            ax2=obj.initaxis('offtime');
                   
            ax3=obj.initaxis('duration');
            ax4=obj.initaxis('onfrac');
            ax5=obj.initaxis('blinks');               
            ax6=obj.initaxis('photdT');
            ax7=obj.initaxis('ondT');

            switch p.showwhat.Value
                case 1
                    layers=find(obj.getPar('sr_layerson'));
                    for k=1:length(layers)
                        files(k)=double(obj.locData.getloc('filenumber','layer',layers(k)).filenumber(1));
                    end
                    leg=obj.getPar('layernames');
                    uselayers=true;
                case 2   
                    files=double(unique(obj.locData.loc.filenumber));
                    uselayers=false;
                    for f=1:length(files)
                        [~, fn]=fileparts(obj.locData.files.file(f).name)
                        leg{f}=strrep(fn,'_',' ');
                    end
            end

            for f=1:length(files)
                obj.locData.files.file(f).name
                exposuret=obj.locData.files.file(f).info.exposure/1000; %s
                maxdTf=p.maxdt/exposuret;
                locd=obj.locData.copy;
                thisf=locd.loc.filenumber==files(f);
                locd.removelocs(~thisf);
                locd.regroup([30 75],maxdTf,2);
                locd.sort('groupindex');
                %%
                gi=locd.loc.groupindex;
                % maxg=gi(end);
                ind1=1;
                gih=gi(1);
                newtrace=false(length(gi),1); newblink=false(length(gi),1);
                ontime=zeros(length(gi),1);offtime=zeros(length(gi),1);
                tracetime=zeros(length(gi),1);numberOfBlinks=zeros(length(gi),1);onfraction=zeros(length(gi),1);
                % dframe=diff(locd.loc.frame);
                for k=1:length(gi) 
                    if gi(k)>gih || k==length(gi) % new track. we lose last localization...
                        indh=ind1:k-1;
                        framesh=locd.loc.frame(indh);
                        newtrace(ind1)=true;
                        tracetime(ind1)=framesh(end)-framesh(1)+1;
                        onfraction(ind1)=length(framesh)/tracetime(ind1);
                        % numberOfBlinks(ind1)=sum(diff(framesh)>1)+1;
                        dfh=[diff(framesh); -1];
                        dfh1=[true;dfh>1];
                        numberOfBlinks(ind1)=sum(dfh1);
                        fdfh=[find(dfh1); length(framesh)+1];
                        for l=1:length(fdfh)-1
                            indbl=indh(fdfh(l));
                            newblink(indbl)=true;
                
                            ontime(indbl)=fdfh(l+1)-fdfh(l);
                            offtime(indbl)=dfh(fdfh(l+1)-1);
                            
                
                        end
                
                        findh=find(indh);
                        gih=gi(k);
                        ind1=k;
                    end
                end
                %%
                % figure(88)
                
                histogram(ax1,ontime(newblink),1:quantile(ontime,0.99),"DisplayStyle","stairs",Normalization="probability")
                xlabel(ax1,'on-times (frames)')
                ax1.YScale=p.yscale.selection;ax1.XLim(1)=0;hold(ax1,'on');
                
               
                offt=offtime(newblink);
                histogram(ax2,offt(offt>0)*exposuret,(0:exposuret:max(offt)*exposuret)+exposuret/20,"DisplayStyle","stairs",Normalization="probability")
                xlabel(ax2,'off-times (s)')
                ax2.YScale=p.yscale.selection;ax2.XLim(1)=0;hold(ax2,'on');
                
                histogram(ax3,tracetime(newtrace)*exposuret,0:exposuret:10,"DisplayStyle","stairs",Normalization="probability")
                xlabel(ax3,'duration of trace (s)')
                ax3.YScale=p.yscale.selection;ax3.XLim(1)=0;ax3.XLim(2)=quantile(tracetime(newtrace)*exposuret,.9);hold(ax3,'on');
          
                histogram(ax4,onfraction(newtrace&numberOfBlinks>1),"DisplayStyle","stairs",Normalization="probability")
                xlabel(ax4,'on fraction')
                hold(ax4,'on');
                
                histogram(ax5,numberOfBlinks(newtrace)-1,0:quantile(numberOfBlinks(newtrace)-1,.99),"DisplayStyle","stairs",Normalization="probability")
                xlabel(ax5,'number of blinks')
                ax5.YScale=p.yscale.selection;ax5.XLim(1)=0;hold(ax5,'on');
                
                
                
                %%
                dT=[0 1 3 10 30 100 300 1000 3000 ];
                % dT=[0 1 10  100  1000 10000];
                % dT=0
                locd2=obj.locData.copy;
                thisf=locd2.loc.filenumber==files(f);
                locd2.removelocs(~thisf);
                ax8=obj.initaxis(['photdT' num2str(f)]);
                ax9=obj.initaxis(['ondT' num2str(f)]);
                
                clear mp dTtxt mon
                for k=1:length(dT)
                    disp(['file ' num2str(f) ', dT ' num2str(dT(k))])
                    locd2.regroup([30 100],dT(k),2);
                    
                    pmax=15000;
                    histogram(ax8,locd2.grouploc.phot,0:100:pmax,"DisplayStyle","stairs",Normalization="probability")
                    xlabel(ax8,'photons')
                    hold(ax8,'on');
                    xlim(ax8,[0 pmax])
                    ax8.YScale=p.yscale.selection;ax8.XLim(1)=0;
                    
                    histogram(ax9,locd2.grouploc.numberInGroup,1:50,"DisplayStyle","stairs",Normalization="probability")
                    xlabel(ax9,'on-time (frames)')
                    hold(ax9,'on')
                    dTtxt{k}=num2str(dT(k));
                    mp(k)=mean(locd2.grouploc.phot);
                    mon(k)=mean(locd2.grouploc.numberInGroup);
                    ax9.YScale=p.yscale.selection;ax9.XLim(1)=0;
                end
                legend(ax8,dTtxt)
                legend(ax9,dTtxt)
               
                plot(ax6,dT,mp)
                ax6.XScale=p.yscale.selection;
                xlabel(ax6,'dT (frames)')
                ylabel(ax6,'mean photon numbers')
                hold(ax6,'on');
                
                plot(ax7,dT,mon)
                ax7.XScale=p.yscale.selection;
                xlabel(ax7,'dT (frames)')
                ylabel(ax7,'mean on time (frames)')
                hold(ax7,'on');
            end

            legend(ax1,leg);legend(ax2,leg);legend(ax3,leg);legend(ax4,leg);legend(ax5,leg);
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end


function pard=guidef(obj)

pard.maxdtt.object=struct('String','maximum off time (s)','Style','text');
pard.maxdtt.position=[1,1];
pard.maxdtt.Width=1.5;
pard.maxdt.object=struct('String','5','Style','edit');
pard.maxdt.position=[1,2.5];
pard.maxdt.Width=0.5;
pard.showwhat.object=struct('String',{{'layers','all files'}},'Style','popupmenu');
pard.showwhat.position=[1,3];
pard.showwhat.Width=1;
% 
pard.yscalet.object=struct('String','scale y-axis','Style','text');
pard.yscalet.position=[2,1];
pard.yscale.object=struct('String',{{'linear','log'}},'Style','popupmenu');
pard.yscale.position=[2,2];
pard.plugininfo.description=sprintf('');
end