classdef keepfiltered<interfaces.DialogProcessor
    methods
        function obj=keepfiltered(varargin)     
            obj@interfaces.DialogProcessor(varargin{:}) ;  
        end
        
        function out=run(obj,p) 
            out=[];
%            switch p.connectmode.selection
                obj.setPar('undoModule','RemoveLocs');
                notify(obj.P,'backup4undo');
                
                layers=find(obj.getPar('sr_layerson'));
                [l,indg]=obj.locData.getloc({'inungrouped'},'layer',layers(1),'Position','all','removeFilter',{'filenumber'});
                indgood=l.inungrouped;
                for k=2:length(layers)
                     [l,indg]=obj.locData.getloc({'inungrouped'},'layer',layers(k),'Position','all','removeFilter',{'filenumber'});
                     indgood=indgood | l.inungrouped;
                end
                indbad=~indgood;
                obj.locData.removelocs(indbad)
                obj.locData.filter;
                obj.locData.regroup;
%                 
%                 
%                case 'connect->unconnect'
%                    disp('if problems, tell Jonas')
%                    indunc=obj.locData.loc.channel==p.channel;
%                    indcon=obj.locData.grouploc.channel==p.channel;
%                    obj.locData.removelocs(indunc);
%                    loccopy=obj.locData.copy;
%                    loccopy.removelocs(~indcon,'grouploc');
%                    loccopy.loc=loccopy.grouploc;
%                    obj.locData.addLocData(loccopy);
%                    obj.locData.regroup;
% 
%                otherwise 
%                    disp('not implemented')
%            end  
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef



pard.textb.object=struct('String','keeps only the filtered localizations','Style','text');
pard.textb.position=[1,1];
% pard.channel.object=struct('String','1','Style','edit');
% pard.channel.position=[1,2];
% 
% 
% 
% pard.connectmode.object=struct('String','connect->unconnect|unconnect->connect','Style','popupmenu','Value',1);
% pard.connectmode.position=[2,1];
pard.plugininfo.type='ProcessorPlugin';
end