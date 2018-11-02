classdef CorrectUnconnectbyConnect<interfaces.DialogProcessor
    methods
        function obj=CorrectUnconnectbyConnect(varargin)     
            obj@interfaces.DialogProcessor(varargin{:}) ;  
        end
        
        function out=run(obj,p) 
            out=[];
            lu=obj.locData.loc;
            lg=obj.locData.grouploc;
            
            [gisort,inds]=sort(lu.groupindex);
%             idx1=1;
            idx2=1;
            for k=1:length(lg.groupindex);
                while idx2<length(gisort)&&gisort(idx2)<lg.groupindex(k)
                    idx2=idx2+1;
                end
                while idx2<length(gisort)&&gisort(idx2)==lg.groupindex(k)
                    obj.locData.loc.xnm(inds(idx2))=lg.xnm(k);
                    obj.locData.loc.ynm(inds(idx2))=lg.ynm(k);
                    if isfield(obj.locData.loc,'znm')
                        obj.locData.loc.znm(inds(idx2))=lg.znm(k);
                    end
                    if isfield(obj.locData.loc,'locprecznm')
                        obj.locData.loc.locprecznm(inds(idx2))=lg.locprecznm(k);
                    end
                    obj.locData.loc.locprecnm(inds(idx2))=lg.locprecnm(k);
                    idx2=idx2+1;
                end
            end
                
             
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef


% 
% pard.textb.object=struct('String','Channel','Style','text');
% pard.textb.position=[1,1];
% pard.channel.object=struct('String','1','Style','edit');
% pard.channel.position=[1,2];
% 


% pard.connectmode.object=struct('String','connect->unconnect|unconnect->connect','Style','popupmenu','Value',1);
% pard.connectmode.position=[2,1];
pard.plugininfo.type='ProcessorPlugin';
end