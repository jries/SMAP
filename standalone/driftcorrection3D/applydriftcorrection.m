function poso=applydriftcorrection(drift,pos)
%pos: localizations with fields xnm ynm znm frame
indg=pos.frame>0;

if isfield(drift,'xy')
    poso.xnm=pos.xnm;
    poso.ynm=pos.ynm;
%     for k=1:length(drit.xy)
        if isfield(drift.xy(1),'x') 
            if strcmp(drift.xy(1).mirror,'horizontal')&&length(drift.xy)>1
                indpos=pos.xnm<drift.xy(1).midpoint;
                poso.xnm(indg&indpos)=pos.xnm(indg&indpos)-drift.xy(1).x(pos.frame(indg&indpos));
                poso.xnm(indg&~indpos)=pos.xnm(indg&~indpos)-drift.xy(2).x(pos.frame(indg&~indpos)); 
            else
                poso.xnm(indg)=pos.xnm(indg)-drift.xy(1).x(pos.frame(indg));
            end
        end
        if isfield(drift.xy(1),'y')
            if strcmp(drift.xy(1).mirror,'vertical')&&length(drift.xy)>1
                indpos=pos.ynm<drift.xy(2).midpoint;
                poso.ynm(indg&indpos)=pos.ynm(indg&indpos)-drift.xy(1).y(pos.frame(indg&indpos));
                poso.ynm(indg&~indpos)=pos.ynm(indg&~indpos)-drift.xy(2).y(pos.frame(indg&~indpos)); 
            else
                poso.ynm(indg)=pos.ynm(indg)-drift.xy(1).y(pos.frame(indg));
            end
        end
end

if isfield(drift,'z')&&isfield(pos,'znm')
    poso.znm=pos.znm;
    poso.znm(indg)=pos.znm(indg)-drift.z(pos.frame(indg));
end