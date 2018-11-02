function h=brightnessdistribution(k,varargin)
% (k,N0,pmature,pblink,c,'blinkmode',bm, 'overlapping',ol)
%k number of events.
% function h=brightnessdistribution(k,N0,pmature,pblink,mode)
global  kx
%input checking
p=inputParser;
defaultMode='Exponential';
validateMode={'Exponential','Poisson','Noblink'};
checkMode = @(x) any(validatestring(x,validateMode));

% defaultc='single';
% validatec={'single','multi'};
% checkc = @(x) any(validatestring(x,validatec));

addRequired(p,'N0',@isnumeric)
addRequired(p,'pmature',@isnumeric)

addOptional(p,'pblink',0,@isnumeric)
addOptional(p,'ncluster',0,@isnumeric)
addParameter(p,'blinkmode',defaultMode,checkMode)
% addParameter(p,'overlapping',defaultc,checkc)

p.KeepUnmatched = true;
parse(p,varargin{:})

par=p.Results;

if par.ncluster<=0
    h=Pcluster(k,par.N0,par.pmature,par.pblink,par.blinkmode);
elseif par.ncluster<3%up to 3
    
    PP11=PP1(1,par.ncluster);    
    pref=Pcluster(k,par.N0,par.pmature,par.pblink,par.blinkmode);
    h=pref.*PP11;
    PP12=PP1(2,par.ncluster);

    for k=0:length(h)-1
        for a1=0:k
            h(k+1)=h(k+1)+PP12.*pref(a1+1).*pref(k-a1+1);
%             h(k)=h(k)+PP1(2,par.ncluster).*Pcluster(a1,par.N0,par.pmature,par.pblink,par.blinkmode).*Pcluster(k-a1,par.N0,par.pmature,par.pblink,par.blinkmode);
        end
    end
    PP13=PP1(3,par.ncluster);
    for k=0:length(h)-1
        for a1=0:k
            for a2=0:k-a1
                h(k+1)=h(k+1)+PP13.*pref(a1+1).*pref(a2+1).*pref(k-a1-a2+1);
            end
        end
    end

else
    disp('ncluster too high for approx with 2')
    
end


function h=Pcluster(k,N0,pmature,pblink,blinkmode)
global sumk
if strcmp(blinkmode,'Poisson')
    h=zeros(size(k));
    for B=0:N0
        h=h+PP(k,B*pblink).*PBN(B,N0,pmature);

    end
elseif strcmp(blinkmode,'Exponential')&&pblink>0&&pblink<1
    if isempty(sumk)
        sumk=makesumk(1000,500);
    end
    h=zeros(size(k));
    for B=0:N0
        h=h+PkB(k,B,pblink).*PBN(B,N0,pmature);
    end
elseif strcmp(blinkmode,'Noblink')||pblink==0;
     h=PBN(k,N0,pmature);
else
    disp('error')
end

function h=PkB(k,B,pblink)
global sumk
k=round(k);
h=pblink^B*(1-pblink).^k.*sumk(k+1,B+1);

function h=PP1(k,mu)
h=poisspdf(k,mu)/(1-exp(-mu));
% h=exp(-mu)*mu.^k./factorial(k)/(1-exp(-mu));

function h=PP(k,mu)
% h=exp(-mu)*mu.^k./factorial(k);
h=poisspdf(k,mu);
function h=PBN(B,N,pmature)
% if pmature<1
    h=zeros(size(B));
    for b=0:min(length(B)-1,N)
%         h(b+1)=nchoosek(N,B(b+1)).*pmature.^B(b+1).*(1-pmature).^(N-B(b+1));
        h(b+1)=binopdf(B(b+1),N,pmature);
    end
% else
%     h=ones(size(B))/length(B);
% end

function sumk=makesumk(kmax,Bmax)
sumk=zeros(kmax+1,Bmax+1);
k=(0:kmax)';
sumk(:,2)=1;
for b=2:Bmax
    sumk(:,b+1)=sumk(:,b).*(1+k/(b-1));
end
