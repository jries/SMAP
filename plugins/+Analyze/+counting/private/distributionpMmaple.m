
% Np=2;

poiss=1;
maturation={0.5, 0.75, 0.5, 0.75};
% maturation=0.75;
form={'bx','b+','ro','rd'}
nblink={0,0,.1,.1};
nblink={0};
rep=25000;

Np=[ 2 3 4 5 7 10 15 20 30 40 50 70 100 150 200 300 400 500 700 1000];
% Np=[10]
% maturation={.15}
% nblink={.1}
% nblink={8}
bfac=1.2;


figure(4)
hold off;


stdn=[];
mn=[];

for c=1:length(nblink)

for n=1:length(Np)
n
loc=zeros(rep,1);
for m=1:rep
    for l=1:Np(n)
    if rand<maturation{c}

    if nblink{c}>0
        if poiss
            nb=poissrnd(nblink{c});
        else
            
            nb=0;
            while rand>= nblink{c}
                nb=nb+1;

            end
        end
    else
        nb=1;
    end
      loc(m)=loc(m)+ nb; %exponential distribution: bleaching
       %loc(m)=loc(m)+ ceil(exprnd(nblink)); %exponential distribution: bleaching, minimum 1
      %loc(m)=loc(m)+ 1;  %activation only once
%       loc(m)=loc(m)+ round(exprnd(nblink{c}))+1; %exponential distribution: bleaching
    end

    end
end
% loc(loc==0)=[];

figure(c+5)
% subplot(4,5,n)
% hist(loc,0:1:max(loc));
[h,x]=hist(loc,0:1:max(loc));


 y=brightnessdistribution(x',Np(n),maturation{c},'blinkmode','NoBlink');

plot(x,h/sum(h(:)));
hold on
plot(x,y,'r')
hold off
title(std(loc)/mean(loc))


stdn(n)=std(loc);
mn(n)=mean(loc);
end
figure(1)
semilogx(Np,stdn./mn,form{c})

hold off
end