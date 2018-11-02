function y=mygaussforfit(param,x)
%a x sigma bg
y=param(1)*exp(-(x-param(2)).^2/(2*param(3)^2))+param(4);