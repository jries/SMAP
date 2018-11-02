function avout=fullfilt(in,window)
avout=zeros(length(in),1);
w2=round(window/2);
av=filter(ones(1,2*w2+1)/(2*w2+1),1,in);
avout(1:w2)=av(2*w2+1);

% length(w2+1:w2+length(av)-2*w2)
% length(2*w2+1:length(av))
avout(w2+1:w2+length(av)-2*w2)=av(2*w2+1:length(av));
avout(length(av)-w2:length(av))=av(length(av));

