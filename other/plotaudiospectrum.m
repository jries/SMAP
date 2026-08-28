%plotaudiospectrum
global path
fmin=30;
fmax=800;
Fres=.1; %Hz

[file path]=uigetfile([path filesep  '*.wav']);
[a,Fs]=audioread([path filesep file]);

l=Fs/Fres;
numbl=floor(length(a)/l);
dl=floor(length(a)/numbl/2)*2;
inds=0:dl:length(a);
if length(inds)<2
    dl=floor(length(a)/2)*2;
    inds=[0 dl]; 
    numbl=1;
end
spall=zeros(dl/2+1,1);
for k=1:numbl
 ah=a(inds(k)+1:inds(k+1));
 L=length(ah);
 spectrum=abs(fft(ah)/L);
 f=Fs*(0:(L/2))/L;
 sp1=spectrum(1:L/2+1);
 spall=sp1+spall;
end
 
 figure(89);
 plot(f,spall/numbl);
 xlim([fmin fmax]);