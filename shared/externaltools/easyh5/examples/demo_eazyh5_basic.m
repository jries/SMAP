%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Demonstration of Basic Utilities of EasyH5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rngstate = rand ('state');
randseed=hex2dec('623F9A9E');
clear data2hdf h52data

opt.releaseid=0;
vers=ver('MATLAB');
if(~isempty(vers))
    opt.releaseid=datenum(vers(1).Date);
end
opt.skipempty=(opt.releaseid<datenum('1-Jan-2015'));

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  a simple scalar value \n')
fprintf(1,'%%=================================================\n\n')

data2hdf=pi
saveh5(data2hdf,'test.h5');
h52data=loadh5('test.h5')
isequaln(data2hdf,h52data.data2hdf)

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  an empty array \n')
fprintf(1,'%%=================================================\n\n')

data2hdf=[]
saveh5(data2hdf,'test.h5');
h52data=loadh5('test.h5')
if(~opt.skipempty)
    isequaln(data2hdf,h52data.data2hdf)
end

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  an ampty string \n')
fprintf(1,'%%=================================================\n\n')

data2hdf=''
saveh5(data2hdf,'test.h5')
h52data=loadh5('test.h5')
if(~opt.skipempty)
    isequaln(data2hdf,h52data.data2hdf)
end

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  a simple row vector \n')
fprintf(1,'%%=================================================\n\n')

data2hdf=1:3
saveh5(data2hdf,'test.h5');
h52data=loadh5('test.h5')
isequaln(data2hdf,h52data.data2hdf)

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  a simple column vector \n')
fprintf(1,'%%=================================================\n\n')

data2hdf=(1:3)'
saveh5(data2hdf,'test.h5');
h52data=loadh5('test.h5')
isequaln(data2hdf,h52data.data2hdf)

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  a string array \n')
fprintf(1,'%%=================================================\n\n')

data2hdf=['AC';'EG']
saveh5(data2hdf,'test.h5');
h52data=loadh5('test.h5')
isequaln(data2hdf,h52data.data2hdf)

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  a string with escape symbols \n')
fprintf(1,'%%=================================================\n\n')

data2hdf=sprintf('AB\tCD\none"two')
saveh5(data2hdf,'test.h5')
h52data=loadh5('test.h5')
isequaln(data2hdf,h52data.data2hdf)

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  a mix-typed cell \n')
fprintf(1,'%%=================================================\n\n')

data2hdf={'a',true,[2;3]}
saveh5(data2hdf,'test.h5');
h52data=loadh5('test.h5','regroup',1)
isequaln(data2hdf,h52data.data2hdf)

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  a 3-D array in nested array form\n')
fprintf(1,'%%=================================================\n\n')

data2hdf=reshape(1:(2*4*6),[2,4,6]);
saveh5(data2hdf,'test.h5')  % nestarray for 4-D or above is not working
h52data=loadh5('test.h5')
isequaln(data2hdf,h52data.data2hdf)

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  a 4-D array in annotated array form\n')
fprintf(1,'%%=================================================\n\n')

data2hdf=reshape(1:(2*4*3*2),[2,4,3,2]);
saveh5(data2hdf,'test.h5')  % nestarray for 4-D or above is not working
h52data=loadh5('test.h5')
isequaln(data2hdf,h52data.data2hdf)

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  a 3-D array in annotated array form (EasyH5 1.9 or earlier)\n')
fprintf(1,'%%=================================================\n\n')

data2hdf=reshape(1:(2*4*6),[2,4,6]);
saveh5(data2hdf,'test.h5')  % nestarray for 4-D or above is not working
h52data=loadh5('test.h5')
isequaln(data2hdf,h52data.data2hdf)


fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  a complex number\n')
fprintf(1,'%%=================================================\n\n')

data2hdf=1+2i
saveh5(data2hdf,'test.h5');
h52data=loadh5('test.h5') 
isequaln(data2hdf,h52data.data2hdf)

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  a complex matrix\n')
fprintf(1,'%%=================================================\n\n')

data2hdf=magic(6);
data2hdf=data2hdf(:,1:3)+data2hdf(:,4:6)*1i
saveh5(data2hdf,'test.h5');
h52data=loadh5('test.h5')
isequaln(data2hdf,h52data.data2hdf)

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  MATLAB special constants\n')
fprintf(1,'%%=================================================\n\n')

data2hdf=[NaN Inf -Inf]
saveh5(data2hdf,'test.h5')  % nestarray for 4-D or above is not working
h52data=loadh5('test.h5')
isequaln(data2hdf,h52data.data2hdf)

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  a real sparse matrix\n')
fprintf(1,'%%=================================================\n\n')

data2hdf=sprand(10,10,0.1)
saveh5(data2hdf,'test.h5')  % nestarray for 4-D or above is not working
h52data=loadh5('test.h5')
isequaln(data2hdf,h52data.data2hdf)

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  a complex sparse matrix\n')
fprintf(1,'%%=================================================\n\n')

data2hdf=data2hdf-data2hdf*1i
saveh5(data2hdf,'test.h5')  % nestarray for 4-D or above is not working
h52data=loadh5('test.h5')
isequaln(data2hdf,h52data.data2hdf)

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  an all-zero sparse matrix\n')
fprintf(1,'%%=================================================\n\n')

data2hdf=sparse(2,3);
saveh5(data2hdf,'test.h5')  % nestarray for 4-D or above is not working
h52data=loadh5('test.h5')
if(~opt.skipempty)
    isequaln(data2hdf,h52data.data2hdf)
end

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  an empty sparse matrix\n')
fprintf(1,'%%=================================================\n\n')

data2hdf=sparse([]);
saveh5(data2hdf,'test.h5')  % nestarray for 4-D or above is not working
h52data=loadh5('test.h5')
if(~opt.skipempty)
    isequaln(data2hdf,h52data.data2hdf)
end

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  an empty 0-by-0 real matrix\n')
fprintf(1,'%%=================================================\n\n')

data2hdf=[];
saveh5(data2hdf,'test.h5')  % nestarray for 4-D or above is not working
h52data=loadh5('test.h5')
if(~opt.skipempty)
    isequaln(data2hdf,h52data.data2hdf)
end

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  an empty 0-by-3 real matrix\n')
fprintf(1,'%%=================================================\n\n')

data2hdf=zeros(0,3);
saveh5(data2hdf,'test.h5')  % nestarray for 4-D or above is not working
h52data=loadh5('test.h5')
if(~opt.skipempty)
    isequaln(data2hdf,h52data.data2hdf)
end

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  a sparse real column vector\n')
fprintf(1,'%%=================================================\n\n')

data2hdf=sparse([0,3,0,1,4]');
saveh5(data2hdf,'test.h5')  % nestarray for 4-D or above is not working
h52data=loadh5('test.h5')
isequaln(data2hdf,h52data.data2hdf)

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  a sparse complex column vector\n')
fprintf(1,'%%=================================================\n\n')

data2hdf=data2hdf-1i*data2hdf;
saveh5(data2hdf,'test.h5')  % nestarray for 4-D or above is not working
h52data=loadh5('test.h5')
isequaln(data2hdf,h52data.data2hdf)

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  a sparse real row vector\n')
fprintf(1,'%%=================================================\n\n')

data2hdf=sparse([0,3,0,1,4]);
saveh5(data2hdf,'test.h5')  % nestarray for 4-D or above is not working
h52data=loadh5('test.h5')
isequaln(data2hdf,h52data.data2hdf)

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  a sparse complex row vector\n')
fprintf(1,'%%=================================================\n\n')

data2hdf=data2hdf-1i*data2hdf;
saveh5(data2hdf,'test.h5')  % nestarray for 4-D or above is not working
h52data=loadh5('test.h5')
isequaln(data2hdf,h52data.data2hdf)

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  a structure\n')
fprintf(1,'%%=================================================\n\n')

data2hdf=struct('name','Think Different','year',1997,'magic',magic(3),...
                 'misfits',[Inf,NaN],'embedded',struct('left',true,'right',false))
saveh5(data2hdf,'test.h5')  % nestarray for 4-D or above is not working
h52data=loadh5('test.h5','regroup',1)
isequaln(data2hdf,h52data.data2hdf)

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  a structure array\n')
fprintf(1,'%%=================================================\n\n')

data2hdf=struct('name','Nexus Prime','rank',9);
data2hdf(2)=struct('name','Sentinel Prime','rank',9);
data2hdf(3)=struct('name','Optimus Prime','rank',9);
saveh5(data2hdf,'test.h5')  % nestarray for 4-D or above is not working
h52data=loadh5('test.h5','regroup',1)
isequaln(data2hdf,h52data.data2hdf)

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  a cell array\n')
fprintf(1,'%%=================================================\n\n')

data2hdf=cell(3,1);
data2hdf{1}=struct('buzz',1.1,'rex',1.2,'bo',1.3,'hamm',2.0,'slink',2.1,'potato',2.2,...
              'woody',3.0,'sarge',3.1,'etch',4.0,'lenny',5.0,'squeeze',6.0,'wheezy',7.0);
data2hdf{2}=struct('Ubuntu',['Kubuntu';'Xubuntu';'Lubuntu']);
data2hdf{3}=[10.04,10.10,11.04,11.10]
saveh5(data2hdf,'test.h5')  % nestarray for 4-D or above is not working
h52data=loadh5('test.h5','regroup',1)
isequaln(data2hdf,h52data.data2hdf)


fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  a function handle\n')
fprintf(1,'%%=================================================\n\n')

data2hdf=@(x) x+1
saveh5(data2hdf,'test.h5')  % nestarray for 4-D or above is not working
h52data=loadh5('test.h5')
isequaln(data2hdf,h52data.data2hdf)

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  a 2D cell array\n')
fprintf(1,'%%=================================================\n\n')

data2hdf={{1,{2,3}},{4,5},{6};{7},{8,9},{10}};
saveh5(data2hdf,'test.h5')  % nestarray for 4-D or above is not working
h52data=loadh5('test.h5','regroup',1)  % only saveh5 works for cell arrays, loadh5 has issues
%isequaln(data2hdf,h52data.data2hdf)

fprintf(1,'\n%%=================================================\n')
fprintf(1,'%%  a 2D struct array\n')
fprintf(1,'%%=================================================\n\n')

data2hdf=repmat(struct('idx',0,'data','structs'),[2,3])
for i=1:6
    data2hdf(i).idx=i;
end
saveh5(data2hdf,'test.h5')  % nestarray for 4-D or above is not working
h52data=loadh5('test.h5','regroup',1)
isequaln(data2hdf,h52data.data2hdf)


if(exist('datetime'))
    fprintf(1,'\n%%=================================================\n')
    fprintf(1,'%%  datetime object \n')
    fprintf(1,'%%=================================================\n\n')

    data2hdf=datetime({'8 April 2015','9 May 2015'}, 'InputFormat','d MMMM yyyy')
    saveh5(data2hdf,'test.h5')  % nestarray for 4-D or above is not working
    h52data=loadh5('test.h5')
    isequaln(data2hdf,h52data.data2hdf)
end

if(exist('containers.Map'))
    fprintf(1,'\n%%=================================================\n')
    fprintf(1,'%%  a container.Maps object \n')
    fprintf(1,'%%=================================================\n\n')

    data2hdf=containers.Map({'Andy','William','Om'},[21,21,22])
    saveh5(data2hdf,'test.h5');
    h52data=loadh5('test.h5')
    isequaln(data2hdf,h52data.data2hdf)
end

if(exist('istable'))
    fprintf(1,'\n%%=================================================\n')
    fprintf(1,'%%  a table object \n')
    fprintf(1,'%%=================================================\n\n')

    Names={'Andy','William','Om'}';
    Age=[21,21,22]';
    data2hdf=table(Names,Age)
    saveh5(data2hdf,'test.h5')  % nestarray for 4-D or above is not working
    h52data=loadh5('test.h5')
    isequaln(data2hdf,h52data.data2hdf)
end

try
    val=zlibencode('test');
    fprintf(1,'\n%%=================================================\n')
    fprintf(1,'%%  a 2-D array in compressed array format\n')
    fprintf(1,'%%=================================================\n\n')

    data2hdf=eye(10);
    data2hdf(20,1)=1;
    saveh5(data2hdf,'test.h5','compression','deflate')  % nestarray for 4-D or above is not working
    h52data=loadh5('test.h5')
    isequaln(data2hdf,h52data.data2hdf)
catch
end

rand ('state',rngstate);

