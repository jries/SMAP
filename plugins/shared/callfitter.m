function varargout=callfitter(fh,varargin)
if contains(func2str(fh),'GPU')
    [varargout{1:nargout}] = feval(fh, varargin{:});
    return
end
%CPU fitter
%P,CRLB,LogL

% p = gcp('nocreate');
% if isempty(p)
%     parpool('threads')
% end


% XXXX bugs in this code, slight difference to direct calling. NO idea why

args=varargin;

% k_block=100;
imfit=args{1};

numlocs=size(imfit,3);
% parpool('local');                          % processes
% imfitc = parallel.pool.Constant(imfit);
% imfitc=imfit;
% starts = 1:k_block:size(imfit,3);
% outparnum=3;

% [varargout0{1:nargout}] = fh(imfitc.Value(:,:,1),args{2:end});
[varargout0{1:nargout}] = fh(imfit(:,:,1,:),args{2:end});
% for k=nargout:-1:1
    out1=zeros(numlocs,size(varargout0{1},2),'single');
    out2=zeros(numlocs,size(varargout0{2},2),'single');
    out3=zeros(numlocs,size(varargout0{3},2),'single');
% end

fhs=func2str(fh);
if contains(fhs,"MultiChannel")
    mode=1;
    % arglocs=[1,3,6,9];
    % locdim=zeros(9,1)-1;
    % locdim(arglocs)=[3,2,3,2];
else
    error('not implemented')
end
if mode==1
    for b=numlocs:-1:1
        argsh{b}={args{1}(:,:,b,:),args{2},args{3}(:,b),args{4},args{5},args{6}(:,:,b),args{7},args{8},args{9}(1,b)};
    end
end

parfor b = 1:numlocs

    [t1,t2,t3] = feval(fh,argsh{b}{:});
    % [t1,t2,t3] = fh(imfit(:,:,b,:),args{2:end});

    out1(b,:)=t1;
    out2(b,:)=t2;
    out3(b,:)=t3;
end
varargout={out1,out2,out3};
end


