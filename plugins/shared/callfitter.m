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
%XXXXXX

fh=@CPUmleFit_LM_MultiChannel2026;
%XXXX
args=varargin;

% k_block=100;
imfit=args{1};
args2=args(2:end);
numlocs=size(imfit,3);
% parpool('local');                          % processes
% imfitc = parallel.pool.Constant(imfit);
imfitc=imfit;
% starts = 1:k_block:size(imfit,3);
% outparnum=3;

% [varargout0{1:nargout}] = fh(imfitc.Value(:,:,1),args{2:end});
[varargout0{1:nargout}] = fh(imfitc(:,:,1,:),args{2:end});
% for k=nargout:-1:1
    out1=zeros(size(imfit,3),size(varargout0{1},2));
    out2=zeros(size(imfit,3),size(varargout0{2},2));
    out3=zeros(size(imfit,3),size(varargout0{3},2));
% end

parfor b = 1:numlocs
    % s = starts(b);
    % e = min(s + k_block - 1, size(imfitc.Value,3));
    % Ablk = imfitc.Value(:,:,b);
    Ablk = imfitc(:,:,b,:);

    [t1,t2,t3] = fh(Ablk,args2{:});

    out1(b,:)=t1;
    out2(b,:)=t2;
    out3(b,:)=t3;

    % out(:,:,s:e) = yblk;
end
varargout={out1,out2,out3};
% [varargout{1:nargout}] = feval(fh, varargin{:});

end


