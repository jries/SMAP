function mysavemovie(img,fsave,varargin)
% varargin: name, value pairs to be passed on to movie constructor
% look at documentation from VideoWriter
% 'profile': standard VideoWriter profile. Default: MPEG-4
% img: frames or array
% fsave: output path

% FrameRate
ind=find(strcmp(varargin,'profile'),1,'first');
if ~isempty(ind)
    profile=varargin{ind+1};
    varargin(ind:ind+1)=[];
else
    profile='MPEG-4';
end

v=VideoWriter(fsave,profile);

for k=1:2:length(varargin)
    set(v,varargin{k},varargin{k+1});
end
if isstruct(img)
    imgp=img;
    %imgp.cdata=pad_to_mod16_rgb(img.cdata);
else
    imgp = pad_to_mod16_rgb(img);
end

open(v);
writeVideo(v,imgp)
close(v);
end


function out = pad_to_mod16_rgb(in)
    [H,W,~,F] = size(in);
    H2 = 16*ceil(H/16);
    W2 = 16*ceil(W/16);
    if H2==H && W2==W
        out = in; return;
    end
    out = zeros(H2, W2, 3, F, 'like', in);
    out(1:H, 1:W, :, :) = in;
    % Replicate last row/col so there's no black border
    if H2 > H
        out(H+1:H2, 1:W, :, :) = repmat(in(H,:,:,:), [H2-H, 1, 1, 1]);
    end
    if W2 > W
        out(1:H2, W+1:W2, :, :) = repmat(out(:,W,:,:), [1, W2-W, 1, 1]);
    end
end
