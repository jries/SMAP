function varargout = polylineEqDist(inSpots, varargin)
% Transform user defined positions along the polyline into equidistant spots 
% inSpots = set of coordinates

defaultMethod = 'linear';
expectedMethod={'linear', 'spline', 'pchip', 'csape'};

p = inputParser;
p.addParameter('method',defaultMethod,@(x) any(validatestring(x,expectedMethod)));
p.parse(varargin{:});
temp=size(inSpots);

try
    [varargout{1:nargout}] = interparc(temp(1),inSpots(:,1),inSpots(:,2),p.Results.method);
catch ME
    warning('Error in transforming user specified polyline into equidistant spots!');
    print(ME);

end
end

