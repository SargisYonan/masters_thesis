function S = variogram(x,y,varargin)

% isotropic and anisotropic experimental (semi-)variogram
%
% Syntax:
%   d = variogram(x,y)
%   d = variogram(x,y,'propertyname','propertyvalue',...)
%
% Description:
%   variogram calculates the experimental variogram in various 
%   dimensions. 
%
% Input:
%   x - array with coordinates. Each row is a location in a 
%       size(x,2)-dimensional space (e.g. [x y elevation])
%   y - column vector with values of the locations in x. 
%
% Propertyname/-value pairs:
%   nrbins - number bins the distance should be grouped into
%            (default = 20)
%   maxdist - maximum distance for variogram calculation
%            (default = maximum distance in the dataset / 2)
%   type -   'gamma' returns the variogram value (default)
%            'binnedcloud' returns the binned variogram cloud
%            'cloud' returns the variogram cloud
%   plot   - true -> plot variogram
%            false -> don't plot (default)
%   subsample - number of randomly drawn points if large datasets are used.
%               scalar (positive integer, e.g. 3000)
%               inf (default) = no subsampling
%   anisotropy - false (default), true (works only in two dimensions)
%   thetastep - if anisotropy is set to true, specifying thetastep 
%            allows you the angle width (default 30°)
%   
%   
% Output:
%   d - structure array with distance and gamma - vector
%   
% Example: Generate a random field with periodic variation in x direction
% 
%     x = rand(1000,1)*4-2;  
%     y = rand(1000,1)*4-2;
%     z = 3*sin(x*15)+ randn(size(x));
%
%     subplot(2,2,1)
%     scatter(x,y,4,z,'filled'); box on;
%     ylabel('y'); xlabel('x')
%     title('data (coloring according to z-value)')
%     subplot(2,2,2)
%     hist(z,20)
%     ylabel('frequency'); xlabel('z')
%     title('histogram of z-values')
%     subplot(2,2,3)
%     d = variogram([x y],z,'plotit',true,'nrbins',50);
%     title('Isotropic variogram')
%     subplot(2,2,4)
%     d2 = variogram([x y],z,'plotit',true,'nrbins',50,'anisotropy',true);
%     title('Anisotropic variogram')
%
%
% See also: KRIGING, VARIOGRAMFIT
%
% Date: 8. August, 2014
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)


% extent of dataset
minx   = min(x,[],1);
maxx   = max(x,[],1);
maxd   = sqrt(sum((maxx-minx).^2));
nrdims = size(x,2);

% Parse inputs
p = inputParser;
p.FunctionName = 'variogram';

addRequired(p,'x',@(x) validateattributes(x,{'numeric'},{}));
addRequired(p,'y',@(y) validateattributes(y,{'numeric'},{'column','nrows',size(x,1)}));

addParamValue(p,'nrbins',20,@(x) validateattributes(x,{'numeric'},{'scalar','integer','>',0}));
addParamValue(p,'maxdist',maxd/2,@(x) validateattributes(x,{'numeric'},{'scalar','>',0}));
addParamValue(p,'type','gamma',@(x) ischar(validatestring(x,{'gamma','binnedcloud1','cloud'})));
addParamValue(p,'plot',true,@(x) isscalar(x));
addParamValue(p,'subsample',inf,@(x) validateattributes(x,{'numeric'},{'scalar','integer','>',0}));
addParamValue(p,'anisotropy',false,@(x) isscalar(x));
addParamValue(p,'thetastep',30,@(x) validateattributes(x,{'numeric'},{'scalar','>',0}));

parse(p,x,y,varargin{:});

% convert inputParser class to struct to allow values to be modified
% further
p = p.Results;

% check for nans
II      = any(isnan(x),2) | isnan(y);
x(II,:) = [];
y(II)   = [];

% check maximum distance
if p.maxdist > maxd;
    warning('Matlab:Variogram',...
            ['Maximum distance exceeds maximum distance \n' ... 
             'in the dataset. maxdist was decreased to ' num2str(maxd) ]);
    p.maxdist  = maxd;
end

% anisotropy
if p.anisotropy && nrdims ~= 2 
    p.anisotropy = false;
    warning('Matlab:Variogram',...
            'Anistropy is only supported for 2D data');
end

% take only a subset of the data;
if ~isinf(p.subsample) && numel(y)>p.subsample;
    IX = randperm(numel(y),p.subsample);
    x  = x(IX,:);
    y  = y(IX,:);
end

% calculate bin tolerance
tol = p.maxdist/p.nrbins;

% calculate distance matrix
iid = distmat(x,p.maxdist);

% calculate squared difference between values of coordinate pairs
lam      = (y(iid(:,1))-y(iid(:,2))).^2;

% anisotropy
if p.anisotropy 
    nrthetaedges = floor(180/(p.thetastep))+1;
  
    % calculate with radians, not degrees
    p.thetastep = p.thetastep/180*pi;

    % calculate angles, note that angle is calculated clockwise from top
    theta    = atan2(x(iid(:,2),1)-x(iid(:,1),1),...
                     x(iid(:,2),2)-x(iid(:,1),2));
    
    % only the semicircle is necessary for the directions
    I        = theta < 0;
    theta(I) = theta(I)+pi;
    I        = theta >= pi-p.thetastep/2;
    theta(I) = 0;
        
    % create a vector with edges for binning of theta
    % directions go from 0 to 180 degrees;
    thetaedges = linspace(-p.thetastep/2,pi-p.thetastep/2,nrthetaedges);
    
    % bin theta
    [ntheta,ixtheta] = histc(theta,thetaedges);
    
    % bin centers
    thetacents = thetaedges(1:end)+p.thetastep/2;
    thetacents(end) = pi; %[];
end

% calculate variogram
switch p.type
    case 'gamma'
        % variogram anonymous function
        fvar     = @(x) 1./2 * mean(x);
        
        % distance bins
        edges      = linspace(0,p.maxdist,p.nrbins+1);
        edges(end) = inf;

        [nedge,ixedge] = histc(iid(:,3),edges);
        
        if p.anisotropy
            S.val      = accumarray([ixedge ixtheta],lam,...
                                 [numel(edges) numel(thetaedges)],fvar,nan);
            S.val(:,end)=S.val(:,1); 
            S.theta    = thetacents;
            S.num      = accumarray([ixedge ixtheta],ones(size(lam)),...
                                 [numel(edges) numel(thetaedges)],@sum,nan);
            S.num(:,end)=S.num(:,1);                 
        else
            S.val      = accumarray(ixedge,lam,[numel(edges) 1],fvar,nan);     
            S.num      = accumarray(ixedge,ones(size(lam)),[numel(edges) 1],@sum,nan);
        end
        S.distance = (edges(1:end-1)+tol/2)';
        S.val(end,:) = [];
        S.num(end,:) = [];

    case 'binnedcloud'
        edges      = linspace(0,p.maxdist,p.nrbins+1);
        edges(end) = inf;
        
        [nedge,ixedge] = histc(iid(:,3),edges);
        
        S.distance = edges(ixedge) + tol/2;
        S.distance = S.distance(:);
        S.val      = lam;  
        if p.anisotropy            
            S.theta   = thetacents(ixtheta);
        end
    case 'cloud'
        S.distance = iid(:,3);
        S.val      = lam;
        if p.anisotropy            
            S.theta   = thetacents(ixtheta);
        end
end


% create plot if desired
if p.plot
    switch p.type
        case {'default','gamma'}
            marker = 'o--';
        otherwise
            marker = '.';
    end
    
    if ~p.anisotropy
        plot(S.distance,S.val,marker);
        axis([0 p.maxdist 0 max(S.val)*1.1]);
        xlabel('h');
        ylabel('\gamma (h)');
        title('(Semi-)Variogram');
    else
        [Xi,Yi] = pol2cart(repmat(S.theta,numel(S.distance),1),repmat(S.distance,1,numel(S.theta)));
        surf(Xi,Yi,S.val)
        xlabel('h y-direction')
        ylabel('h x-direction')
        zlabel('\gamma (h)')
        title('directional variogram')
%         set(gca,'DataAspectRatio',[1 1 1/30])
    end
end
        
end


% subfunction distmat

function iid = distmat(X,dmax)

% constrained distance function
%
% iid -> [rows, columns, distance]
 

n     = size(X,1);
nrdim = size(X,2);
if size(X,1) < 1000;
    [i,j] = find(triu(true(n)));
    if nrdim == 1;
        d = abs(X(i)-X(j));
    elseif nrdim == 2;
        d = hypot(X(i,1)-X(j,1),X(i,2)-X(j,2));
    else
        d = sqrt(sum((X(i,:)-X(j,:)).^2,2));
    end
    I = d<=dmax;
    iid = [i(I) j(I) d(I)];
else
    ix = (1:n)';
    if nrdim == 1;
        iid = arrayfun(@distmatsub1d,(1:n)','UniformOutput',false);
    elseif nrdim == 2;
        % if needed change distmatsub to distmatsub2d which is numerically
        % better but slower
        iid = arrayfun(@distmatsub,(1:n)','UniformOutput',false);
    else
        iid = arrayfun(@distmatsub,(1:n)','UniformOutput',false);
    end
    nn  = cellfun(@(x) size(x,1),iid,'UniformOutput',true);  
    I   = nn>0;
    ix  = ix(I);
    nn  = nn(I);
    nncum = cumsum(nn);
    c     = zeros(nncum(end),1);
    c([1;nncum(1:end-1)+1]) = 1;
    i = ix(cumsum(c));
    iid = [i cell2mat(iid)];
    
end

function iid = distmatsub1d(i) 
    j  = (i+1:n)'; 
    d  = abs(X(i)-X(j));
    I  = d<=dmax;
    iid = [j(I) d(I)];
end

function iid = distmatsub2d(i)  %#ok<DEFNU>
    j  = (i+1:n)'; 
    d = hypot(X(i,1) - X(j,1),X(i,2) - X(j,2));
    I  = d<=dmax;
    iid = [j(I) d(I)];
end
    
function iid = distmatsub(i)
    j  = (i+1:n)'; 
    d = sqrt(sum(bsxfun(@minus,X(i,:),X(j,:)).^2,2));
    I  = d<=dmax;
    iid = [j(I) d(I)];
end
end
