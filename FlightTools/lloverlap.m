function overlapFlag = lloverlap(lat,lon,leg,legs2use,method,radius)
% function overlapFlag = lloverlap(lat,lon,leg,legs2use,method,radius)
% defines area of overlap among multiple flight legs and returns flag for points within that area.
% 
% INPUTS:
% lat: vector latitude
% lon: vector longitude
% leg: index for flight legs. Same length as lat/lon
% legs2use: array specifying which legs to use when looking for overlap.
% method: string specifying how to define the overlap region. Options are
%   'box': area determined by miminum latitude and longitude range spanned by all legs.
%   'circle': area determined by finding all points in a legs that lay within a specified distance
%             from at least one point in every other leg.
% radius: radius of circle over which points will be included. Only required for "circle" method.
%
% OUTPUTS:
% overlapFlag: logical flag, same size as lat/lon, identifying which point are within the area of
%   overlap. If no overlapping points are identified, this will return a flag for all points
%   corresponding to legs within the "legs2use" input.
%
% 20170813 GMW, JRW

% initialize some stuff
n = length(legs2use);
N = length(lat);
ileg = ismember(leg,legs2use); %index for legs to use

switch method
    case 'box'
        % get lat-lon limits to determine box of overlap among all legs
        latrange=nan(n,2); lonrange=nan(n,2);
        for i=1:n
            k = leg==legs2use(i);
            latrange(i,:) = range(lat(k)); %min and max
            lonrange(i,:) = range(lon(k));
        end
        overlapFlag = ileg & ...
            lat>=max(latrange(:,1)) & lat<=min(latrange(:,2)) & ...
            lon>=max(lonrange(:,1)) & lon<=min(lonrange(:,2));
        
    case 'circle'
        % initialize overlap index for each leg
        overlap = cell(n,1);
        for i=1:n
            overlap{i} = true(sum(leg==legs2use(i)),1);
        end
        
        % filter points to be within given radius
        % first loop indexes all points within "radius" of any other point in a leg
        % second loop refines index to points that are all within "radius" of one another
        for loopdeloo = 1:2
            overlap2 = overlap;
            for i=1:(n-1)
                k = leg==legs2use(i);
                lat1 = lat(k); 
                lon1 = lon(k); 
                
                %subsample (only relevant for second loop)
                lat1 = lat1(overlap{i});
                lon1 = lon1(overlap{i});
                len1 = sum(overlap{i});
                
                for j=(i+1):n
                    k = leg==legs2use(j);
                    lat2 = lat(k);
                    lon2 = lon(k);
                    
                    %subsample (only relevant for second loop)
                    lat2 = lat2(overlap{j});
                    lon2 = lon2(overlap{j});
                    len2 = sum(overlap{j});
                    
                    % create all lat-lon pairs
                    lat1big = repmat(lat1,len2,1);
                    lon1big = repmat(lon1,len2,1);
                    lat2big = repmat(lat2',len1,1); lat2big = lat2big(:);
                    lon2big = repmat(lon2',len1,1); lon2big = lon2big(:);
                    
                    % get distances
                    d = lldistkm([lat1big lon1big],[lat2big lon2big]);
                    d = reshape(d,len1,len2);
                    
                    % flag points within radius of one another
                    overlap2{i}(overlap{i}) = overlap2{i}(overlap{i}) & any(d<=radius,2);
                    overlap2{j}(overlap{j}) = overlap2{j}(overlap{j}) & any(d<=radius,1)';
                end
            end
            overlap = overlap2; %refined overlap index
        end
        
        % make big overlap index
        overlapFlag = false(N,1);
        for i=1:n
            overlapFlag(leg==legs2use(i)) = overlap{i};
        end
        
    otherwise
        error(['lloverlap: method "' method '" not recognized. Valid methods are "box" or "circle."'])
end

%if no overlap, assume it wasn't intended that they overlap
if ~sum(overlapFlag)
    overlapFlag = ileg;
end


