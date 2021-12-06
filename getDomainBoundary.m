function [domainBoundary,internalDomain] = getDomainBoundary(pivData)

    % Convert domain into 2D grid
    fluidDomain = array2grid(pivData,pivData.domain);
    
    % Extract boundary and internal domain
    binaryDomain = imbinarize(int8(fluidDomain));
    [boundaryCoords,~] = bwboundaries(binaryDomain,8,'noholes');
    domainBoundary = false(size(fluidDomain));

    for i = 1:length(boundaryCoords{1}(:,1))
        domainBoundary(boundaryCoords{1}(i,1),boundaryCoords{1}(i,2)) = true;
    end

    % Remove the internal boundary from the fluid domain
    internalDomain = logical(int8(fluidDomain) - int8(domainBoundary));
    
    % Convert to array form
   	domainBoundary = logical(grid2array(domainBoundary));
    internalDomain = logical(grid2array(internalDomain));

end