function outputValue = kNNInterp(xField,yField,fieldValues,xOutput,yOutput,numNeighbours)

    % Map field to output space using nearest neighbour interpolation
    % 
    
    % Distance between field points and output points
    dist = sqrt((xOutput - xField).^2 + (yOutput - yField).^2);
    
    % Choose numNeighbours closest points to take average
    [minDist,minIndx] = mink(dist,numNeighbours);
    
    % Inverse distance is used as a weight for the average
    % Small error is introduced to ensure that infinity is never calculated
    err = 1e-7;
    alpha = 1./(minDist + err);
    
    % Make weighted average
    closestFieldValues = fieldValues(minIndx);
    outputValue = sum(alpha.*closestFieldValues)./sum(alpha);
    
end