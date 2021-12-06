function grid = array2grid(pivData,array)
    
    % Measure lengths of grid
    arrayLength = length(pivData.x);
    Lx = length(unique(pivData.x));
    Ly = arrayLength/Lx;
    
    % Initialise grid
    grid = zeros(Ly,Lx);
    
    % Arrange array elements into grid
    for i = 1:Ly 
        grid(Ly+1-i,:) = array((i-1)*Lx+1:i*Lx);
    end

end