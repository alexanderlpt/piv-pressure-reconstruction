function array = grid2array(grid)

    % Measure grid
    [Ly,Lx] = size(grid);
    
    % Initialise array
    array = zeros(1,Lx*Ly);
    
    for i = 1:Ly
        array((i-1)*Lx+1:i*Lx) = grid(Ly+1-i,:);
    end
end