function pressureGrad = getPressureGradient(param,pivData,initialParticles,particlePaths,currentTime)
    
    % 
    
    % Map material acceleration and velocity Laplacian vectors onto
    % original grid to determine pressure gradient
    % grad(P) = - Du/Dt + grad^2(u)
    
    % Initialise grid value arrays
    gridArrayLength = length(pivData.x);
    
    % Velocity distributions
    ux = pivData.ux{currentTime};
    uy = pivData.uy{currentTime};
    
    % For x component of pressure gradient
    DuxDt = zeros(1,gridArrayLength);
    ux_xx = DuxDt;
    ux_yy = DuxDt;
    
    % For y component
    DuyDt = DuxDt;
    uy_xx = DuxDt;
    uy_yy = DuxDt;
    
    % Get internal border of valid fluid domain
    [~,internalDomain] = getDomainBoundary(pivData);
    
    % x width of 2D grid
    Lx = length(unique(pivData.x));
    
    % Non-dimensional steps
    dx = 1;
    dy = 1;
    
    % Loop through original grid
    numNeighbours = param.kNN;
    for k = 1:gridArrayLength
        
        % Interpolate material acceleration from the virtual particles
        % Perform computations inside valid domain
        if pivData.domain(k)
        
            DuxDt(k) = kNNInterp(initialParticles.x,initialParticles.y,particlePaths.DuxDt,pivData.x(k),pivData.y(k),numNeighbours);
            DuyDt(k) = kNNInterp(initialParticles.x,initialParticles.y,particlePaths.DuyDt,pivData.x(k),pivData.y(k),numNeighbours);
        end
        
        % Calculate the velocity derivatives
        % Perform computations inside internal domain (valid points without
        % boundary)
        % Central difference scheme
        % 2D grid to array key
        % ij -> k | i+1,j -> k + 1 | i-1,j -> k - 1
        % ij -> k | i,j+1 -> k + Lx | i,j-1 -> k - Lx
        if internalDomain(k)
            
            % x component velocity derivatives
            ux_xx(k) = (ux(k+1) - 2*ux(k) + ux(k-1))/dx^2;
            ux_yy(k) = (ux(k+Lx) - 2*ux(k) + ux(k-Lx))/dy^2;
            
            % y component velocity derivatives
            uy_xx(k) = (uy(k+1) - 2*uy(k) + uy(k-1))/dx^2;
            uy_yy(k) = (uy(k+Lx) - 2*uy(k) + uy(k-Lx))/dy^2;
            
        end
        
    end
    
    % Calculate pressure gradient components (additional factor is due to
    % nondimensionalisation)
    pressureGrad.x = (1/param.Re0)*(ux_xx + ux_yy) - DuxDt;
    pressureGrad.y = (1/param.Re0)*(uy_xx + uy_yy) - DuyDt;
    
end