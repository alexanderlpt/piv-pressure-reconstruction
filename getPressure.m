function pressureField = getPressure(param,pivData,pressureGrad,currentTime)

    % 
    % Calculate lagrangian pressure by Poisson solver
    
    %% CALCULATE PRESSURE LAPLACIAN FOR POISSON SOLVER
    
    % Length of array representing 2D grid
    arrayLength = length(pivData.x);
    
    % x and y width of 2D grid
    Lx = length(unique(pivData.x));
    
    % Non-dimensional steps
    dx = 1;
    dy = dx;
    
    % Obtain domain boundaries
    [~,internalDomain] = getDomainBoundary(pivData);
    
    % Calculate pressure laplacian
    grad2p = zeros(1,arrayLength);
    
    % Loop through grid
    for k = 1:arrayLength

        % Bottom left starting point
        % Perform computations inside bounded domain
        if internalDomain(k)

            % Compute derivatives of pressure gradient (central difference scheme)
            d2p_d2x = (pressureGrad.x(k+1)-pressureGrad.x(k-1))/(2*dx);
            d2p_d2y = (pressureGrad.y(k+Lx)-pressureGrad.y(k-Lx))/(2*dy);

            % Compute pressure laplacian
            grad2p(k) = d2p_d2x + d2p_d2y;

        end

    end
    
    % Convert laplacian to grid format
    grad2p = array2grid(pivData,grad2p);
    
    %% POISSON SOLVER
    pressureField = SORPoissonSolver(param,pivData,grad2p,currentTime);
    
end