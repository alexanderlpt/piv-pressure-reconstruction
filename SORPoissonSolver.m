function pressureField = SORPoissonSolver(param,pivData,sourceTerm,currentTime)

    %
    
    %% SETTINGS
    
    % Length of array representing 2D grid
    arrayLength = length(pivData.x);
    
    % x and y width of 2D grid
    Lx = length(unique(pivData.x));
    Ly = arrayLength/Lx;
    
    % Non-dimensional steps
    dx = 1; % = dy
    
    % SOR parameters
    maxIter = param.maxIter;                    % Maximum iterations
    minError = param.minError;                  % Convergence criteria
    if strcmp(param.omega,'opt')
        xi = (0.5*(cos(pi/Lx) + cos(pi/Ly)))^2;
        omegaOpt = 2*(1-sqrt(1-xi))/xi;             % Optimal SOR parameter for fastest convergence
        omega = omegaOpt;                           % Choice of SOR parameter
    else
        omega = param.omega;
    end
    
    % Initial pressure field
    pressureField = zeros(1,arrayLength);

    %% BOUNDARY CONDITIONS
    
    % Set Dirichlet boundary conditions (assuming farfield Bernoulli
    % pressure)
    % Obtain domain boundaries
    [domainBoundary,internalDomain] = getDomainBoundary(pivData);
    
    % Obtain velocity at boundaries
    uxBC = pivData.ux{currentTime}(domainBoundary);
    uyBC = pivData.uy{currentTime}(domainBoundary);
    
    % Calculate Bernoulli pressure at boundaries
    pressureBC = 0.5*(1 - (uxBC.^2 + uyBC.^2));
    pressureField(domainBoundary) = pressureBC;
    
    % Convert to 2D grid
    pressureField = array2grid(pivData,pressureField);
    internalDomain = array2grid(pivData,internalDomain);
    
    % Iterate algorithm towards convergence
    for iter = 1:maxIter
        
        % Record previous field before calculation
        oldField = pressureField;
    
        % Loop through grid
        for j = 1:Ly
            for i = 1:Lx

                % Bottom left starting point
                % Perform computations inside bounded domain
                if internalDomain(j,i)

                    % Calculate pressure change
                    dp = 0.25*(pressureField(j,i-1) + pressureField(j,i+1)...
                        + pressureField(j-1,i) + pressureField(j+1,i) - sourceTerm(j,i)*dx^2);

                    % Calculate new pressure field
                    pressureField(j,i) = (1-omega)*pressureField(j,i) + omega*dp;

                end
            end

        end
        
        % Check convergence (take minimum of difference and relative error)
        error = max(abs(pressureField - oldField)./abs(oldField),[],'all');
%         diff_error = max(abs(pressureField - oldField),[],'all');
        if error < minError %|| diff_error < minError
            break;
        end
    
    end
    
    if iter == maxIter
        disp('Solution has not converged');
    end
    
    % Return pressure coefficient
    pressureField = 2*pressureField;

end