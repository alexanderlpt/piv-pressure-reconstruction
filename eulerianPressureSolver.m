function pressureField = eulerianPressureSolver(param,pivData,currentTime)

    % Calculate Eulerian pressure by Poisson solver
    
    %% CALCULATE VELOCITY SOURCE TERM FOR POISSON SOLVER
    
    % Length of array representing 2D grid
    arrayLength = length(pivData.x);
    
    % x and y width of 2D grid
    Lx = length(unique(pivData.x));
    
    % Non-dimensional steps
    dx = 1;
    dy = dx;
    
    % Obtain domain boundaries
    [~,internalDomain] = getDomainBoundary(pivData);
    
    % Calculate velocity source terms
    ux = pivData.ux{currentTime};
    uy = pivData.uy{currentTime};
    
    % First derivatives
    dux_dx = zeros(1,arrayLength);
    dux_dy = dux_dx;
    duy_dx = dux_dx;
    duy_dy = dux_dx;
    
    % Compute first spatial velocity derivatives
    for k = 1:arrayLength
            
        % Check if point is inside fluid domain
        % If so perform calculation
        if internalDomain(k)

            % Compute first derivatives (central difference scheme)
            dux_dx(k) = (ux(k+1)-ux(k-1))/(2*dx);
            dux_dy(k) = (ux(k+Lx)-ux(k-Lx))/(2*dy);
            duy_dx(k) = (uy(k+1)-uy(k-1))/(2*dx);
            duy_dy(k) = (uy(k+Lx)-uy(k-Lx))/(2*dy);

        end
            
    end
    
    % Check if data is 2D or 3D
    % 2D implies divxy = 0
    if strcmp(pivData.type,'2D')
        
        % Compute source term
        velSource = -(dux_dx.^2 + 2*dux_dy.*duy_dy + duy_dy.^2);
        
    else
        % 3D data
        % divxy is non zero
        
        % Higher order terms
        d_divxy_dx = zeros(1,arrayLength);
        d_divxy_dy = d_divxy_dx;
        d2_divxy_d2x = d_divxy_dx;
        d2_divxy_d2y = d_divxy_dx;

        % Look behind or forward one step to obtain time derivative of
        % divergence
        timesteps = length(pivData.ux);
        if currentTime < timesteps
            j = 1;
        else
            j = -1;
        end

        % Compute divergence for next or previous timestep
        ux_next = pivData.ux{currentTime + j};
        uy_next = pivData.uy{currentTime + j};
        dux_dx_next = dux_dx;
        duy_dy_next = dux_dx_next;

        % Divergence terms
        divxy = dux_dx + duy_dy;

        % Spatial divergence derivatives
        for k = 1:arrayLength

            % Check if point is inside fluid domain
            % If so perform calculation
            if internalDomain(k)

                % Compute divergence terms for next or previous step
                dux_dx_next(k) = (ux_next(k+1)-ux_next(k-1))/(2*dx);
                duy_dy_next(k) = (uy_next(k+Lx)-uy_next(k-Lx))/(2*dy);

                % Compute first derivatives of divergence (central difference scheme)
                d_divxy_dx(k) = (divxy(k+1)-divxy(k-1))/(2*dx);
                d_divxy_dy(k) = (divxy(k+Lx)-divxy(k-Lx))/(2*dy);

                % Second derivatives
                d2_divxy_d2x(k) = (divxy(k+1)-2*divxy(k)+divxy(k-1))/(dx^2);
                d2_divxy_d2y(k) = (divxy(k+Lx)-2*divxy(k)+divxy(k-Lx))/(dy^2);

            end

        end

        % Time derivative
        divxy_next = dux_dx_next + duy_dy_next;
        dt = param.dt/param.t0;
        d_divxy_dt = j*(divxy_next - divxy)/dt;

        % Compute source term
        velSource = -(dux_dx.^2 + 2*dux_dy.*duy_dy + duy_dy.^2)...
            - (ux.*d_divxy_dx + uy.*d_divxy_dy)...
            + (param.L0/param.Re0)*(d2_divxy_d2x + d2_divxy_d2y)...
            - d_divxy_dt;
    end
     
    % Convert source term to grid format
    velSource = array2grid(pivData,velSource);
    
    %% POISSON SOLVER
    pressureField = SORPoissonSolver(param,pivData,velSource,currentTime);
    
end