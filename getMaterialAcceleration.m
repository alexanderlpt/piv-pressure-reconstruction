function particlePaths = getMaterialAcceleration(particlePaths,order)

    % 
    
    % Extract particlePaths features
    time = particlePaths.time;
    ux = particlePaths.ux;
    uy = particlePaths.uy;
    
    % Initialise material acceleration arrays 
    DuxDt = zeros(1,length(ux));
    DuyDt = DuxDt;
    
    % Loop over each particle
    for p = 1:length(ux)
        
        % Fit a polynomial of time vs ux or uy to each particle trajectory
        % Linear slope (2nd polyCoeff) represents the material acceleration
        % at t = 0 (current time)
        
        % x component material acceleration
        [~, UXCoeffs, ~] = leastSquaresPolynomial(time,ux{p},order);
        DuxDt(p) = UXCoeffs(2);
        
        % y component material acceleration
        [~, UYCoeffs, ~] = leastSquaresPolynomial(time,uy{p},order);
        DuyDt(p) = UYCoeffs(2);
        
    end
    
    % Store material acceleration values in particlePaths structure
    particlePaths.DuxDt = DuxDt;
    particlePaths.DuyDt = DuyDt;

end