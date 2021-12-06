function [initialParticles,particlePaths] = createParticlePaths(param,pivData,numParticles,currentTime,desiredSteps,numTimesteps)

    % createParticlePaths
    % Particles start at specific timestep, and are tracked using the 
    % first order Euler integration method for their evolving coordinates
    % 
    % Author: Alexander Le Poer Trench
    % 
    % INPUTS:
    % folderPath        = path to folder where PIV data files are stored (contains
    %                       M files)
    %
    % OUTPUTS:
    % particlePaths           = (struct) contains all formatted PIV data from provided
    %                       files
    % pivData.X         = (1 x M cell) x coordinate matrix for each time step
    % pivData.Y         = y coordinate matrix
    % pivData.UX        = x component velocity field matrix 
    % pivData.UY        = y component velocity field matrix 
    % pivData.domain    = (logical matrix) matrix describing flow domain
    
    % Create particle cell structures
    % Provide particles with initial values
    particlePaths = randomPopulation(pivData,numParticles,currentTime);
    initialParticles.x = cell2mat(particlePaths.x);
    initialParticles.y = cell2mat(particlePaths.y);

    % Determine timesteps that must be accessed relative to currentTime
    stepBack = floor((desiredSteps-1)/2);
    stepForward = ceil((desiredSteps-1)/2);
    distForward = numTimesteps-currentTime;
    distBack = currentTime-1;
    
    if stepBack <= distBack && stepForward <= distForward
        relativeSteps = -stepBack:stepForward;
    elseif stepBack <= distBack && stepForward > distForward
        relativeSteps = -(stepBack+(stepForward-distForward)):distForward;
    else
        relativeSteps = -distBack:(stepForward + (stepBack - distBack));
    end
    
    timesteps = relativeSteps + currentTime;
    currentIndx = find(timesteps == currentTime);
    
    % Create time array inside particle struct
    % Nondim time-step
    dt = param.dt/param.t0;
    particlePaths.time = relativeSteps*dt;
    
    % Loop through timesteps centered around currentTime for each particle
    % Perform Euler integration backwards and forwards in time with dt
    numNeighbours = param.kNN;
    for p = 1:numParticles
        
        % currentTime values
        x0 = particlePaths.x{p};
        y0 = particlePaths.y{p};
        ux0 = particlePaths.ux{p};
        uy0 = particlePaths.uy{p};
        
        % Preallocated temporary arrays
        particleX = zeros(1,desiredSteps);
        particleY = particleX;
        particleUx = particleX;
        particleUy = particleX;
        
        particleX(currentIndx) = x0;
        particleY(currentIndx) = y0;
        particleUx(currentIndx) = ux0;
        particleUy(currentIndx) = uy0;
        
        % Forwards in time
        for t = 1:length(timesteps(timesteps > currentTime))
            
            % Position of next data element in preallocated array
            tForward = t + currentIndx;
            
            % Next particle position
            nextX = particleX(tForward-1) + particleUx(tForward-1)*dt;
            nextY = particleY(tForward-1) + particleUy(tForward-1)*dt;
            
            % Next particle position velocities
            nextFieldUx = pivData.ux{timesteps(tForward)};
            nextFieldUy = pivData.uy{timesteps(tForward)};
            
            % Interpolate velocity from next field at next particle position
            nextUx = kNNInterp(pivData.x,pivData.y,nextFieldUx,nextX,nextY,numNeighbours);
            nextUy = kNNInterp(pivData.x,pivData.y,nextFieldUy,nextX,nextY,numNeighbours);
            
            % Store new particle properties
            particleX(tForward) = nextX;
            particleY(tForward) = nextY;
            particleUx(tForward) = nextUx;
            particleUy(tForward) = nextUy;
            
            % Need at least two points in particle trajectory
            % Now check if particle has left domain, if so break loop
            if ~inDomain(pivData,nextX,nextY)
                break;
            end
            
        end
        
        % Backwards in time
        for t = 1:length(timesteps(timesteps < currentTime))
            
            % Position of next data element in preallocated array
            % (backwards)
            tBack = currentIndx - t;
            
            % Next particle position (backwards in time)
            nextX = particleX(tBack+1) - particleUx(tBack+1)*dt;
            nextY = particleY(tBack+1) - particleUy(tBack+1)*dt;
            
            % Next particle position velocities
            nextFieldUx = pivData.ux{timesteps(tBack)};
            nextFieldUy = pivData.uy{timesteps(tBack)};
            
            % Interpolate velocity from next field at next particle position
            numNeighbours = 4;
            nextUx = kNNInterp(pivData.x,pivData.y,nextFieldUx,nextX,nextY,numNeighbours);
            nextUy = kNNInterp(pivData.x,pivData.y,nextFieldUy,nextX,nextY,numNeighbours);
            
            % Store new particle properties
            particleX(tBack) = nextX;
            particleY(tBack) = nextY;
            particleUx(tBack) = nextUx;
            particleUy(tBack) = nextUy;
            
            % Need at least two points in particle trajectory
            % Now check if particle has left domain, if so break loop
            if ~inDomain(pivData,nextX,nextY)
                break;
            end
            
        end
        
        % Store filled arrays in particle cell structures
        particlePaths.x{p} = particleX;
        particlePaths.y{p} = particleY;
        particlePaths.ux{p} = particleUx;
        particlePaths.uy{p} = particleUy;
        
    end
    
    % Cut off particle paths when too long (distance metric)
end

% Populate domain with particles randomly
function particleDist = randomPopulation(pivData,numParticles,currentTime)

    % Indices of points inside domain
    domainIndx = find(pivData.domain);
    
    % Choose numParticles number of positions randomly inside the domain
    randIndx = randsample(domainIndx,numParticles);
    
    % Store initial randomly chosen x and y coordinates and velocities for particles
    % Convert arrays to cells
    particleDist.x = num2cell(pivData.x(randIndx));
    particleDist.y = num2cell(pivData.y(randIndx));
    
    currentUx = pivData.ux{currentTime};
    currentUy = pivData.uy{currentTime};
    
    particleDist.ux = num2cell(currentUx(randIndx));
    particleDist.uy = num2cell(currentUy(randIndx));

end

function in = inDomain(pivData,targetX,targetY)

    % Check if point is out of bounds of domain
    
    % Bounds
    xMax = max(pivData.x);
    xMin = min(pivData.x);
    yMax = max(pivData.y);
    yMin = min(pivData.y);
    
    if targetX < xMin || targetX > xMax || targetY < yMin || targetY > yMax
        in = false;
    else
        in = true;
    end

end