% PRESSURE FIELD RECONSTRUCTION
% Author: Alexander Le Poer Trench
% Last Change: 06/12/2021

%% SETUP 
clear; clc; close all
plottingParameters();
set(groot,'defaultFigurePosition',[100 100 600 400]); % For large images

%% PARAMETERS

% Experimental Conditions
param.rho = 1.225;                                      % Ambient Air Density [kg/m^3]
param.visc = 1.81e-5;                                   % Ambient Air Viscosity [Pa.s]
param.Uinf = 1;                                         % Free-stream velocity [m/s] 
param.Pinf = 0;                                         % Free-stream static pressure[Pa]
param.L0 = 0.001;                                       % Length per pixel [m]
param.dt = 0.01;                                        % Time-step [s]

% Nondim parameters
param.P0 = param.rho*param.Uinf^2;                      % Scale pressure [Pa]
param.t0 = param.L0/param.Uinf;                         % Scale time [s]
param.Re0 = param.rho*param.Uinf*param.L0/param.visc;   % Scale Reynolds number

%% READ DATA

% Type of data to be read
datatype = 'piv'; % 'piv', 'cfd', 'pot_hemisphere', 'pot_jet'

if strcmp('piv',datatype) || strcmp('cfd',datatype)

    disp('Choose a velocity field data folder.');
    folderPath = uigetdir();
    
    % Type of reconstruction approach (Eulerian or Lagrangian)
    approach = 'lag'; % 'eul' or 'lag
    
    if strcmp('piv',datatype)
        % Read PIV data
        pivData = extractPIVData(param,folderPath);
    
    elseif strcmp('cfd',datatype)
        % Read CFD data
        pivData = extractCFDData(param,folderPath);

    end

elseif strcmp('pot',datatype(1:3))
    
    % Read potential flow data
    Lx = 200; % Horizontal length [pixels]
    Ly = 200; % Vertical length [pixels]
    pivData = getPotentialFlow(param,Lx,Ly,datatype(5:end));
    
    % Type of reconstruction approach
    approach = 'eul'; % only 'eul'
end

%% SETTINGS

% General settings
maxPoints = length(pivData.domain(pivData.domain==true)); 

param.kNN = 7;          % Nearest Neighbour (<= maxPoints)
param.maxIter = 1e4;    % SOR max iter (>=1)
param.minError = 1e-6;  % SOR min error
param.omega = 'opt';    % SOR over-relaxation factor (between 1 and 2, or 'opt')

% Lagrangian settings
numParticles = 500;     % Number of particles (<= maxPoints)
order = 1;              % Order of polynomial fit (>= 1)
desiredSteps = 5;       % Particle track length (>= order + 1)

% Time-step settings
startingTimestep = 1;   % Time-step in folder to start process (1 for potential flow)
numTimesteps = 5;       % Number of time-steps to have pressure reconstructed

% Where to save velocity and reconstructed pressure plots
disp('Choose a save folder.');
saveFolder = uigetdir();

%% PRESSURE RECONSTRUCTION

% Loop through timesteps in chosen folder
for currentTime = startingTimestep:numTimesteps
    
    if strcmp(approach,'eul')
        
        % Perform Eulerian pressure reconstruction
        pressureField = eulerianPressureSolver(param,pivData,currentTime);
        
    elseif strcmp(approach,'lag')
        
        % Perform Lagrangian pressure reconstruction
        [initialParticles,particlePaths] = createParticlePaths(param,pivData,numParticles,currentTime,desiredSteps,numTimesteps);
        
        particlePaths = getMaterialAcceleration(particlePaths,order);

        pressureGrad = getPressureGradient(param,pivData,initialParticles,particlePaths,currentTime);

        pressureField = getPressure(param,pivData,pressureGrad,currentTime);
    end

    %% PLOTTING
    
    % Coordinate matrices
    X = array2grid(pivData,pivData.x);
    Y = array2grid(pivData,pivData.y);
    
    % Total velocity matrix
    ux = pivData.ux{currentTime};
    uy = pivData.uy{currentTime};
    u = sqrt(ux.^2 + uy.^2);
    U = array2grid(pivData,u);
    
    % Remove all non valid points
    gridDomain = array2grid(pivData,pivData.domain);
    U(~gridDomain) = NaN;
    pressureField(~gridDomain) = NaN;
    
    % Plot total velocity contour
    figure(1);
    contourf(X,Y,U,30,'LineColor', 'none');
    colorbar('TickLabelInterpreter','latex');
    xlabel('$x$'); ylabel('$y$');
    legend('$U/U_{\infty}$');
    axis equal;

    % Save plot
    plotName = ['piv_total_vel_contour',num2str(currentTime)];
    savePDF(saveFolder,plotName);
    
    % Plot pressure coefficient contour
    figure(2);
    contourf(X,Y,pressureField,30,'LineColor', 'none');
    colorbar('TickLabelInterpreter','latex');
    xlabel('$x$'); ylabel('$y$');
    legend('$C_{p}$');
    axis equal;

    % Save plot
    plotName = ['piv_cp_contour',num2str(currentTime)];
    savePDF(saveFolder,plotName);
   
    % If CFD or potential flow data plot pressure coefficient error field
    if strcmp('cfd',datatype) || strcmp('pot',datatype(1:3))
        
        % Pressure coefficient from CFD
        Cp = pivData.cp{currentTime};
        Cp = array2grid(pivData,Cp);
        Cp(~gridDomain) = NaN;
        
        % Pressure coefficient error
        error = abs(pressureField-Cp)./abs(Cp);
        
        % Change to logarithmic basis
        error = log10(error);
        
        % Plot pressure coefficient error contour
        figure(3);
        contourf(X,Y,error,30,'LineColor', 'none');
        colorbar('TickLabelInterpreter','latex');
        xlabel('$x$'); ylabel('$y$');
        legend('log$_{10}(\epsilon_{C_{p}})$');
        axis equal;

        % Save plot
        plotName = ['piv_cp_error_contour',num2str(currentTime)];
        savePDF(saveFolder,plotName);
        
    end
    
    % End loop if potential flow data is chosen
    if strcmp('pot',datatype(1:3))
        break;
    end

end
