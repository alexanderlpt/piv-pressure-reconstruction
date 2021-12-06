function potData = getPotentialFlow(param,xPoints,yPoints,flow)

    
    if strcmp(flow,'jet')
        % Get potential flow of impinging jet
        % Strength of jet
        Uinf = param.Uinf;
        C = param.Uinf/yPoints;

        % Jet impinges along centreline of domain
        x = -floor(xPoints/2):floor(xPoints/2)-1;
        y = 0:yPoints-1;

        % Create grid of points
        [X,Y] = meshgrid(x,y);

        % Convert to array format
        xArray = grid2array(X);
        yArray = grid2array(Y);

        % Calculate velocity field
        ux = C*xArray/Uinf;
        uy = -C*yArray/Uinf;

        % Calculate streamlines
        psi = -C.*xArray.*yArray;
        
        % Valid domain
        domain = ones(1,length(xArray));
        
    elseif strcmp(flow,'hemisphere')
        
        % Get potential flow of hemispere
        % Strength of doublet
        Uinf = param.Uinf;       % Free-stream velocity
        r = xPoints/8; % radius of hemisphere
        M = 2*pi*Uinf*r^2;

        % Hemisphere exists along centreline of domain
        x = -floor(xPoints/2):floor(xPoints/2)-1;
        y = 0:yPoints-1;

        % Create grid of points
        [X,Y] = meshgrid(x,y);

        % Calculate velocity field
        ux = Uinf -M*(X.^2 - Y.^2)./(2*pi*(X.^2 + Y.^2).^2);
        uy = -M*X.*Y./(pi*(X.^2 + Y.^2).^2);

        % Calculate streamlines
        psi = Uinf*Y -M*Y./(2*pi*(X.^2 + Y.^2));
        
        % Valid domain
        domain = ones(yPoints,xPoints);
        domain(X.^2 + Y.^2 <= r^2) = 0;
        
        % Convert to array format
        xArray = grid2array(X);
        yArray = grid2array(Y);
        ux = grid2array(ux)/Uinf;
        uy = grid2array(uy)/Uinf;
        psi = grid2array(psi);
        domain = grid2array(domain);
        
    end
    
    % Calculate pressure coefficient field
    Cp = 1 - (ux.^2 + uy.^2);
    
    % Domain
    potData.domain = domain;
    
    % Store data in struct format for solver
    potData.type = '2D';
    potData.x = xArray;
    potData.y = yArray;
    
    potData.ux = cell(1,1);
    potData.uy = potData.ux;
    potData.cp = potData.ux;
    potData.psi = potData.ux;
    
    potData.ux{1} = ux;
    potData.uy{1} = uy;
    
    potData.cp{1} = Cp;
    potData.psi{1} = psi;
    
    
end