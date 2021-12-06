% AERO4560: Flight Mechanics 2 Assignment 1
% SID: 470427349

% Set desired plotting parameters

function plottingParameters()

    % Set Plotting parameters
    myred           = [216 30 49]/255;
    myblue          = [27 99 157]/255;
    myblack         = [0 0 0]/255;
    mygreen         = [0 128 0]/255;
    mycyan          = [2 169 226]/255;
    myyellow        = [251 194 13]/255;
    mygray          = [89 89 89]/255;
    set(groot,'defaultAxesColorOrder',[myblack;myblue;myred;mygreen;myyellow;mycyan;mygray]);
    alw             = 1;                        % AxesLineWidth
    fsz             = 13;                       % Fontsize
    lw              = 1.5;                        % LineWidth
    msz             = 7;                       % MarkerSize
    set(0,'defaultLineLineWidth',lw);           % set the default line width to lw
    set(0,'defaultLineMarkerSize',msz);         % set the default line marker size to msz
    set(0,'defaultAxesLineWidth',alw);           % set the default line width to lw
    set(0,'defaultAxesFontSize',fsz);         % set the default line marker size to msz
    set(0,'defaultFigureColor','w');
    set(0,'defaultAxesColor','w');

    % Set default interpreter to latex
    set(groot,'defaultAxesTickLabelInterpreter','latex'); 
    set(groot,'defaulttextinterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');

%     set(gca,'GridLineStyle','-');
%     set(gca,'MinorGridLineStyle','-');
%     set(gca,'GridColor','k');
%     set(gca,'MinorGridColor','k');
%     colorbar('TickLabelInterpreter','latex');

    % Figure size and position
%     set(groot,'defaultFigurePosition',[100 100 400 300]); % For small images
    set(groot,'defaultFigurePosition',[100 100 600 300]); % For large images
    
end