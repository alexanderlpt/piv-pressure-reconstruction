function pivData = extractPIVData(param,folderPath)

    % extractPIVData
    % Read all .dat files associated with PIV data and
    % commit to memory as struct
    %
    % Author: Alexander Le Poer Trench
    % 
    % INPUTS:
    % folderPath        = path to folder where PIV data files are stored (contains
    %                       M files)
    %
    % OUTPUTS:
    % pivData           = (struct) contains all formatted PIV data from provided
    %                       files
    % pivData.x         = (int array) x coordinate array for velocity measurements
    % pivData.y         = (int array) y coordinate array for velocity measurements
    % pivData.ux        = (1 x M cell) x component velocity field array
    % pivData.uy        = (1 x M cell) y component velocity field array
    % pivData.domain    = (logical array) describes valid points inside flow domain

    % How many files in folder?
    folderObj = dir([folderPath '/*.dat']);
    timesteps = size(folderObj,1);

    % Read PIV data files for all timesteps
    % Create a struct for the data
    % Store cells with each changing data variable: ux,uy
    pivData.ux = cell(1,timesteps);
    pivData.uy = pivData.ux;
    
    % First file
    filePath = [folderPath,'\',folderObj(1).name];
    fileID = fopen(filePath);
    
    % Check if data is 2D or 3D
    line = textscan(fileID,'%s',1,'Delimiter','\t', 'Headerlines',3);
    fclose(fileID);
    numCols = length(cell2mat(strfind(line{1},' '))) + 1;
    
    % Read data
    % numCols = 5 -> 2D PIV data (x,y,Ux,Uy,domain)
    % numCols = 6 -> 3D PIV data (x,y,Ux,Uy,Uz,domain)
    fileID = fopen(filePath);
    dataString = repmat('%f ',1,numCols);
    dataString(end) = [];
    data = textscan(fileID,dataString,'HeaderLines',3,'Delimiter','\t');
    fclose(fileID);
    
    % Label pivData with type of data
    if numCols == 5
        pivData.type = '2D';
    else
        pivData.type = '3D';
    end
    
    % X,Y,domain data does not vary for same PIV case
    pivData.x = data{1}';
    pivData.y = data{2}';
    pivData.domain = logical(data{numCols}');
    
    % Fill up struct
    for step = 1:timesteps
        
        % Access file
        filePath = [folderPath,'\',folderObj(step).name];
        fileID = fopen(filePath);
        
        % Read data
        data = textscan(fileID,dataString,'HeaderLines',3,'Delimiter','\t');
        fclose(fileID);
        
        % Store velocity arrays in struct
        pivData.ux{step} = data{3}'/param.Uinf;
        pivData.uy{step} = data{4}'/param.Uinf;

    end
    
end
