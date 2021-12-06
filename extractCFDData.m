function cfdData = extractCFDData(param,folderPath)

    % 
    % How many files in folder?
    folderObj = dir([folderPath '/*.dat']);
    timesteps = size(folderObj,1);

    % Read data files for all timesteps
    % Create a struct for the data
    % Store cells with each changing data variable: ux,uy
    cfdData.type = '3D';
    cfdData.ux = cell(1,timesteps);
    cfdData.uy = cfdData.ux;
    cfdData.ps = cfdData.ux;
    cfdData.pt = cfdData.ux;
    
    dataString = repmat('%f ',1,8);
    dataString(end) = [];
    
    % Fill up struct
    for step = 1:timesteps
        
        % Access file
        filePath = [folderPath,'\',folderObj(step).name];
        fileID = fopen(filePath);
        
        % Read data
        data = textscan(fileID,dataString,'HeaderLines',3,'Delimiter','\t');
        fclose(fileID);
        
        % X,Y,domain data does not vary for same PIV case
        if step == 1
            cfdData.x = data{1}';
            cfdData.y = data{2}';
            cfdData.domain = logical(data{8}');
        end
        
        % Store velocity arrays in struct
        cfdData.ux{step} = data{3}'/param.Uinf;
        cfdData.uy{step} = data{4}'/param.Uinf;
        cfdData.cp{step} = data{6}'/param.P0;
        cfdData.pt{step} = data{7}'/param.P0;

    end
    

end