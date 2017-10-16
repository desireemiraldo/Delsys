%function newData = DeleteNullData(Path,Files, deleteTime, deltaT)

deleteTime = 3; deltaT = 2.75;
Path = '.\Piloto\RNW\';

Files = {'RNW_Calçado_Confortavel_Rep_1.1',...
    'RNW_Calçado_Confortavel_Rep_2.2',...
    'RNW_Calçado_Confortavel_Rep_3.3',...
    'RNW_Calçado_Confortavel_Rep_4.4',...
    'RNW_Calçado_Confortavel_Rep_5.5',...
    'RNW_Delcalco_Confortavel_Rep_1.17',...
    'RNW_Delcalco_Confortavel_Rep_2.18',...
    'RNW_Delcalco_Confortavel_Rep_4.20',...
    'RNW_Delcalco_Confortavel_Rep_5.21',...
    'RNW_Delcalco_Confortavel_Rep_6.22'};

Trials = 0;

for i = 1: length(Files)
    
    FileData = importdata([Path,Files{i},'.xls']);
    
    
    %% delete 'EMG' and 'MAG'
    j = 2;
    while j <= length(FileData.colheaders)
       if contains(string(FileData.colheaders(j)),'EMG') ||...
                contains(string(FileData.colheaders(j)),'Mag')
            % disp(FileData.colheaders(j))
            FileData.data(:,j-1:j) = [];
            FileData.colheaders(:,j-1:j) = [];
            j = j - 2;
       end
        j = j + 2;
    end
    
    %% delete null data
    while FileData.data(end,2)==0
        FileData.data(end,:) = [];
    end
    
    % if strcmp(Channel,'EMG')
    %     Fs = 1/(FileData.data(2,1));
    %     linesNumber = ceil(deltaT*Fs);
    % elseif strcmp(Channel,'Mag X') || strcmp(Channel,'Mag Y') || strcmp(Channel,'Mag Z')
    %     Fs = 1/(FileData.data(2,15));
    %     linesNumber = ceil(deltaT*Fs);
    % else
    Fs = 1/(FileData.data(2,3));
    linesNumber = ceil(deltaT*Fs);
    % end
    
    
    % Delete first and last 'deleteTime' seconds
    del = floor(Fs*deleteTime);
    FileData.data(1:del,:) = []; FileData.data(end-del:end,:) = [];
    
    
    t = (FileData.data(end,1)-FileData.data(1,1));
    Win = floor(t/deltaT);
    
    newData = zeros(linesNumber,size(FileData.data,2),Win);
    for w = 1: Win
        newData(:,:,w) = FileData.data((w-1)*linesNumber+1:linesNumber*w,:);
       
    end
    
    
    Data = zeros(linesNumber,size(FileData.data,2),Trials + Win);
    if Trials ~=0
        Data(:,:,1:Trials) = Temp;
        Data(:,:,1:Trials+1:end) = newData;
    else
        Data = newData;
    end
    Temp = Data;
    Trials = Trials + Win;

end

%end
