function FilesData = BuildingTrials2(Path,Files, deleteTime, deltaT)

Trials = 0;

for i = 1: length(Files)
    
    FilesData = importdata([Path,Files{i},'.xls']);
    % disp(Files{i})
    
    
    %% delete 'EMG' and 'MAG'
    j = 2;
    while j <= length(FilesData.colheaders)
       if contains(string(FilesData.colheaders(j)),'EMG') ||...
                contains(string(FilesData.colheaders(j)),'Mag')
            % disp(FileData.colheaders(j))
            FilesData.data(:,j-1:j) = [];
            FilesData.colheaders(:,j-1:j) = [];
            j = j - 2;
       end
        j = j + 2;
    end
    
    %% delete null data
    while FilesData.data(end,2)==0
        FilesData.data(end,:) = [];
    end
    
    % if strcmp(Channel,'EMG')
    %     Fs = 1/(FileData.data(2,1));
    %     linesNumber = ceil(deltaT*Fs);
    % elseif strcmp(Channel,'Mag X') || strcmp(Channel,'Mag Y') || strcmp(Channel,'Mag Z')
    %     Fs = 1/(FileData.data(2,15));
    %     linesNumber = ceil(deltaT*Fs);
    % else
    Fs = 1/(FilesData.data(2,3));
    linesNumber = ceil(deltaT*Fs);
    % end
    
    
    % Delete first and last 'deleteTime' seconds
    del = floor(Fs*deleteTime);
    FilesData.data(1:del,:) = []; FilesData.data(end-del:end,:) = [];
    
    
    t = (FilesData.data(end,1)-FilesData.data(1,1));
    Win = floor(t/deltaT);
    
    newData = zeros(linesNumber,size(FilesData.data,2),Win);
    for w = 1: Win
        newData(:,:,w) = FilesData.data((w-1)*linesNumber+1:linesNumber*w,:);
       
    end
    
    
    Data = zeros(linesNumber,size(FilesData.data,2),Trials + Win);
    if Trials ~=0
        Data(:,:,1:Trials) = Temp;
        Data(:,:,Trials+1:end) = newData;
    else
        Data = newData;
    end
    Temp = Data;
    Trials = Trials + Win;
    clear newData; 

end
FilesData.data = Data;
clear Temp; clear Data;
end
