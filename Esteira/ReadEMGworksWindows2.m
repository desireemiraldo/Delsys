function Var = ReadEMGworksWindows2(Data,Sensor,Signal, deltaT, deleteTime)

FileData = importdata(FilePath);

if Sensor==11 || Sensor==12
    label = strcat("Trigno FSR Adapter ",num2str(Sensor),":"," ",Channel," ", num2str(Sensor));
else
    label = strcat("Trigno IM sensor ",num2str(Sensor),":"," ",Channel," ", num2str(Sensor)," (IM)");
end

if strcmp(Channel,'EMG')
    Fs = 1/(FileData.data(2,1)); 
    linesNumber = ceil(deltaT*Fs);
elseif strcmp(Channel,'Mag X') || strcmp(Channel,'Mag Y') || strcmp(Channel,'Mag Z') 
    Fs = 1/(FileData.data(2,15)); 
    linesNumber = ceil(deltaT*Fs);
else
    Fs = 1/(FileData.data(2,3)); 
    linesNumber = ceil(deltaT*Fs);
end

% delete null data
while FileData.data(end,2)==0
    FileData.data(end,:) = [];
end
% Delete first and last 'deleteTime' seconds
del = floor(Fs*deleteTime);
FileData.data(1:del,:) = []; FileData.data(end-del:end,:) = [];


t = FileData.data(end,1);
Win = floor(t/deltaT);

Var = zeros(linesNumber,2,Win);



for i = 1: length(FileData.colheaders)
    disp(FileData.colheaders(i));
    disp(label);
    if strcmp(string(FileData.colheaders(i)),label)
        for w = 1: Win
            Var(:,:,w) = FileData.data((w-1)*linesNumber+1:linesNumber*w,i-1:i);
        end
    end
end

end




