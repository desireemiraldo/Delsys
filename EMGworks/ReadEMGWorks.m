function Data = ReadEMGWorks(FilePath,Sensor,Channel,deleteTime)

FileData = importdata(FilePath);

if Sensor==11 || Sensor==12
    label = strcat("Trigno FSR Adapter ",num2str(Sensor),":"," ",Channel," ", num2str(Sensor));
else
    label = strcat("Trigno IM sensor ",num2str(Sensor),":"," ",Channel," ", num2str(Sensor)," (IM)");
end

if strcmp(Channel,'EMG')
    Fs = 1/(FileData.data(2,1)); 
elseif strcmp(Channel,'Mag X') || strcmp(Channel,'Mag Y') || strcmp(Channel,'Mag Z') 
    Fs = 1/(FileData.data(2,15)); 
else
    Fs = 1/(FileData.data(2,3)); 
end


% delete null data
while FileData.data(end,2)==0
    FileData.data(end,:) = [];
end
t = FileData.data(end,1);

% Delete first and last 'deleteTime' seconds
first = ceil(deleteTime*Fs)+1;
last = ceil((t-deleteTime)*Fs)+1;


for i = 1: length(FileData.colheaders)
    if strcmp(string(FileData.colheaders(i)),label)
        Data = FileData.data(first:last,i-1:i);
    end
end

end





