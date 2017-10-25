function Data = ReadEMGWorks2(FilePath,Sensor,Channel)

FileData = importdata(FilePath);

if Sensor==11 || Sensor==12
    label = strcat("Trigno FSR Adapter ",num2str(Sensor),":"," ",Channel," ", num2str(Sensor));
else
    label = strcat("Trigno IM sensor ",num2str(Sensor),":"," ",Channel," ", num2str(Sensor)," (IM)");
end


for i = 1: length(FileData.colheaders)
    if strcmp(string(FileData.colheaders(i)),label)
        Data = FileData.data(:,i-1:i);
    end
end

end





