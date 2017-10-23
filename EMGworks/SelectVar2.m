function Var = SelectVar2(Data,Sensor,Channel)

if Sensor==11 || Sensor==12
    label = strcat("Trigno FSR Adapter ",num2str(Sensor),":"," ",Channel," ", num2str(Sensor));
else
    label = strcat("Trigno IM sensor ",num2str(Sensor),":"," ",Channel," ", num2str(Sensor)," (IM)");
end


for i = 1: length(Data.colheaders)
    if strcmp(string(Data.colheaders(i)),label)
            Var = Data.data(:,i-1:i,:);
    end
end

end




