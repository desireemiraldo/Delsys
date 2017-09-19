function Acc = ReadDelsys1(FilePath, ChannelType, Signal,numSensor)

% fileID = fopen(FilePath,'r');
% formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';
% Data = textscan(fileID, formatSpec, 'Delimiter',',',  'ReturnOnError', false);

Data =  readtable(FilePath);

Acc = []; t = [];
column = numSensor + 4;
for i = 1: size(Data,1)
   if  strcmp(Data.Var1(i),ChannelType) && strcmp(Data.Var2(i),Signal)
       Acc = [Acc; Data(i,column)];
       t = [t; str2num(cell2mat(Data.Var4(i)))];
   end
end
Acc = table2array(Acc);
% t = str2num(cell2mat(t.Var4));

Acc = [t,Acc];
end