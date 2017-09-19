function Acc = ReadDelsys2(FilePath, ChannelType, Signal)

% fileID = fopen(FilePath,'r');
% formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';
% Data = textscan(fileID, formatSpec, 'Delimiter',',',  'ReturnOnError', false);

Data =  readtable(FilePath);

Acc = []; t = [];

for i = 1: height(Data)
   if  strcmp(Data.Var1(i),ChannelType) && strcmp(Data.Var2(i),Signal)
       Acc = [Acc; Data(i,5:end-1)];
       t = [t; str2num(cell2mat(Data.Var4(i)))];
   end
end
Acc = table2array(Acc);
% t = str2num(cell2mat(t.Var4));

Acc = [t,Acc];
end