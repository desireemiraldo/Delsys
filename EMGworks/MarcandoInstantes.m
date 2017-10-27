clc; close all; clear all;

Path = '.\Piloto\';

Folder = {'S001\','S002\','S003\','S004\','S005\'};

right = {[1,2,11],[4,3,11],[3,4,11],[3,4,11],[3,4,11]};

left = {[1,2,11],[6,5,12],[5,6,12],[5,6,12],[5,6,12]};

Shank = {'right','left'};


% --

Sub = 4;


% for Sub = 1: length(Folder)
Files = dir([Path,Folder{Sub},'*Rep*.xlsx']);
Files = {(Files(:).name)}';

for f = 3:length(Files)
    FilePath = [Path,Folder{Sub},Files{f}];
    
    for s =1:2
        if s ==1
        Sensor = eval(Shank{s});
        
        tic
        Ftoe = ReadEMGWorks2(FilePath,Sensor{Sub}(3),'ACC Y');
        %FheelL = ReadEMGWorks2(FilePath,Sensor{Sub}(3),'ACC X');
        toc
        
        figure(s);
        plot(Ftoe(:,1),Ftoe(:,2));
        
        %     plot(FtoeR(:,1),FtoeR(:,2),'b', FtoeL(:,1),FtoeL(:,2),'r',...
        %         FheelR(:,1),FheelR(:,2),'b--', FheelL(:,1),FheelL(:,2),'r--')
        title([Files{f},' - ',Shank{s}]);
        keyboard;
        end
    end

    
end
