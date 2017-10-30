clc; close all; clear all;

Path = '.\Piloto\';

Folder = {'S001\','S002\','S003\','S004\','S005\','S006\','S007\','S008\'};

right = {[1,2,11],[4,3,11],[3,4,11],[3,4,11],[3,4,11],[3,4,11],[3,4,11],[3,4,11]};

left = {[1,2,11],[6,5,12],[5,6,12],[5,6,12],[5,6,12],[5,6,12],[5,6,12],[5,6,12]};

shank = {'right','left'};

side = randperm(2,1);

side = 1;

Sensor = eval(shank{side});


% --

Sub = 8;


% for Sub = 1: length(Folder)
Files = dir([Path,Folder{Sub},'*Rep*.xlsx']);
Files = {(Files(:).name)}';

for f = 9:length(Files)
    FilePath = [Path,Folder{Sub},Files{f}];
    
    tic
    Ftoe = ReadEMGWorks2(FilePath,Sensor{Sub}(3),'ACC Y');
    % Fheel = ReadEMGWorks2(FilePath,Sensor{Sub}(3),'ACC X');
    toc
    
    figure(side);
    plot(Ftoe(:,1),Ftoe(:,2),'b');
    
    
    % plot(FtoeR(:,1),FtoeR(:,2),'b', FtoeL(:,1),FtoeL(:,2),'r',...
    %    FheelR(:,1),FheelR(:,2),'b--', FheelL(:,1),FheelL(:,2),'r--')
    title([Files{f},' - ',shank{side}]);
    keyboard;
    
end
