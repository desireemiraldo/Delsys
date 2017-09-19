%clear all; clc; close all
FsFP = 300; %Hz

Path = 'C:\Users\BMClab\Downloads\Desiree\Acelerometro\Acelerometro-GitHub\Esteira\Dropfoot\';
File = '\MD\Descalco3_Dual_170831_1';
FilePath = [Path,File];
csv = '-Delsys 1.csv';
numSensors = [11,12];

% Cortex data
    Forces = importdata([FilePath,'.forces']);
    Fy(:,1) = (Forces.data(:,1))/FsFP;
    Fy(:,2) = Forces.data(:,10)/2;
    
%     Contact1 = ReadDelsys1([FilePath,csv], 'EMG', 'EMG',numSensors);
%     Contact2 = ReadDelsys1([FilePath,csv], 'ACC', 'ACC Pitch',numSensors);
%     Contact3 = ReadDelsys1([FilePath,csv], 'ACC', 'ACC Roll',numSensors);
%     Contact4 = ReadDelsys1([FilePath,csv], 'ACC', 'ACC Yaw',numSensors);
   
    
    figure(2);
    plot(Fy(:,1),Fy(:,2),'k')
    
%     hold on
%     
%     plot(Contact1(:,1),Contact1(:,2),'b')
%     plot(Contact1(:,1),Contact1(:,3),'r')
%     
%     plot(Contact2(:,1),Contact2(:,2),'b')
%     plot(Contact2(:,1),Contact2(:,3),'r')
%  
%     plot(Contact3(:,1),Contact3(:,2),'b')
%     plot(Contact3(:,1),Contact3(:,3),'r')
%     
%     plot(Contact4(:,1),Contact4(:,2),'b')
%     plot(Contact4(:,1),Contact4(:,3),'r')