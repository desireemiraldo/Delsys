clear all; clc; close all
% FsFP = 300; %Hz

% Path = 'C:\Users\desir\Desktop\Acelerometro - Aug17\Delsys\Esteira\Piloto\';
% File = 'RNW\RNW_Delcalco_Confortavel_Rep_1.17';
% ext = '.xls';
% 
% 
% FilePath = [Path,File,ext];

FilePath = '.\Piloto\RNW\RNW_Descalco_Rapido_Rep_4.26.xls';

x= ReadEMGWorks(FilePath,1,'ACC X',3);

%Fheel = ReadEMGWorks(FilePath,11,'ACC X',3);
Ftoe = ReadEMGWorks(FilePath,11,'ACC Y',3);




figure();
plot(Ftoe(:,1),Ftoe(:,2),'k')%, Fheel(:,1),Fheel(:,2),'b')
xlim([3 30]);
title(FilePath);


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
