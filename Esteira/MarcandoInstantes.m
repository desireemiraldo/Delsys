%clear all; clc; close all
FsFP = 300; %Hz

Path = 'C:\Users\BMClab\Downloads\Desiree\Drop Foot Project\Identificacao de evento\Processamento\Delsys\Esteira\Piloto\';
File = 'RNW\RNW_Calçado_Confortavel_Rep_2.2';
ext = '.xls';


FilePath = [Path,File,ext];



Fheel = ReadEMGWorks(FilePath,11,'ACC X',3);
Ftoe = ReadEMGWorks(FilePath,11,'ACC Y',3);




figure(2);
plot(Ftoe(:,1),Ftoe(:,2),'k', Fheel(:,1),Fheel(:,2),'b')
xlim([3 25])


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