Name = 'DELAY_Plot_and_Store_Rep_1.0';
FilePath =  'C:\Users\BMClab\Downloads\Desiree\Acelerometro\Delsys\acc\';
Sensor = 2;

AccX = ReadEMGWorks([FilePath,Name,'.xls'],Sensor,'ACC X');
AccY = ReadEMGWorks([FilePath,Name,'.xls'],Sensor,'ACC Y');
AccZ = ReadEMGWorks([FilePath,Name,'.xls'],Sensor,'ACC Z');
% fsr1 = ReadEMGWorks([FilePath,Name,'.xls'],11,'EMG');
fsr2 = ReadEMGWorks([FilePath,Name,'.xls'],11,'ACC X');
% fsr3 = ReadEMGWorks([FilePath,Name,'.xls'],11,'ACC Y');
% fsr4 = ReadEMGWorks([FilePath,Name,'.xls'],11,'ACC Z');

AccDelsys = [AccX(:,1),sqrt(AccX(:,2).^2 + AccY(:,2).^2 + AccZ(:,2).^2)];
AccDelsys(:,2) = AccDelsys(:,2);%/max(AccDelsys(:,2));

% fsr1(:,2) = fsr1(:,2)/max(fsr1(:,2));
fsr2(:,2) = fsr2(:,2)/max(fsr2(:,2));
% fsr3(:,2) = fsr3(:,2)/max(fsr3(:,2));
% fsr4(:,2) = fsr4(:,2)/max(fsr4(:,2));

figure;
plot(AccDelsys(:,1), AccDelsys(:,2),fsr2(:,1),(fsr2(:,2)+1))
xlim([0 10]); %ylim([0 10]);
legend('ACC Delsys','FSR - Channel 2')
title(['Sensor ',num2str(Sensor)])
