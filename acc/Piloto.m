FilePath = 'C:\Users\BMClab\Downloads\Desiree\Acelerometro\Delsys\Piloto\';
Name = 'Piloto_profRenato__Plot_and_Store_Rep_1.5';

Sensor = 1;

AccX = ReadEMGWorks([FilePath,Name,'.xls'],Sensor,'ACC X');
AccY = ReadEMGWorks([FilePath,Name,'.xls'],Sensor,'ACC Y');
AccZ = ReadEMGWorks([FilePath,Name,'.xls'],Sensor,'ACC Z');

EMG = ReadEMGWorks([FilePath,Name,'.xls'],Sensor,'EMG');
% -- Filter
EMGenvelope = abs(EMG(:,2));
Fs = 1111.11;
[t,Wn] = buttord(6/(Fs/2),20/(Fs/2),3,60);
[b,a] = butter(t, Wn);

EMGenvelope = filtfilt(b,a,EMGenvelope);



% fsr1 = ReadEMGWorks([FilePath,Name,'.xls'],11,'EMG');
fsr2 = ReadEMGWorks([FilePath,Name,'.xls'],11,'ACC X');
fsr3 = ReadEMGWorks([FilePath,Name,'.xls'],11,'ACC Y');

Sensor = 2;
AccX2 = ReadEMGWorks([FilePath,Name,'.xls'],Sensor,'ACC X');
AccY2 = ReadEMGWorks([FilePath,Name,'.xls'],Sensor,'ACC Y');
AccZ2 = ReadEMGWorks([FilePath,Name,'.xls'],Sensor,'ACC Z');

AccDelsys = [AccX(:,1),sqrt(AccX(:,2).^2 + AccY(:,2).^2 + AccZ(:,2).^2)];
AccDelsys(:,2) = AccDelsys(:,2);%/max(AccDelsys(:,2));

figure;
subplot(2,1,1)
plot(fsr2(:,1),fsr2(:,2),fsr3(:,1),fsr3(:,2), ...
    AccX(:,1),AccX(:,2),AccY(:,1),AccY(:,2),AccZ(:,1),AccZ(:,2),...
    AccX2(:,1),AccX2(:,2),AccY2(:,1),AccY2(:,2),AccZ2(:,1),AccZ2(:,2));
xlim([0 20]);
legend('Heel','Toe','AccX1','AccY1','AccZ1','AccX2','AccY2','AccZ2')
subplot(2,1,2)
plot(EMG(:,1),EMGenvelope)
xlim([0 20]);


AccDelsys2 = [AccX2(:,1),sqrt(AccX2(:,2).^2 + AccY2(:,2).^2 + AccZ2(:,2).^2)];
AccDelsys2(:,2) = AccDelsys2(:,2);%/max(AccDelsys(:,2));
figure;
plot(AccDelsys(:,1),AccDelsys(:,2))


DIF = AccDelsys2 - AccDelsys;