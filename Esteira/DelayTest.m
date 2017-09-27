Name = 'acc51';


FilePath =  'C:\Users\BMClab\Downloads\Desiree\Acelerometro\Delsys\acc\';

F = importdata([FilePath,Name,'.forces']);
Acc = importdata([FilePath,Name,'.anc'],'\t',11);
Contato1 = ReadDelsys([FilePath,Name,'-Delsys 1.csv'], 'EMG', 'EMG');
Contato2 = ReadDelsys([FilePath,Name,'-Delsys 1.csv'], 'ACC', 'ACC Pitch');
AccPitch = ReadDelsys([FilePath,Name,'-Delsys 1.csv'], 'AUX', 'IM ACC Pitch');
AccRoll = ReadDelsys([FilePath,Name,'-Delsys 1.csv'], 'AUX', 'IM ACC Roll');
AccYaw = ReadDelsys([FilePath,Name,'-Delsys 1.csv'], 'AUX', 'IM ACC Roll');

AccDelsys = [AccPitch(:,1),AccPitch(:,3),AccRoll(:,3),AccYaw(:,3)];

AccK = [Acc.data(:,1),Acc.data(:,29:31)];
Fy = [(F.data(:,1)-1)/1200,F.data(:,24)];


Acc = [AccK(:,1),sqrt(AccK(:,2).^2+AccK(:,3).^2+AccK(:,4).^2)];
Acc(:,2) = Acc(:,2)/max(Acc(:,2));

AccD = [AccDelsys(:,1),sqrt(AccDelsys(:,2).^2+AccDelsys(:,3).^2+AccDelsys(:,4).^2)];
AccD(:,2) = AccD(:,2);%/max(AccD(:,2));

Contato2(:,12) = Contato2(:,12)/max(Contato2(:,12));

figure;
plot(Acc(:,1), Acc(:,2:end),AccD(:,1), 6*AccD(:,2:end),Contato2(:,1),(Contato2(:,12)))
legend('kistler','delsys','Contato 2')

% figure;
% plot(AccD(:,1), 6*AccD(:,2:end),Fy(:,1),1*Fy(:,2),Contato2(:,1),1*(Contato2(:,12)))
% legend('AccPitch','AccRoll','AccYaw','FP Kistler','Contato 2')


2.645-2.687
4.6	-4.634
6.223-6.257