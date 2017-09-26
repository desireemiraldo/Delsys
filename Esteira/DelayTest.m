FilePath =  'C:\Users\BMClab\Downloads\Desiree\Acelerometro\Delsys\acc\acc51-Delsys 1.csv';

Contato1 = ReadDelsys(FilePath, 'EMG', 'EMG');
Contato2 = ReadDelsys(FilePath, 'ACC', 'ACC Pitch');

F = importdata('C:\Users\BMClab\Downloads\Desiree\Acelerometro\Delsys\acc\acc51.forces');
Acc = importdata('C:\Users\BMClab\Downloads\Desiree\Acelerometro\Delsys\acc\acc51.anc','\t',11);

AccK = [Acc.data(:,1),Acc.data(:,29:31)];
Fy = [(F.data(:,1)-1)/1200,F.data(:,24)];

plot(AccK(:,1), AccK(:,2:end))
hold on
plot(Fy(:,1),10000*Fy(:,2))
plot(Contato2(:,1),5000*Contato2(:,12)+100000)
