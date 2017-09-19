clc; clear all; close all
% load('RCdelay_9var_Desiree.mat')
% len = 511;

load('RCdelay_13var_Desiree.mat')
len = 8191;

[Sensib,Specif] = ROCcurve(RC);

figure; cmap = colormap(jet(len));
xlim([0 1]); ylim = ([0 1]);
xlabel('Taxa de falso positivo');
ylabel('Taxa de verdadeiro positivo');
hold on

Q = zeros(len,1);
for i = 1:len
    ind = i:len:len*101;
    Q(i) = trapz(1-Specif(ind),Sensib(ind));
    plot(1-Specif(ind),Sensib(ind),'color',cmap(i,:))%,'.','MarkerSize',12);
%     
end

line([0 1], [0 1],'LineStyle','--','color','k','linewidth',2);
colormap(jet(len))
c = colorbar('XtickLabel',string(1:floor(len/10):len));
c.Label.String = 'Índice das combinações';
c.Label.FontSize = 12;


load('FeaturesComb1.mat');

Order = repmat(struct('Index',{0},'Area',{0},'Features',{0}),len,1);

[values,orderInd] = sort(abs(Q));

n=0;
for i = len:-1:1
   n = n+1;
   Order(n).Index = orderInd(i);
   Order(n).Area = values(i);
   Order(n).Features = Combinations(orderInd(i)).Features;
    
end

    
