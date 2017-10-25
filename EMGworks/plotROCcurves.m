%clc; clear all; close all
% % % load([Path1,'RCdelay_9var_Desiree.mat'])

len = 511;
Path1 = 'C:\Users\BMClab\Downloads\Desiree\Drop Foot Project\Identificacao de evento\Processamento\Resultados\';
load([Path1,'FeaturesComb.mat']);

pct = -0.5: 0.01 : 1;
% %
% % load([Path1,'RCdelay_13var_Desiree.mat'])
% % len = 8191;
% % load([Path1,'FeaturesComb1.mat']);


for sensor = 1:2
    if sensor == 1
        load('Resultados\RNWsensor1_velocities2.mat')
    else
        load('Resultados\RNWsensor2_velocities2.mat')
    end
    [Sensib,Specif] = ROCcurve(RC);
    
    figure(sensor); cmap = colormap(jet(len));
    xlim([0 1]); ylim = ([0 1]);
    xlabel('Taxa de falso positivo');
    ylabel('Taxa de verdadeiro positivo');
    hold on
    
    Q = zeros(len,1);
    for i = 1:len
        ind = i:len:len*length(pct);
        Q(i) = trapz(1-Specif(ind),Sensib(ind));
        plot(1-Specif(ind),Sensib(ind),'color',cmap(i,:))%,'.','MarkerSize',12);
        %
    end
    
    line([0 1], [0 1],'LineStyle','--','color','k','linewidth',2);
    colormap(jet(len))
    c = colorbar('XtickLabel',string(1:floor(len/10):len));
    c.Label.String = 'Índice das combinações';
    c.Label.FontSize = 12;
    
    
    Order = repmat(struct('Index',{0},'Area',{0},'Features',{0}),len,1);
    
    [values,orderInd] = sort(abs(Q));
    
    n=0;
    for i = len:-1:1
        n = n+1;
        Order(n).Index = orderInd(i);
        Order(n).Area = values(i);
        Order(n).Features = Combinations(orderInd(i)).Features;
        
    end
    
    if sensor == 1
        Order1 = Order;
    else
        Order2 = Order;
    end
%     load('RNWsensor2.mat')
end

%load([Path1,'MelhoresCombinacoes_9var.mat']);


