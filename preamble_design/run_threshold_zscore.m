% Data
y = [1 1 1.1 1 0.9 1 1 1.1 1 0.9 1 1.1 1 1 0.9 1 1 1.1 1 1,...
    1 1 1.1 0.9 1 1.1 1 1 0.9 1 1.1 1 1 1.1 1 0.8 0.9 1 1.2 0.9 1,...
    1 1.1 1.2 1 1.5 1 3 2 5 3 2 1 1 1 0.9 1,...
    1 3 2.6 4 3 3.2 2 1 1 0.8 4 4 2 2.5 1 1 1];

% Settings
lag = 30;
threshold = 5;
influence = 0;

% Get results
[signals,avg,dev] = threshold_zscore(y,lag,threshold,influence);

figure; subplot(2,1,1); hold on;
x = 1:length(y); ix = lag+1:length(y);
area(x(ix),avg(ix)+threshold*dev(ix),'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
area(x(ix),avg(ix)-threshold*dev(ix),'FaceColor',[1 1 1],'EdgeColor','none');
plot(x(ix),avg(ix),'LineWidth',1,'Color','cyan','LineWidth',1.5);
plot(x(ix),avg(ix)+threshold*dev(ix),'LineWidth',1,'Color','green','LineWidth',1.5);
plot(x(ix),avg(ix)-threshold*dev(ix),'LineWidth',1,'Color','green','LineWidth',1.5);
plot(1:length(y),y,'b');
subplot(2,1,2);
stairs(signals,'r','LineWidth',1.5); ylim([-1.5 1.5]);