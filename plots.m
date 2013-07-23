set(gcf,'Visible','off');
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultaxesfontSize',8)

nyqFrq = 512;
halfDftPoints = 256/2;

spectro = csvread('spectro.csv');
phase = csvread('phase.csv');
signal = csvread('signal.csv');

plot(signal(:,1), signal(:,2));
ylabel('amplitude');
xlabel('time(sec)');
title('Input signal');
set(gcf,'color','w');
myaa([4 2],'signal.png');

imagesc(spectro);
set(gca,'Xtick',1:10:size(spectro,2),'XTickLabel',round((0:10:size(spectro,2)) / halfDftPoints * nyqFrq * 10)/10);
set(gca,'Ytick',1:35:size(spectro,1),'YTickLabel',round((0:35:size(spectro,1))/size(spectro,1) * 100)/100);
xlabel('frequency bin');
ylabel('time (sec)');
title('Spectrogram');
colorbar
set(gcf,'color','w');
myaa([4 2],'spectrogram.png');

imagesc(phase);
set(gca,'Xtick',1:10:size(spectro,2),'XTickLabel',round((0:10:size(spectro,2)) / halfDftPoints * nyqFrq * 10)/10);
set(gca,'Ytick',1:35:size(spectro,1),'YTickLabel',round((0:35:size(spectro,1))/size(spectro,1) * 100)/100);
xlabel('frequency bin');
ylabel('time (sec)');
title('Phase');
colorbar
set(gcf,'color','w');
myaa([4 2],'phase.png');

