function [] = Main_STP_bar(sensor)

% ------ 1,5, and 20 Hz only ----------------

% sensor='XXX';

B.pr_1HZ = load('csv/pr1HZ.csv');
B.prBLOCKER_1HZ = load('csv/prBLOCKER1HZ.csv');

B.pr_5HZ = load('csv/pr5HZ.csv');
B.prBLOCKER_5HZ = load('csv/prBLOCKER5HZ.csv');

B.pr_20HZ = load('csv/pr20HZ.csv');
B.prBLOCKER_20HZ = load('csv/prBLOCKER20HZ.csv');

B.pr_50HZ = load('csv/pr50HZ.csv');
B.prBLOCKER_50HZ = load('csv/prBLOCKER50HZ.csv');


B.pr_1HZ_raw = load('csv/pr1HZ_raw.csv');
B.prBLOCKER_1HZ_raw = load('csv/prBLOCKER1HZ_raw.csv');

B.pr_5HZ_raw = load('csv/pr5HZ_raw.csv');
B.prBLOCKER_5HZ_raw = load('csv/prBLOCKER5HZ_raw.csv');

B.pr_10HZ_raw = load('csv/pr10HZ_raw.csv');
B.prBLOCKER_10HZ_raw = load('csv/prBLOCKER10HZ_raw.csv'); 

B.pr_20HZ_raw = load('csv/pr20HZ_raw.csv');
B.prBLOCKER_20HZ_raw = load('csv/prBLOCKER20HZ_raw.csv');

B.pr_50HZ_raw = load('csv/pr50HZ_raw.csv');
B.prBLOCKER_50HZ_raw = load('csv/prBLOCKER50HZ_raw.csv');

fsz=10;

%  First Spike Pr -----------------------------------------

mean_ACSF_1    = (B.pr_1HZ_raw(1) + B.pr_5HZ_raw(1) + B.pr_20HZ_raw(1))/3;

mean_BLOCKER_1 = (B.prBLOCKER_1HZ_raw(1) + B.prBLOCKER_5HZ_raw(1) + B.prBLOCKER_20HZ_raw(1))/3;

fprintf('\n Mean Pr for 1st spikes: ACSF %0.2f,  BLOCKER %0.2f \n', mean_ACSF_1, mean_BLOCKER_1); 

fprintf('\n Mean Pr for 1st spike at 1,5,20 Hz: ACSF %0.2f,  %0.2f ,  %0.2f \n', B.pr_1HZ_raw(1), B.pr_5HZ_raw(1), B.pr_20HZ_raw(1) );

%%%%%%%%%%%%%%%%%%   2010   %%%%%%%%%%%%%%%%%%
h1=figure(2010);
clf;
hold on;

mainTitle = sprintf('%s Calcium Sensor Model', sensor);
sgtitle({mainTitle ' '  ' '});

% annotation('textbox', [0 0.9 1 0.1], ...
%     'String', mainTitle, ...
%     'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'center',...
%     'FontSize', 12, 'FontWeight', 'bold');

%%%%%%%%%%%%%%%%%%%% 1 Hz %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subplot(3,1,1);
hold on;  grid on;
t=title('1 Hz,  (1000 isi)');
set(t, 'FontSize', fsz);
bar_handle=bar(1:length(B.pr_1HZ),[B.pr_1HZ' B.prBLOCKER_1HZ'],  'EdgeColor','none');
set(bar_handle(1),'FaceColor',[0,0,1])
set(bar_handle(2),'FaceColor',[1,0,0])

ylabel('Normalised Pr'); 
% ylim([0 1]);
maxY=max([ 400 B.pr_1HZ B.prBLOCKER_1HZ] );
ylim([0,maxY]);
xlim([0, length(B.pr_1HZ)+1]);


%%%%%%%%%%%%%%%%%%%%% 5 Hz
subplot(3,1,2);
% p = get(h, 'position');
% p(2) = p(2) - 0.02;
% set(h, 'position', p);
grid on;   hold on;
t=title('5 Hz,  (200 isi)');
set(t, 'FontSize', fsz);

bar_handle=bar(1:length(B.pr_5HZ),[B.pr_5HZ' B.prBLOCKER_5HZ'],  'EdgeColor','none');
set(bar_handle(1),'FaceColor',[0,0,1])
set(bar_handle(2),'FaceColor',[1,0,0])

ylabel('Normalised Pr'); 
maxY=max([ 400 B.pr_5HZ B.prBLOCKER_5HZ] );
% ylim([0 1]);
ylim([0,maxY]);
xlim([0, length(B.pr_5HZ)+1]);


%%%%%%%%%%%%%%%%%%%%%  20 Hz
subplot(3,1,3);
% p = get(h, 'position');
% p(2) = p(2) - 0.055;
% set(h, 'position', p);
hold on; grid on;
t=title('20 Hz,  (50 isi)');
set(t, 'FontSize', fsz);
bar_handle=bar(1:length(B.pr_20HZ),[B.pr_20HZ' B.prBLOCKER_20HZ'],  'EdgeColor','none');
set(bar_handle(1),'FaceColor',[0,0,1])
set(bar_handle(2),'FaceColor',[1,0,0])

ylabel('Normalised Pr'); 
% ylim([0 1]);
maxY=max([ 400 B.pr_20HZ B.prBLOCKER_20HZ] );
ylim([0,maxY]);
xlim([0, length(B.pr_20HZ)+1]);
xlabel('Spike number');  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%   2000   %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h2=figure(2000);
clf;

mainTitle = sprintf('%s Calcium Sensor Model', sensor);
sgtitle({mainTitle ' '  ' '});

% annotation('textbox', [0 0.9 1 0.1], ...
%     'String', mainTitle, ...
%     'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'center',...
%     'FontSize', 12, 'FontWeight', 'bold');

%%%%%%%%%%%%%%%%%%%% 1 Hz 
subplot(3,3,1);
hold on;  grid on;
t=title('1 Hz,  (1000 isi)');
set(t, 'FontSize', fsz);

bar_handle=bar(1:length(B.pr_1HZ),[B.pr_1HZ' B.prBLOCKER_1HZ'],  'EdgeColor','none');
set(bar_handle(1),'FaceColor',[0,0,1])
set(bar_handle(2),'FaceColor',[1,0,0])

ylabel('Normalised Pr'); 
% ylim([0 1]);
maxY=max([ 400 B.pr_1HZ B.prBLOCKER_1HZ] );
ylim([0,maxY]);
xlim([0, length(B.pr_1HZ)+1]);
xlabel('Spike number');


%%%%%%%%%%%%%%%%%%%%% 5 Hz
subplot(3,3,2);
% p = get(h, 'position');
% p(2) = p(2) - 0.02;
% set(h, 'position', p);

grid on;   hold on;
t=title('5 Hz,  (200 isi)');
set(t, 'FontSize', fsz);

bar_handle=bar(1:length(B.pr_5HZ),[B.pr_5HZ' B.prBLOCKER_5HZ'],  'EdgeColor','none');
set(bar_handle(1),'FaceColor',[0,0,1])
set(bar_handle(2),'FaceColor',[1,0,0])

ylabel('Normalised Pr'); 
maxY=max([ 400 B.pr_5HZ B.prBLOCKER_5HZ] );
% ylim([0 1]);
ylim([0,maxY]);
xlim([0, length(B.pr_5HZ)+1]);
xlabel('Spike number');


%%%%%%%%%%%%%%%%%%%%%  20 Hz
subplot(3,3,3);
% p = get(h, 'position');
% p(2) = p(2) - 0.055;
% set(h, 'position', p);

hold on; grid on;
t=title('20 Hz,  (50 isi)');
set(t, 'FontSize', fsz);

bar_handle=bar(1:length(B.pr_20HZ),[B.pr_20HZ' B.prBLOCKER_20HZ'],  'EdgeColor','none');
set(bar_handle(1),'FaceColor',[0,0,1])
set(bar_handle(2),'FaceColor',[1,0,0])

ylabel('Normalised Pr'); 
% ylim([0 1]);
maxY=max([ 400 B.pr_20HZ B.prBLOCKER_20HZ] );
ylim([0,maxY]);
xlim([0, length(B.pr_20HZ)+1]);
xlabel('Spike number');


%%%%%%%%%%%%%   Pr(50) / Pr(1)
h=subplot(3,3,4);
% p = get(h, 'position');
% p(2) = p(2) - 0.06;
% set(h, 'position', p);

hold on;  grid on;
t=title({'';'Pr(5)/Pr(1)'});
set(t, 'FontSize', fsz);

plot(2, B.pr_1HZ(5)/B.pr_1HZ(1) , 'sb', 'MarkerFaceColor', 'b');  hold on;
% errorbar(2,B.pr_1HZ(5)/B.pr_1HZ(1), stdErrorPr1HZ, 'b', 'LineWidth',1 );


plot(2, B.prBLOCKER_1HZ(5)/B.prBLOCKER_1HZ(1), '*r', 'MarkerFaceColor', 'r');
%errorbar(2,B.prBLOCKER_1HZ(5)/B.prBLOCKER_1HZ(1), stdErrorPrBLOCKER1HZ, 'r', 'LineWidth',1 );


plot(3, B.pr_5HZ(5)/B.pr_5HZ(1) , 'sb', 'MarkerFaceColor', 'b');  hold on;
% errorbar(3,B.pr_5HZ(5)/B.pr_5HZ(1), stdErrorPr5HZ, 'b', 'LineWidth',1 );
% 
plot(3, B.prBLOCKER_5HZ(5)/B.prBLOCKER_5HZ(1), '*r', 'MarkerFaceColor', 'r');
% errorbar(3,B.prBLOCKER_5HZ(5)/B.prBLOCKER_5HZ(1),stdErrorPrBLOCKER5HZ, 'r', 'LineWidth',1 );
% 
% %%%%%%%%%%%
% 
plot(4, B.pr_20HZ(5)/B.pr_20HZ(1) , 'sb', 'MarkerFaceColor', 'b');  hold on;
% % errorbar(4,B.pr_20HZ(5)/B.pr_20HZ(1),stdErrorPr20HZ, 'b', 'LineWidth',1 );
% 
plot(4, B.prBLOCKER_20HZ(5)/B.prBLOCKER_20HZ(1) , '*r', 'MarkerFaceColor', 'r');  hold on;
% errorbar(4,B.prBLOCKER_20HZ(5)/B.prBLOCKER_20HZ(1),stdErrorPrBLOCKER20HZ, 'r', 'LineWidth',1 );
% 
plot(1, 0);
plot(5, 0);

x = [0 1 5 20];
set(gca,'XTick', [1,2,3,4]);
set(gca,'XTickLabel',x);
xlabel('Hz');
%  
% plot(5, B.pr_20HZ(5)/B.pr_20HZ(1) , 'sk', 'MarkerFaceColor', 'b');  hold on;
% plot(5, B.prBLOCKER_20HZ(5)/B.prBLOCKER_20HZ(1), 'sk', 'MarkerFaceColor', [.5 .5 .5]);

ylim([0 4]);


%%%%%%%%%%%%%%%%%%%%  Amount of facilitation
subplot(3,3,5);
hold on; grid on;
title('Mean Facilitation for spikes 2:10');

stdErrorPr1HZ=std(B.pr_1HZ(2:end))/mean(B.pr_1HZ(2:end));
stdErrorPrBLOCKER1HZ=std(B.prBLOCKER_1HZ(2:end))/mean(B.prBLOCKER_1HZ(2:end));

stdErrorPr5HZ=std(B.pr_5HZ(2:end))/mean(B.pr_5HZ(2:end));
stdErrorPrBLOCKER5HZ=std(B.prBLOCKER_5HZ(2:end))/mean(B.prBLOCKER_5HZ(2:end));

m1a =  mean(B.pr_1HZ(2:end));
m1b = mean(B.prBLOCKER_1HZ(2:end));

m5a =  mean(B.pr_5HZ(2:end));
m5b = mean(B.prBLOCKER_5HZ(2:end));

m20a = mean(B.pr_20HZ(2:end));
m20b = mean(B.prBLOCKER_20HZ(2:end));

m50a = mean(B.pr_50HZ(2:end));
m50b = mean(B.prBLOCKER_50HZ(2:end));

bar_handle=bar([1 2 3 4],[ [m1a m5a m20a m50a]' [m1b m5b m20b m50b]' ] ,  'EdgeColor','none');
set(bar_handle(1),'FaceColor',[0,0,1])
set(bar_handle(2),'FaceColor',[1,0,0])

x=[1 5 20 50];
set(gca,'XTick',[1 2 3 4] );
set(gca,'XTickLabel',x);
xlabel('Hz');
%ylim([0 1]);
ylim([0 400]);


stdErrorPr20HZ=std(B.pr_20HZ(2:end))/mean(B.pr_20HZ(2:end));
stdErrorPrBLOCKER20HZ=std(B.prBLOCKER_20HZ(2:end))/mean(B.prBLOCKER_20HZ(2:end));



%%%%%%%%%%%%%%%%%%%%  Mean Pr for spikes 1 to n.
subplot(3,3,6);
hold on; grid on;
title('Mean Pr for spikes 8:10');

m1a = mean(B.pr_1HZ_raw([8,9,10]));
% m1b = mean(B.prBLOCKER_1HZ_raw([8,9,10]));

m5a = mean(B.pr_5HZ_raw([8,9,10]));
% m5b = mean(B.prBLOCKER_5HZ_raw([8,9,10]));

m10a = mean(B.pr_10HZ_raw([8,9,10]));
% m10b = mean(B.prBLOCKER_10HZ_raw([8,9,10]));

m20a = mean(B.pr_20HZ_raw([8,9,10]));
% m20b = mean(B.prBLOCKER_20HZ_raw([8,9,10]));

m50a = mean(B.pr_50HZ_raw([8,9,10]));
% m50b = mean(B.prBLOCKER_50HZ_raw([8,9,10]));

% bar_handle=bar([1 2 3 4],[ [m1a m5a m20a m50a]' [m1b m5b m20b m50b]' ] ,  'EdgeColor','none');
% bar_handle=bar([1 2 3 4],[ [m1a m5a m10a m20a m50a]' ] ,  'EdgeColor','none');
% set(bar_handle(1),'FaceColor',[0,0,1])
% set(bar_handle(2),'FaceColor',[1,0,0])
x=[1 5 20 50];
bar(x, [m1a m5a m20a m50a] );

% set(gca,'XTick',[1 2 3 4] );
% set(gca,'XTickLabel',x);
xlabel('Hz');
ylabel('Pr');
xlim([0 52]);
ylim([0 1]);



subplot(3,3,7)
title('Mean Pr at each 1 Hz Spike');  hold on; grid on;
plot(B.pr_1HZ_raw);  
plot(B.prBLOCKER_1HZ_raw, '--r');
ylim([0 1]);
ylabel('Pr','rot',0);    % set(get(gca,'YLabel'),'Rotation',270);
xlabel('Spike number');

subplot(3,3,8);
title('Mean Pr at each 5 Hz Spike');  hold on; grid on;
plot(B.pr_5HZ_raw);   
plot(B.prBLOCKER_5HZ_raw, '--r');
ylim([0 1]);
xlabel('Spike number');

subplot(3,3,9);
title('Mean Pr at each 20 Hz Spike');  hold on; grid on;
plot(B.pr_20HZ_raw);   
plot(B.prBLOCKER_20HZ_raw, '--r');
ylim([0 1]);
xlabel('Spike number');


%%%% Save Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = sprintf('png/%d-sensor-%s.png', 2010, sensor);

saveas(h1,filename);

filename = sprintf('png/%d-sensor-%s.png', 2000, sensor);

saveas(h2,filename);

% print('-dpdf','-r300', filename);

% cmd = sprintf('pdfcrop -margins 10 %s', filename);
% 
% system(cmd);  





