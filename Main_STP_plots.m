function [] = Main_STP_plots(isi, ex, sensor)
% ---- plots for one stimulation frequency -----
% http://dgleich.github.io/hq-matlab-figs/
% 
% isi=99;
% ex='test';
% sensor = 'Test';
h=figure(isi);
clf;
mainTitle = sprintf('%d isi, %s, %s Calcium Sensor', isi, ex, sensor);

% annotation('textbox', [0 0.9 1 0.1], ...
%     'String', mainTitle, ...
%     'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'center',...
%     'FontSize', 12, 'FontWeight', 'bold');
sgtitle({mainTitle ' '  ' '});

B.v    = load('csv/bv.csv');

tme=1:length(B.v);
tme=tme./(1000/0.05);

B.c    = load('csv/bc.csv');
B.Gsyn = load('csv/bg.csv');

B.ca_vgcc = load('csv/b_ca_VGCC.csv');
B.ca_nmdaR = load('csv/b_ca_PreNMDAR.csv');
B.ca_nmdaR_mean = load('csv/b_ca_PreNMDAR_mean.csv');

B.ca_md = load('csv/b_ca_MD.csv');
B.ca_md_BLOCKER = load('csv/b_ca_MD_BLOCKER.csv');

% B.ca_md_vgcc  = load('csv/b_ca_MD_vgcc.csv');
% B.ca_md_ryr   = load('csv/b_ca_MD_ryr.csv');


B.ca_ryr = load('csv/b_ca_RyR.csv');
B.cer = load('csv/b_cer.csv');

% from the current simulation at X Hz

% at spike times, for bar chart
B.pr = load('csv/pr.csv');
B.prBLOCKER = load('csv/prBLOCKER.csv');

B.ca_vgcc_ryr = load('csv/b_ca_vgcc_ryr.csv');

% continuous
B.ves.Prel         = load('csv/b_ves_P_release.csv');
B.ves.Prel_BLOCKER = load('csv/b_ves_P_release_BLOCKER.csv');

FREQ = fix(1000/isi);

file = sprintf('csv/pr%dHZ_raw.csv', FREQ);
B.pr_mean = load(file);
file = sprintf('csv/prBLOCKER%dHZ_raw.csv', FREQ);
B.prBLOCKER_mean = load(file);

fsz = 8;      % Fontsize

xmax=10*isi/1000;   % 10 + 2 spikes worth of padding

subplot(4,3,1);  hold on;
t=title('Bouton Vm');
set(t, 'FontSize', fsz);

plot(tme, B.v, 'k');

% xlim([0 xmax]);
ylabel('mV'); 

% B.nmdaR.s1a=load('csv/bs1a.csv');
% B.nmdaR.s1b=load('csv/bs1b.csv');
% 
subplot(4,3,2);   hold on;
t=title('Bouton [Ca^{2+}] global');
set(t, 'FontSize', fsz);

plot(tme, B.c/1000, 'b');
% xlim([0 xmax]);
ylabel('uM');
grid on;


subplot(4,3,3);    hold on;
t=title('Synaptic Glutamate');
set(t, 'FontSize', fsz);

plot(tme, B.Gsyn, 'g');
% xlim([0 xmax]);
ylabel('mM');


subplot(4,3,4);    hold on;
t=title('VGCC [Ca^{2+}] at vesicles');
set(t, 'FontSize', fsz);

plot(tme, B.ca_vgcc/1000, 'b');
% xlim([0 xmax]);
ylabel('uM');
grid on;

subplot(4,3,5);   hold on;
t=title('PreNMDAR [Ca^{2+}] at vesicles');
set(t, 'FontSize', fsz);

plot(tme, B.ca_nmdaR/1000, 'b');
% xlim([0 xmax]);
ylabel('uM');
grid on;

subplot(4,3,6);   hold on;
t=title('[Ca^{2+}] at vesicles (ND)');
set(t, 'FontSize', fsz);

plot(tme, B.ca_md/1000, 'b');
plot(tme, B.ca_md_BLOCKER/1000, '--r');
% xlim([0 xmax]);
ylabel('uM');
grid on;

subplot(4,3,7);   hold on;
t=title('RyR [Ca^{2+}] at vesicles');
set(t, 'FontSize', fsz);

plot(tme, B.ca_ryr/1000, 'b');
% xlim([0 xmax]);
ylabel('uM');
ylim([0, 2]);
grid on;

subplot(4,3,8);   hold on;
t=title('Mean preNMDAR [Ca^{2+}] for vesicles');
set(t, 'FontSize', fsz);

plot(tme, B.ca_nmdaR_mean/1000, 'b');
% xlim([0 xmax]);
ylabel('uM');
grid on;

subplot(4,3,9);  hold on;
title('Sensor activation'); 
set(t, 'FontSize', fsz);

plot(tme, B.ves.Prel, '.-k');  
plot(tme, B.ves.Prel_BLOCKER, '--r'); 

% set(t, 'FontSize', fsz);
% ylabel('mM');
ym= max( max(B.ves.Prel), 1.0);
ylim([0 ym]);
% xlim([0 xmax]);
grid on;


subplot(4,3,10);  hold on;
t=title('[Ca^{2+}] in ER');
set(t, 'FontSize', fsz);

plot(tme, B.cer/1000, 'k');
% xlim([0 xmax]);
xlabel('time (s)');
ylabel('uM');

subplot(4,3,11);  hold on;
title(' Mean Pr '); 
set(t, 'FontSize', fsz);

plot(B.pr_mean, '-k');
plot(B.pr_mean, 'o b');
plot(B.prBLOCKER_mean, '--r');

xlabel('spike number');
ylim([0 1]);
grid on;

subplot(4,3,12);   hold on;
t=title('Mean For Each Spike'); 
set(t, 'FontSize', fsz);

bar_handle=bar(1:length(B.pr),[B.pr' B.prBLOCKER'],  'EdgeColor','none');
set(bar_handle(1),'FaceColor',[0,0,1])
set(bar_handle(2),'FaceColor',[1,0,0])

ylabel('Normalised Pr'); 
maxY=max([ 400 B.pr B.prBLOCKER] );
ylim([0,maxY]);

xlim([0, length(B.pr)+1]);
xlabel('spike number');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = sprintf('png/%d-sensor-%s.png', isi, sensor);

saveas(h,filename);

vg=max(B.ca_vgcc);
ry=max(B.ca_ryr);
vg_ry=max(B.ca_vgcc_ryr);

pre=max(B.ca_nmdaR);
md=mean(B.ca_md);
mx=max(B.ca_md);

% fprintf('max: VGCC=%0.1f,  Ry=%0.1f,  vg+ry=%0.1f \n', vg,ry,vg_ry);
% fprintf('max: preNMDAR=%f \n', pre);
% fprintf('max: MD=%f,  mean MD=%f \n', mx,md);
