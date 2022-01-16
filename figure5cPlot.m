clear;
clc
close all

load('figure5a40Hz.mat')


aa=seis_record;
imagesc([1:1609]*2.5,[150:850],seis_record(:,150:850),[-10^-3 10^-3])
colormap gray
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')

load('figure5b40Hz.mat')
figure;imagesc([1:1609]*2.5,[150:850],seis_record(:,150:850),[-10^-3 10^-3])
colormap gray
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')


figure;imagesc([1:1609]*2.5,[150:850],seis_record-aa,[-10^-3 10^-3])
colormap gray
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')

figure;plot([1:1609]*2.5,-aa(:,400),'k','linewidth',1);
hold on;plot([1:1609]*2.5,-seis_record(:,400),'r','linewidth',1);
legend("E-SGFD Scheme","INB-SGFD Scheme");
axis([1 1609*2.5 -1.855*10^-3 3.2*10^-3])
grid on

xlabel('\it Travel time(ms)')
ylabel('\it Amp')