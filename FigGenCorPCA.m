%% wrapper to generate figures

clear
clc
close all

load r_5.mat
load r_10.mat
load r_30.mat
load r_50.mat

figure
subplot(221)
plot(AlRange_r5, mean(FinalSubspaceError_r5, 1));
hold
plot(AlRange_r5, max(FinalSubspaceError_r5, [], 1), 'r');
plot(AlRange_r5, SE_theory_r5, 'g')
axis tight
xlabel('alpha')
ylabel('SE')
title('r = 5')


subplot(222)
plot(AlRange_r10, mean(FinalSubspaceError_r10, 1));
hold
plot(AlRange_r10, max(FinalSubspaceError_r10, [], 1), 'r');
plot(AlRange_r10, SE_theory_r10, 'g')
axis tight
xlabel('alpha')
ylabel('SE')
title('r = 10')

subplot(223)
plot(AlRange_r30, mean(FinalSubspaceError_r30, 1));
hold
plot(AlRange_r30, max(FinalSubspaceError_r30, [], 1), 'r');
plot(AlRange_r30, SE_theory_r5, 'g')
axis tight
xlabel('alpha')
ylabel('SE')
title('r = 30')


subplot(224)
plot(AlRange_r50, mean(FinalSubspaceError_r50, 1));
hold
plot(AlRange_r50, max(FinalSubspaceError_r50, [], 1), 'r');
plot(AlRange_r50, SE_theory_r50, 'g')
axis tight
xlabel('alpha')
ylabel('SE')
title('r = 50')
