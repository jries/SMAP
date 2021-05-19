% close all
% clear all
% clc

load('test_gpu_1.mat')
nIterations =int32([
0;
0;
]);

disp('test gpu...');
tic

dme_cuda(coords, crlb, framenum, numspots, maxit, drift, framesperbin, gradientStep, maxdrift, scores, flags, maxneighbors, nIterations);
toc

initial_drift = drift;
initial_scores  = scores ;
initial_nIterations=nIterations(1);

load('test_gpu_2.mat')

disp('test gpu...');
tic

dme_cuda(coords, crlb, framenum, numspots, maxit, drift, framesperbin, gradientStep, maxdrift, scores, flags, maxneighbors, nIterations);
toc

Estimated_drift_DME = drift;
Estimated_scores_DME  = scores ;
Estimated_nIterations_DME=nIterations(1);

% load('drift_trace_1.mat')
load('drift_trace_test.mat')
load('estimated_drift_rcc_test.mat')

figure('name','gpu test reslut')
x_esti_dme = Estimated_drift_DME(1,:)+0.2;
y_esti_dme = Estimated_drift_DME(2,:)+0.2;
z_esti_dme = Estimated_drift_DME(3,:)+0.2;

h_axis = 1:length(x_esti_dme);
x_true = drift_trace(1,:);
y_true = drift_trace(2,:);
z_true = drift_trace(3,:);

x_esti_rcc = estimated_drift_rcc(1,:)-0.2;
y_esti_rcc = estimated_drift_rcc(2,:)-0.2;
z_esti_rcc = estimated_drift_rcc(3,:)-0.2;

subplot(3,1,1);
plot (h_axis,x_true,'b',h_axis,x_esti_dme,'r',h_axis,x_esti_rcc,'g');
ylabel('Drift [px]'),title('x');
legend('True drift','Estimated drift(DME)','Estimated drift(RCC)','Location','northwest');
subplot(3,1,2);
plot (h_axis,y_true,'b',h_axis,y_esti_dme,'r',h_axis,y_esti_rcc,'g');
ylabel('Drift [px]'),title('y');
subplot(3,1,3);
plot (h_axis,z_true,'b',h_axis,z_esti_dme,'r',h_axis,z_esti_rcc,'g');
ylabel('Drift [um]'),title('z');


