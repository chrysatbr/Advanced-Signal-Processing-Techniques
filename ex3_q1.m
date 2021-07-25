clear
clc
close all

% Load real discrete signal X(k)
load('SignalX.mat');


% Signal Generation
    N = 2048;
    v = exprnd(1,1,N);
    v = v - mean(v);
    q=[1 .93 .85 .72 .59 -.1];
    x = filter(q,1,v);

% Parametric Estimation of Skewness to check for Non-Gaussianity
    skew = 0;
    m = mean(v);
    s = std(v);
    for i=1:N
       skew = skew + ((v(i)-m)^3) ;
    end

    skew = skew/((N-1)*(s^3));
    error = (abs(skew - skewness(v))/abs(skewness(v)))*100;
    str = ['v[k] is Non-Gaussian with estimated skewness = ' ,num2str(skew)];
    str2 = ['Approximation Error = ',num2str(error),'%'];

    if skew ~= 0
        disp(str);
        disp(str2);
    else
        disp('v[k] is Gaussian');
    end
    
% Calculate the skewness of V(k)
sk=skewness(V);

% 3rd order cumulants of x[k] using the indirect method 
K = 32;
M = 64;
L = 20;
c3=Cumulants3(X,L,K,M);

% 3d plot
axisX=-20:20;
axisY=-20:20;
figure
contour(axisX,axisY,c3)
title('3rd Order Cummulants hosa');
figure;
surf(axisX,axisY,c3)
title('3rd order cumulants of x[k]')

% Estimate the impulse response of the MA system using the Giannakis formula
q=5;
hest=Giannakis(c3,q);

% Sub-Sup estimate using Giannakis formula
qsub=q-2;
hsub=Giannakis(c3,qsub);

qsup=q+3;
hsup=Giannakis(c3,qsup);

% Estimate the MA-q system output and calculate the NRMSE
[nrmse,Xest]=NRMSE(hest,V,N,X);

% Repeat the latter, for sub-estimation and over-estimation
[nrmsesub,Xsubest]=NRMSE(hsub,V,N,X);
[nrmsesup,Xsupest]=NRMSE(hsup,V,N,X);

% Repeat the above but instead of x[k] use the noise infected output y[k]
snr=(30:-5:-5); % SNR values
y=zeros(length(snr),N);
hesty=zeros(length(snr),q+1);
nrmsey=zeros(1,length(snr));

for i=1:length(snr)
    y(i,:)=awgn(X,snr(i),'measured');
    c3y=Cumulants3(y(i,:),L,K,M);
    hesty(i,:)=Giannakis(c3y,q);
    set(0,'DefaultFigureVisible','off');
    [nrmsey(i),~]=NRMSE(hesty(i,:),V,N,y(i,:));
end

set(0,'DefaultFigureVisible','on')
figure;
plot(snr,nrmsey);
title('NRMSE of y vs SNR range')
xlabel('SNR(dB)')
ylabel('NRMSE')
    
