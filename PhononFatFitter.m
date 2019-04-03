In=double(Linescan_SubtractLG);
SizeIn = size(In);
Out = In;
Number = (1:SizeIn(1));
n = [8,17,26,31];
PlotBuffer=10;
Windows = [Channel0+55,Channel0+200];

% Linescan_HAADF=double(J04_LS8_HAADF);
% Linescan_HAADFSub=Linescan_HAADF(1,:)-min(Linescan_HAADF(1,:));
% Linescan_HAADFNorm=Linescan_HAADFSub(1,:)/max(Linescan_HAADFSub(1,:));

Amplitude = zeros(sizex,1);
Position = zeros(sizex,1);
Variance = zeros(sizex,1);

for m = 1:sizex

% % Comment/Uncomment this to switch on/off a 1 pixel median filter    
% sizeen=int16(sizeen);
% In(m,1)= median(In(m,1:2));
% for n = 2:sizeen-1
%     In(m,n)= median([In(m,n),In(m,n-1),Int(m,n+1)]);% Pick out the mth spectrum for the curve fitting. 
% end
% In(m,sizeen)= median(In(m,sizeen-1:sizeen));
%     
% Comment/Uncomment this to switch to no oversampling    

fittingwin=In(m,Windows(1):Windows(2));

% Rawwin1=LinescanRaw(m,Windows(1):Windows(2));
% Rawwin2=LinescanRaw(m,Windows(3):Windows(4));
% Rawwin = cat(2,Rawwin1,Rawwin2); %Concatenates the values from the two windows. 

fittingEN = EN(Windows(1):Windows(2));

%Fit Gaussian plus Lorentzian to the peaks.
func= @(C,x)C(1)*exp(-(((x-C(2))/C(3)).^2))+C(4);
C1=[  max(In(m,:)), 0.060, 0.020, 0]; % Initial guesses for the fit parameters
lb=[             0, 0.040, 0.008, min(In(m,:))-20*abs(min(In(m,:)))]; % Lower bounds for the fit parameters
ub=[2*max(In(m,:)), 0.080, 0.040, 20*abs(min(In(m,:)))+10]; % Upper bounds for the fit parameters 
%FR = linspace(fitEN(1),fitEN(end));
ft = lsqcurvefit(func,C1,fittingEN,fittingwin, lb, ub, options); % Fit Curve to data.

Out(m,:) = ft(1)*exp(-(((EN(1,:)-ft(2))/ft(3)).^2))+ft(4);
Amplitude(m,1) = ft(1);
Position(m,1) = ft(2);
Variance(m,1) = ft(3);


end
close all;

NormInt = IntegratedZLP(:,1)-min(IntegratedZLP(:,1));
NormInt = NormInt(:,1)/max(NormInt(:,1));

figure(1);
subplot(2,2,1)
plot(EN(1,:),In(n(1),:),'ro',EN(1,:),Out(n(1),:),'b-');
xlim([EN(Windows(1)-PlotBuffer) EN(Windows(2)+PlotBuffer)]);
subplot(2,2,2)
plot(EN(1,:),In(n(2),:),'ro',EN(1,:),Out(n(2),:),'b-');
xlim([EN(Windows(1)-PlotBuffer) EN(Windows(2)+PlotBuffer)]);
subplot(2,2,3)
plot(EN(1,:),In(n(3),:),'ro',EN(1,:),Out(n(3),:),'b-');
xlim([EN(Windows(1)-PlotBuffer) EN(Windows(2)+PlotBuffer)]);
subplot(2,2,4)
plot(EN(1,:),In(n(4),:),'ro',EN(1,:),Out(n(4),:),'b-');
xlim([EN(Windows(1)-PlotBuffer) EN(Windows(2)+PlotBuffer)]);


figure(2);
subplot(3,1,1)
plot(Number,NormInt(:,1)*max(Amplitude(:,1)),'g-',Number,Linescan_HAADFNorm*max(Amplitude(:,1)),'b-');
title('Integrated ZLP Intensity vs HAADF');
subplot(3,1,2)
plot(Number,Position(:,1),'r-');
title('Centre position (meV)');
subplot(3,1,3)
plot(Number,2.355*Variance(:,1),'r');
title('Fitted FWHMs (meV)');
