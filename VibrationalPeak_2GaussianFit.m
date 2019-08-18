In=double(Linescan_SubtractLG);
Scale = 0.192;
SizeIn = size(In);
Out = In;
% startEN = -0.400;
% dispersion2 = 0.0005;
% EN = (startEN:dispersion2:startEN+dispersion2*(SizeIn(2)-1));
% Channel02=-startEN/dispersion2;
Number = (1:SizeIn(1));
Angstroms = (Number-1)*Scale;
Intensities = zeros(SizeIn(1),2);
Positions = zeros(SizeIn(1),2);
Variance = zeros(SizeIn(1),2);
n = [2,9,19,25,27,28];
PlotBuffer=2;
Windows = [Channel0+14,Channel0+48];

Linescan_HAADFSub=Linescan_HAADF(1,:)-min(Linescan_HAADF(1,:));
Linescan_HAADFNorm=Linescan_HAADFSub(1,:)/max(Linescan_HAADFSub(1,:));
NormIntZLP(:,1) = IntegratedZLP(:,1)./max(IntegratedZLP(:,1));

for m = 1:SizeIn(1)

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

%Fit Two Gaussians to the Data.
func= @(C,x)C(1)*exp(-(((x-C(2))/C(3)).^2))+C(4)*exp(-(((x-C(5))/C(6)).^2))+C(7);
C1=[  max(In(m,:)), 0.045, 0.012,   max(In(m,:)), 0.060, 0.012, 0]; % Initial guesses for the fit parameters
lb=[             0, 0.040, 0.008,              0, 0.058, 0.008, min(In(m,:))-20*abs(min(In(m,:)))]; % Lower bounds for the fit parameters
ub=[2*max(In(m,:)), 0.048, 0.020, 2*max(In(m,:)), 0.065, 0.020, 20*abs(min(In(m,:)))+10]; % Upper bounds for the fit parameters 
%FR = linspace(fitEN(1),fitEN(end));
ft = lsqcurvefit(func,C1,fittingEN,fittingwin, lb, ub, options); % Fit Curve to data.

Out(m,:) = ft(1)*exp(-(((EN(1,:)-ft(2))/ft(3)).^2))+ft(4)*exp(-(((EN(1,:)-ft(5))/ft(6)).^2))+ft(7);
Intensities(m,1) = ft(1);
Intensities(m,2) = ft(4);
Positions(m,1) = ft(2);
Positions(m,2) = ft(5);
Variance(m,1) = ft(3);
Variance(m,2) = ft(6);

end
close all;

mEN = EN*1000;

figure(1);
set(gcf, 'color', [1 1 1])
subplot(2,3,1)
plot(mEN(1,:),In(n(1),:),'ro',mEN(1,:),Out(n(1),:),'b-');
xlim([mEN(Windows(1)-PlotBuffer) mEN(Windows(2)+PlotBuffer)]);
ylim([-1000 20000]);
set(gca,'FontSize',16);
set(gca,'YTickLabel',[]);
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 20:5:100;
subplot(2,3,2)
plot(mEN(1,:),In(n(2),:),'ro',mEN(1,:),Out(n(2),:),'b-');
xlim([mEN(Windows(1)-PlotBuffer) mEN(Windows(2)+PlotBuffer)]);
ylim([-1000 20000]);
set(gca,'FontSize',16);
set(gca,'YTickLabel',[]);
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 20:5:100;
subplot(2,3,3)
plot(mEN(1,:),In(n(3),:),'ro',mEN(1,:),Out(n(3),:),'b-');
xlim([mEN(Windows(1)-PlotBuffer) mEN(Windows(2)+PlotBuffer)]);
ylim([-1000 20000]);
set(gca,'FontSize',16);
set(gca,'YTickLabel',[]);
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 20:5:100;
subplot(2,3,4)
plot(mEN(1,:),In(n(4),:),'ro',mEN(1,:),Out(n(4),:),'b-');
xlim([mEN(Windows(1)-PlotBuffer) mEN(Windows(2)+PlotBuffer)]);
ylim([-1000 20000]);
set(gca,'FontSize',16);
set(gca,'YTickLabel',[]);
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 20:5:100;
subplot(2,3,5)
plot(mEN(1,:),In(n(5),:),'ro',mEN(1,:),Out(n(5),:),'b-');
xlim([mEN(Windows(1)-PlotBuffer) mEN(Windows(2)+PlotBuffer)]);
ylim([-1000 20000]);
set(gca,'FontSize',16);
set(gca,'YTickLabel',[]);
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 20:5:100;
subplot(2,3,6)
plot(mEN(1,:),In(n(6),:),'ro',mEN(1,:),Out(n(6),:),'b-');
xlim([mEN(Windows(1)-PlotBuffer) mEN(Windows(2)+PlotBuffer)]);
ylim([-1000 20000]);
set(gca,'FontSize',16);
set(gca,'YTickLabel',[]);
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 20:5:100;

% figure(2);
% set(gcf, 'color', [1 1 1])
% subplot(2,2,1)
% plot(Number,1.7725*Variance(:,1).*Intensities(:,1)./(dispersion*IntegratedZLP(1,1)),'r-');
% title('45 ish meV Intensity');
% subplot(2,2,2)
% plot(Number,Positions(:,1),'r-');
% title('45 ish meV Actual Fitted Energies');
% subplot(2,2,3)
% plot(Number,1.7725*Variance(:,2).*Intensities(:,2)/(dispersion*IntegratedZLP(1,1)),'r-');
% title('60 ish meV Intensity');
% subplot(2,2,4)
% plot(Number,Positions(:,2),'b-');
% title('60 ish meV Actual Fitted Energies');

% figure(2);
% NormInt = IntegratedZLP(:,1)-min(IntegratedZLP(:,1));
% NormInt = NormInt(:,1)/max(NormInt(:,1));

% subplot(2,2,1)
% plot(Number,Intensities(:,1),'r-',Number,Linescan_HAADFNorm*max(Intensities(:,1)),'b-');
% title('45 meV Intensity vs HAADF');
% subplot(2,2,2)
% plot(Number,Intensities(:,1),'r-',Number,NormInt(:,1)*max(Intensities(:,1)),'g-');
% title('45 meV Intensity vs EELS');
% subplot(2,2,3)
% plot(Number,Intensities(:,2),'r-', Number,Linescan_HAADFNorm*max(Intensities(:,2)),'b-');
% title('60 meV Intensity vs HAADF');
% subplot(2,2,4)
% plot(Number,Intensities(:,2),'r-',Number,NormInt(:,1)*max(Intensities(:,2)),'g-');
% title('60 meV Intensity vs EELS');
Phonon45(:,1)=1.7725*Variance(:,1).*Intensities(:,1);
Phonon60(:,1)=1.7725*Variance(:,2).*Intensities(:,2);

Phonon60Sub=Phonon60(:,1)-min(Phonon60(:,1));
Phonon60Norm=Phonon60Sub(:,1)/max(Phonon60Sub(:,1));

Phonon45Sub=Phonon45(:,1)-min(Phonon45(:,1));
Phonon45Norm=Phonon45Sub(:,1)/max(Phonon45Sub(:,1));

figure(2);
set(gcf, 'color', [1 1 1]);
subplot(3,1,1)
plot(Angstroms/10,Phonon45,'r-', Angstroms/10,Phonon60,'b-', 'Linewidth', 2);% Angstroms/10,Linescan_HAADFNorm*max(Phonon60),'k--', 'Linewidth', 2);
% title('Intensities (45ish Red, 60ish Blue, Scaled HAADF Black )');
set(gca,'FontSize',16);
set(gca,'YTickLabel',[]);
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 0:0.1:10;
xlim([0 Angstroms(end)/10]);
subplot(3,1,2)
plot(Angstroms/10,Positions(:,1),'r-',Angstroms/10,Positions(:,2),'b-');
title('Actual Fitted Energies (45ish Red, 60ish Blue)');
xlim([0 Angstroms(end)/10]);
subplot(3,1,3)
plot(Angstroms/10,1.665*Variance(:,1),'r-', Angstroms/10,1.665*Variance(:,2),'b-');
title('Fitted FWHMs (45ish Red, 60ish Blue)' );
xlim([0 Angstroms(end)/10]);

