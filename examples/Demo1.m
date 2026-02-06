%----------------------------------------------------------------------------------------------------
%     "What we observe is not nature itself, but nature exposed to our method of questioning."
%                                                        ——— Werner Heisenberg
%----------------------------------------------------------------------------------------------------
%  *[Overview]
%   This script demonstrates the high-frequency resolution of SRSA.   
%    
%   Since this experiment is performed under a noise-free condition, K_used can be chosen over a  ...
%   wide range (e.g., 50–240, adjust the optimal inspection range accordingly). A value larger    ...
%   than the known K_true is intentionally used to keep the noise subspace as pure as possible.
%   *[Caption for Figure 1]
%       subplot (a): Simulation Sequence y_1.
%       subplot (b): Simulation Sequence y_2.
%       subplot (c): y_1 Fourier spectrum of simulated sequence.
%       subplot (d): y_1 Local 2D SRS of Simulated Sequence.
%       subplot (e): Take the logarithm of the y_1 SR spectrum
%       subplot (f): y_2 Fourier spectrum of simulated sequence.
%       subplot (g): y_2 Local 2D SRS of Simulated Sequence.
%       subplot (h): Take the logarithm of the y_2 SR spectrum
%----------------------------------------------------------------------------------------------------
% Author      :   Zhifeng Chen
% E-mail      :   zfchen@whu.edu.cn
% Created     :   2026-01-20  
%----------------------------------------------------------------------------------------------------
%  Data simulation (real sequences)
clc;
clear;
%
DeltaT  =   30/365.25;
t       =   1982.5:DeltaT:2023.5;    
N       =   length(t);      
%           Amp             frequency
y1      =   1 * cos(2 * pi * 1.180 * t + 1.1)+...
            1 * cos(2 * pi * 1.185 * t + 0.0)+...
            1 * cos(2 * pi * 1.190 * t + 0.0)+...
            1 * cos(2 * pi * 1.195 * t + 0.1)+...
            1 * cos(2 * pi * 1.200 * t + 1.1);         
y2      =   1 * cos(2 * pi * 2.712 * t)+...
            1 * cos(2 * pi * 2.723 * t)+...
            1 * cos(2 * pi * 2.734 * t)+...
            1 * cos(2 * pi * 2.745 * t)+...
            1 * cos(2 * pi * 2.756 * t)+...
            1 * cos(2 * pi * 2.767 * t)+...
            1 * cos(2 * pi * 2.778 * t)+...
            1 * cos(2 * pi * 2.789 * t)+...
            1 * cos(2 * pi * 2.800 * t);% \Deltaf = 0.011
% Perform a Fourier transform to obtain a Fourier spectrum
R                   =   240;
f0                  =   (0:R*N-1)/DeltaT/N/R;
FFTf                =   f0-f0(fix(R*N/2)+1);
FFTp1               =   abs(fftshift(fft((y1.').*hann(N),R*N)))*2/N*2; 
FFTp2               =   abs(fftshift(fft((y2.').*hann(N),R*N)))*2/N*2; 
% Perform SRSA to obtain SR spectrum
Deltaf              =   1/(DeltaT*N*R); % The resolution is oversampled by a factor of 240.
bins1               =   1.17:Deltaf:1.215;      
[Bins1,LVEC1,SRSp1] =   DoSRSA(y1,DeltaT,bins1,10,100);   
pp1                 =   log(SRSp1(11,:));
% Intercept the same interval as shown in 2D SRS to superimpose it on the original image.
[~, idx1]           =   min(abs(Bins1 - 1.17)); 
[~, idx2]           =   min(abs(Bins1 - 1.21)); 
Bins_1              =   Bins1(idx1:idx2);
pp1                 =   pp1(idx1:idx2);
% normalization
pp1                 =   (pp1-min(pp1))./(max(pp1)-min(pp1));
bins2               =   2.6:Deltaf:3.0;         
[Bins2,LVEC2,SRSp2] =   DoSRSA(y2,DeltaT,bins2,10,220);   
pp2                 =   log(SRSp2(25,:));
% Intercept the same interval as shown in 2D SRS to superimpose it on the original image.
[~, idx1]           =   min(abs(Bins2 - 2.70)); 
[~, idx2]           =   min(abs(Bins2 - 2.81)); 
Bins_2              =   Bins2(idx1:idx2);
pp2                 =   pp2(idx1:idx2);
% normalization
pp2                 =   (pp2-min(pp2))./(max(pp2)-min(pp2));
% Plot
figure(1);tiledlayout(2,6);
nexttile(1,[1,3])
plot(t,y1,'k-');
xlabel('Time (year)');
ylabel('Amplitude');
legend('Simulation Sequence-1 (y_1)','box','off');
text(0.02,0.92,'(a)','Units','normalized','FontWeight','bold','FontSize',12);
nexttile(7)
plot(FFTf,FFTp1,'k-');hold on
xlim([1.06 1.31]);ylim([0,5]); 
xlabel('Frequency (Hz)');
ylabel('Amplitude');
legend('y_1 Fourier Spectrum','box','off');
text(0.06,0.92,'(c)','Units','normalized','FontWeight','bold','FontSize',12);
nexttile(8)
pcolor(Bins1,LVEC1(:,2),log(SRSp1+1));
shading interp;
colormap('Pink');
colorbar("northoutside");
xlim([1.17 1.21]);
ylim([51 61]);
yticklabels([]);
ylabel('2D SRS');
legend('y_1 2D SRS','box','off','textcolor','w');
text(0.06,0.92,'(d)','Units','normalized','FontWeight','bold','FontSize',12,'Color','w');
nexttile(9)
plot(Bins_1,pp1,'k-');
xlim([1.17 1.21]);ylim([0,1.2]);
ylabel('Norm SRS')
legend('y_1 logarithmic SR spectrum','box','off');
text(0.06,0.92,'(e)','Units','normalized','FontWeight','bold','FontSize',12);
nexttile(4,[1,3])
plot(t,y2,'k-');
xlabel('Time (year)');
ylabel('Amplitude');
legend('Simulation Sequence-2 (y_2)','box','off');
text(0.02,0.92,'(b)','Units','normalized','FontWeight','bold','FontSize',12);
nexttile(10)
plot(FFTf,FFTp2,'k-');hold on
xlim([2.5 3]);ylim([0,5]); 
xlabel('Frequency (Hz)');
ylabel('Amplitude');
legend('y_2 Fourier Spectrum','box','off');
text(0.06,0.92,'(f)','Units','normalized','FontWeight','bold','FontSize',12);
nexttile(11)
pcolor(Bins2,LVEC2(:,2),log(SRSp2+1));
shading interp;
colormap('Pink');
colorbar("northoutside");
xlim([2.7 2.81]);
ylim([230 250]);
yticklabels([]);
ylabel('2D SRS');
text(0.06,0.92,'(g)','Units','normalized','FontWeight','bold','FontSize',12,'Color','w');
legend('y_2 2D SRS','box','off','textcolor','w');
nexttile(12)
plot(Bins_2,pp2,'k-');
xlim([2.7 2.81]);ylim([0,1.2]);
ylabel('Norm SRS')
legend('y_2 logarithmic SR spectrum','box','off');
text(0.06,0.92,'(h)','Units','normalized','FontWeight','bold','FontSize',12);




