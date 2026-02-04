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
%----------------------------------------------------------------------------------------------------
% Author      :   Zhifeng Chen
% E-mail      :   zfchen@whu.edu.cn
% Created     :   2026-01-20  
%----------------------------------------------------------------------------------------------------

%  Data simulation (real sequences)
clc;
clear;
close all;
%
DeltaT  =   30/365.25;
t       =   1982.5:DeltaT:2023.5;    
N       =   length(t);      
%           Amp             frequency
y       =   1 * cos(2 * pi * 2.712 * t)+...
            1 * cos(2 * pi * 2.723 * t)+...
            1 * cos(2 * pi * 2.734 * t)+...
            1 * cos(2 * pi * 2.745 * t)+...
            1 * cos(2 * pi * 2.756 * t)+...
            1 * cos(2 * pi * 2.767 * t)+...
            1 * cos(2 * pi * 2.778 * t)+...
            1 * cos(2 * pi * 2.789 * t)+...
            1 * cos(2 * pi * 2.800 * t);% \Deltaf = 0.011
% Perform a Fourier transform to obtain a Fourier spectrum
R               =   480;
f0              =   (0:R*N-1)/DeltaT/N/R;
FFTf            =   f0-f0(fix(R*N/2)+1);
FFTp            =   abs(fftshift(fft((y.').*hann(N),R*N)))*2/N*2; 
% Perform SRSA to obtain SR spectrum
Deltaf          =   1/(DeltaT*N*R); % The resolution is oversampled by a factor of 240.
bins            =   2.6:Deltaf:3.0;         
[Bins,LVEC,SRSp]=   DoSRSA(y,DeltaT,bins,10,220);     
pp6             =   log(SRSp(25,:));
% Intercept the same interval as shown in 2D SRS to superimpose it on the original image.
[~, idx1]       =   min(abs(Bins - 2.70)); 
[~, idx2]       =   min(abs(Bins - 2.81)); 
Bins_6          =   Bins(idx1:idx2);
pp6             =   pp6(idx1:idx2);
% normalization
pp6             =   (pp6-min(pp6))./(max(pp6)-min(pp6));
% Plot
figure(3);tiledlayout(2,2);
nexttile(1)
plot(t,y,'k-');
xlabel('Time (year)');
ylabel('Amplitude');
legend('Beat Sequence','box','off');
nexttile(2)
plot(FFTf,FFTp,'k-');hold on
xlim([2.5 3]);ylim([0,5]); 
xlabel('Frequency (Hz)');
ylabel('Amplitude');
legend('Fourier Spectrum','box','off');
nexttile(3)
pcolor(Bins,LVEC(:,2),log(SRSp+1));
shading interp;
colormap('Pink');
colorbar("northoutside");
xlim([2.7 2.81]);
ylim([230 250]);
yticklabels([]);
ylabel('2D SRS');
nexttile(4)
plot(Bins_6,pp6,'k-');
xlim([2.7 2.81]);ylim([0,1.2]);
ylabel('Norm SRS')
legend('logarithmic SR spectrum','box','off');

