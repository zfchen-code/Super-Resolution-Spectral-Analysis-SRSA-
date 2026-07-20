%----------------------------------------------------------------------------------------------------
%                    "Information is the resolution of uncertainty."
%                                      — Paraphrased from Claude E. Shannon
%----------------------------------------------------------------------------------------------------
% *[Overview]
%  This script evaluates the frequency-resolution capability of SRSA in the presence of Gaussian
%  white noise. Each simulated time series contains two unit-amplitude sinusoidal components at 3 cpy
%  and 3 + Delta_f cpy. Their phases are selected so that the two components are in phase at the
%  centre of the observation interval.
%
%  The three frequency separations are Delta_f = 0.040, 0.020, and 0.0133 cpy. For the nominal 50-
%  year observation interval, 2*(1/T) = 0.040 cpy is adopted as the nominal main-lobe-width reference.The 
%  three separations therefore correspond approximately to the full nominal main-lobe width, one-half 
%  of that width, and one-third of that width, respectively.
%
%  Different noise levels are added to generate nine groups of simulated time series:
%
%                         { {y11}, {y12}, {y13};
%                           {y21}, {y22}, {y23};
%                           {y31}, {y32}, {y33} }
%
%  The first subscript identifies a progressively smaller frequency separation, whereas the second
%  identifies a progressively higher signal-to-noise ratio.
%
%  The Gaussian white-noise standard deviation is calibrated using the one-sided Fourier amplitude
%  normalization. At any fixed nonzero frequency bin, the pointwise Rayleigh 99.99th-percentile
%  noise amplitude is prescribed to equal 1, 1/5, 1/9, 1/13, or 1/17 of the theoretical amplitude
%  of a single tone. For the nine simulation groups, these ratios are arranged as 1, 1/5, and 1/9
%  for the first frequency separation; 1/5, 1/9, and 1/13 for the second; and 1/9, 1/13, and 1/17
%  for the third.
%
%  One hundred independent Monte Carlo realizations are generated for each group using a fixed
%  random seed for reproducibility. Figure (2) displays the Fourier spectra of all realizations,
%  whereas the mean 2-D SRSA responses and mean Fourier amplitude spectra used in Figures (3)
%  and (4) are obtained by ensemble averaging over the 100 realizations. This averaging reduces
%  realization-specific noise fluctuations and facilitates assessment of whether the two closely
%  spaced components remain consistently distinguishable.
%
% *[Caption for Figures]
%
%  Figure (1): Singular-value spectra of the Hankel matrices constructed from the 100 Monte Carlo
%              realizations in each of the nine simulation groups.
%
%  Figure (2): Fourier amplitude spectra of all 100 Monte Carlo realizations in each
%              simulation group, together with the prescribed pointwise Rayleigh 99.99th-percentile
%              noise-amplitude threshold and the theoretical amplitude of a single tone.
%
%  Figure (3): Ensemble-mean 2-D SRSA coherence responses of the nine simulation groups.
%
%  Figure (4): Detailed comparison between the ensemble-mean SRSA and Fourier
%              spectra.
%
%      Panels (a)-(c): Selected Q_L ranges of the ensemble-mean 2-D SRSA coherence responses for
%                      y11, y12, and y13, respectively.
%      Panels (d)-(f): Comparisons between the normalized stacked 1-D SRSA spectra and normalized
%                      ensemble-mean Fourier amplitude spectra for y11, y12,
%                      and y13, respectively.
%
%      Panels (g)-(i): Selected Q_L ranges of the ensemble-mean 2-D SRSA coherence responses for
%                      y21, y22, and y23, respectively.
%      Panels (j)-(l): Comparisons between the normalized stacked 1-D SRSA spectra and normalized
%                      ensemble-mean Fourier amplitude spectra for y21, y22,
%                      and y23, respectively.
%
%      Panels (m)-(o): Selected Q_L ranges of the ensemble-mean 2-D SRSA coherence responses for
%                      y31, y32, and y33, respectively.
%      Panels (p)-(r): Comparisons between the normalized stacked 1-D SRSA spectra and normalized
%                      ensemble-mean Fourier amplitude spectra for y31, y32,
%                      and y33, respectively.
%
% *[Computational Environment]
%  The calculations were performed using MATLAB R2024b Update 3 on a workstation equipped with an
%  Intel Core i9-14900K processor (24 cores and 32 threads), 128 GB of RAM, and an NVIDIA GeForce
%  RTX 4090 D GPU with 24 GB of memory. The process-based parallel pool used eight workers. The
%  present implementation does not invoke GPU computing; therefore, the GPU did not contribute to
%  the reported runtimes.
%
% *[Representative Runtime]
%  The expected wall-clock time for the complete Monte Carlo experiment is approximately 70 min on
%  the stated workstation. The computational cost is dominated by the singular-value analyses and
%  the rSVD-based SRSA calculations for the 900 simulated realizations.
%----------------------------------------------------------------------------------------------------
% Author       :   Zhifeng Chen
% E-mail       :   zfchen@whu.edu.cn
% Created      :   2026-01-20
% Last updated :   2026-06-20
%====================================================================================================

%% Two-tone simulations with pointwise Rayleigh-0.9999 SNR

clear;
close all;
clc;

rng(20260025,'twister'); 

N      = 6001;
t      = linspace(0,50,N);
DeltaT = t(2)-t(1);
tc     = (t(1)+t(end))/2;

dfList = [0.040,0.020,0.0133];

% Each row corresponds to one frequency separation:
% df=0.0400: SNR9999=[1,5,9]
% df=0.0200: SNR9999=[5,9,13]
% df=0.0133: SNR9999=[9,13,17]
SNR9999Grid = [1,  5,  9; ...
               5,  9, 13; ...
               9, 13, 17];

nMC   = 100;
q9999 = 0.9999;

win = hann(N);
cg  = mean(win);

A_signal = 1;

% Rayleigh scale of the one-sided amplitude spectrum for unit-variance
betaUnit = 2*sqrt(sum(win.^2)/2)/(N*cg);

% Pointwise Rayleigh-0.9999 noise amplitude for unit-variance white noise.
Anoise9999 = betaUnit*sqrt(-2*log(1-q9999));

% y(:,r,j,i): time, realization, SNR index, and frequency separation.
y = zeros(N,nMC,size(SNR9999Grid,2),numel(dfList));

noiseStd9999 = zeros(size(SNR9999Grid));

for i = 1:numel(dfList)

    df = dfList(i);

    phaseShift = -2*pi*df*tc;

    s = cos(2*pi*3*t) + ...
        cos(2*pi*(3+df)*t + phaseShift);

    s = s(:);

    for j = 1:size(SNR9999Grid,2)

        SNR9999 = SNR9999Grid(i,j);
        noiseStd = A_signal/(SNR9999*Anoise9999);

        noiseStd9999(i,j) = noiseStd;
        y(:,:,j,i) = s + noiseStd*randn(N,nMC);

    end
end

% % Original-style aliases

y11 = num2cell(y(:,:,1,1),1).';  % df=0.0400, SNR9999=1
y12 = num2cell(y(:,:,2,1),1).';  % df=0.0400, SNR9999=5
y13 = num2cell(y(:,:,3,1),1).';  % df=0.0400, SNR9999=9

y21 = num2cell(y(:,:,1,2),1).';  % df=0.0200, SNR9999=5
y22 = num2cell(y(:,:,2,2),1).';  % df=0.0200, SNR9999=9
y23 = num2cell(y(:,:,3,2),1).';  % df=0.0200, SNR9999=13

y31 = num2cell(y(:,:,1,3),1).';  % df=0.0133, SNR9999=9
y32 = num2cell(y(:,:,2,3),1).';  % df=0.0133, SNR9999=13
y33 = num2cell(y(:,:,3,3),1).';  % df=0.0133, SNR9999=17

data9 = {y11,y12,y13, ...
         y21,y22,y23, ...
         y31,y32,y33};

%% Singular-value spectra

P = fix(N/2);
L = N+1-P;

SV9 = cell(size(data9));

for k = 1:numel(data9)

    % Column 1 contains singular-value indices.
    % Columns 2:(nMC+1) contain the Monte Carlo realizations.
    S = zeros(P,nMC+1);
    S(:,1) = (1:P).';

    for r = 1:nMC

        x = data9{k}{r}(:);
        H = hankel(x(1:L),x(L:N));

        S(:,r+1) = svd(H,'econ');

        fprintf('Group %d/%d | Realization %d/%d completed\n',k,numel(data9),r,nMC);

    end

    SV9{k} = S;
end

%% Figure (1): singular-value spectra

figure(1);
clf;

set(gcf,'Units','centimeters','Position',[2,2,18,16]);

tiledlayout(3,3,'TileSpacing','tight','Padding','tight');

ax = gobjects(9,1);

for k = 1:numel(SV9)

    ax(k) = nexttile;

    plot(SV9{k}(:,1),SV9{k}(:,2:end));
    xlim([1,50]);

    if ismember(k,[1,4,7])
        ylabel('Singular Value','FontName','Times New Roman','FontSize',10);
    end

    if ismember(k,[7,8,9])
        xlabel('Singular-Value Index','FontName','Times New Roman','FontSize',10);
    end
end

% Add panel labels (a)-(i).
for k = 1:numel(ax)

    text(ax(k),0.025,0.97,sprintf('(%c)','a'+k-1), ...
        'Units','normalized', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'FontName','Times New Roman', ...
        'FontSize',10, ...
        'FontWeight','bold', ...
        'Color','k');

end

%% One-sided Fourier amplitude spectra

Rfft      = 30;
Nfft      = Rfft*N;
nPositive = Nfft/2+1;

freq = (0:Nfft/2).'/(Nfft*DeltaT);

FFT9 = cell(size(data9));

for k = 1:numel(data9)

    X = [data9{k}{:}];


    F   = fft(X.*win,Nfft,1);
    Amp = abs(F(1:nPositive,:))/sum(win);
    Amp(2:end-1,:) = 2*Amp(2:end-1,:);

    % Column 1 contains frequency; the remaining columns contain spectra.
    FFT9{k} = [freq,Amp];

end

%% Figure (2): Fourier spectra of all Monte Carlo realizations

figure(2);
clf;

set(gcf,'Units','centimeters','Position',[2,2,18,14]);

tiledlayout(3,3,'TileSpacing','tight','Padding','tight');

ax = gobjects(9,1);

for k = 1:numel(FFT9)

    row = ceil(k/3);
    col = mod(k-1,3)+1;

    ax(k) = nexttile;

    plot(FFT9{k}(:,1),FFT9{k}(:,2:end));
    ylim([0,1.2]);

    yline(1/SNR9999Grid(row,col),'m--', ...
        'noise-amplitude threshold', ...
        'LineWidth',1, ...
        'FontSize',8);

    yline(1,'g--', ...
        'Theoretical amplitude of a single tone', ...
        'LineWidth',1, ...
        'FontSize',8);

    if ismember(k,[1,4,7])
        ylabel('Fourier Amplitude','FontName','Times New Roman','FontSize',10);
    end

    if ismember(k,[7,8,9])
        xlabel('Frequency (cpy)','FontName','Times New Roman','FontSize',10);
    end
end

% Add panel labels (a)-(i).
for k = 1:numel(ax)

    text(ax(k),0.025,0.97,sprintf('(%c)','a'+k-1), ...
        'Units','normalized', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'FontName','Times New Roman', ...
        'FontSize',10, ...
        'FontWeight','bold', ...
        'Color','k');

end

%% Compute SRSA spectra for all nine Monte Carlo groups
tic
% Use the same frequency spacing as the 30-times zero-padded FFT.
bins       = 2.9:1/(Rfft*N*DeltaT):3.15;

SRSp9     = cell(size(data9));
MeanSRSp9 = cell(size(data9));

for k = 1:numel(data9)

    nRealizations = numel(data9{k});

    for ii = 1:nRealizations

        fprintf('Group %d/%d | Realization %d/%d\n',k,numel(data9),ii,nRealizations);

        [Bins,Q_L,IM] = DoSRSA(data9{k}{ii},DeltaT,bins,4,'AlgorithmMode','rsvd');

        % Allocate after the first calculation to follow the actual output size.
        if ii == 1
            SRSp9{k} = zeros(size(IM,1),size(IM,2),nRealizations,'like',IM);
        end

        SRSp9{k}(:,:,ii) = IM;

    end

    MeanSRSp9{k} = mean(SRSp9{k},3);

end
toc;
%% Original-style aliases

[SRSp11,SRSp12,SRSp13, ...
 SRSp21,SRSp22,SRSp23, ...
 SRSp31,SRSp32,SRSp33] = SRSp9{:};

[MeanSRSp11,MeanSRSp12,MeanSRSp13, ...
 MeanSRSp21,MeanSRSp22,MeanSRSp23, ...
 MeanSRSp31,MeanSRSp32,MeanSRSp33] = MeanSRSp9{:};

% All groups use the same frequency vector and Q_L settings.
[Bins11,Bins12,Bins13, ...
 Bins21,Bins22,Bins23, ...
 Bins31,Bins32,Bins33] = deal(Bins);

[Q_L11,Q_L12,Q_L13, ...
 Q_L21,Q_L22,Q_L23, ...
 Q_L31,Q_L32,Q_L33] = deal(Q_L);

% local rows used to construct the stacked 1-D SRSA spectra.
stackRows = 22:29;
stackQL   = Q_L(stackRows([1,end]),2);

NS = @(S) rescale(S,0,1);

%% Figure (3): full mean 2-D SRSA spectra

figure(3);
clf;

set(gcf,'Units','centimeters','Position',[2,2,16,18]);

tiledlayout(3,3,'TileSpacing','tight','Padding','tight');

colormap('pink');

ax = gobjects(9,1);

% % Row 1: df = 0.0400

% y11: df=0.0400, SNR9999=1
ax(1) = nexttile(1);

pcolor(Bins11,Q_L11(:,2),log(MeanSRSp11+1));
xlim([2.95,3.09]);
shading interp;
colorbar('northoutside','Color','k');

ylabel('$Q_L$', ...
       'Interpreter','latex', ...
       'FontName','Times New Roman', ...
       'FontSize',12);

% y12: df=0.0400, SNR9999=5
ax(2) = nexttile(2);

pcolor(Bins12,Q_L12(:,2),log(MeanSRSp12+1));
xlim([2.95,3.09]);
shading interp;
colorbar('northoutside','Color','k');

% y13: df=0.0400, SNR9999=9
ax(3) = nexttile(3);

pcolor(Bins13,Q_L13(:,2),log(MeanSRSp13+1));
xlim([2.95,3.09]);
shading interp;
colorbar('northoutside','Color','k');

% % Row 2: df = 0.0200

% y21: df=0.0200, SNR9999=5
ax(4) = nexttile(4);

pcolor(Bins21,Q_L21(:,2),log(MeanSRSp21+1));
xlim([2.95,3.07]);
shading interp;
colorbar('northoutside','Color','k');

ylabel('$Q_L$', ...
       'Interpreter','latex', ...
       'FontName','Times New Roman', ...
       'FontSize',12);

% y22: df=0.0200, SNR9999=9
ax(5) = nexttile(5);

pcolor(Bins22,Q_L22(:,2),log(MeanSRSp22+1));
xlim([2.95,3.07]);
shading interp;
colorbar('northoutside','Color','k');

% y23: df=0.0200, SNR9999=13
ax(6) = nexttile(6);

pcolor(Bins23,Q_L23(:,2),log(MeanSRSp23+1));
xlim([2.95,3.07]);
shading interp;
colorbar('northoutside','Color','k');

% % Row 3: df = 0.0133

% y31: df=0.0133, SNR9999=9
ax(7) = nexttile(7);

pcolor(Bins31,Q_L31(:,2),log(MeanSRSp31+1));
xlim([2.95,3.0633]);
shading interp;
colorbar('northoutside','Color','k');

xlabel('Frequency (cpy)','FontName','Times New Roman','FontSize',12);

ylabel('$Q_L$', ...
       'Interpreter','latex', ...
       'FontName','Times New Roman', ...
       'FontSize',12);

% y32: df=0.0133, SNR9999=13
ax(8) = nexttile(8);

pcolor(Bins32,Q_L32(:,2),log(MeanSRSp32+1));
xlim([2.95,3.0633]);
shading interp;
colorbar('northoutside','Color','k');

xlabel('Frequency (cpy)','FontName','Times New Roman','FontSize',12);

% y33: df=0.0133, SNR9999=17
ax(9) = nexttile(9);

pcolor(Bins33,Q_L33(:,2),log(MeanSRSp33+1));
xlim([2.95,3.0633]);
shading interp;
colorbar('northoutside','Color','k');

xlabel('Frequency (cpy)','FontName','Times New Roman','FontSize',12);

% Add white panel labels (a)-(i).
for k = 1:numel(ax)

    text(ax(k),0.025,0.97,sprintf('(%c)','a'+k-1), ...
        'Units','normalized', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'FontName','Times New Roman', ...
        'FontSize',10, ...
        'FontWeight','bold', ...
        'Color','w');

end

%% Figure (4): stacked 1-D SRSA, and Fourier spectra

figure(4);
clf;

tiledlayout(3,6);
colormap('pink');

ax = gobjects(18,1);

% % Row 1: df = 0.0400

% y11: local range
ax(1) = nexttile(1);

pcolor(Bins11,Q_L11(:,2),log(MeanSRSp11+1));
xlim([2.95,3.09]);
ylim([2500,3500]);
shading interp;
colorbar('northoutside','Color','k');

ylabel('$Q_L$', ...
    'Interpreter','latex', ...
    'FontName','Times New Roman', ...
    'FontSize',12);

% y12: local range
ax(2) = nexttile(2);

pcolor(Bins12,Q_L12(:,2),log(MeanSRSp12+1));
xlim([2.95,3.09]);
ylim([2500,3500]);
shading interp;
colorbar('northoutside','Color','k');

% y13: local range
ax(3) = nexttile(3);

pcolor(Bins13,Q_L13(:,2),log(MeanSRSp13+1));
xlim([2.95,3.09]);
ylim([2500,3500]);
shading interp;
colorbar('northoutside','Color','k');

% y11: stacked 1-D SRSA and mean Fourier spectrum
ax(4) = nexttile(4);

p1 = plot(Bins11,NS(mean(MeanSRSp11(stackRows,:),1)));
hold on;
p2 = plot(FFT9{1}(:,1),NS(mean(FFT9{1}(:,2:end),2)));
hold off;

xlim([2.95,3.09]);
ylim([0,1.2]);
xline([3.000,3.040],'m--');

ylabel('Normalized Spectra','FontName','Times New Roman','FontSize',12);

legend([p1,p2], ...
    'Normalized 1-D SRSA', ...
    'Normalized Mean Fourier Spectrum', ...
    'Box','off');

text(0.32,0.10, ...
    sprintf('Stacked $Q_L \\in [%d,%d]$',stackQL(1),stackQL(2)), ...
    'Units','normalized', ...
    'Color',[0,0.4470,0.7410], ...
    'Interpreter','latex');

% y12: stacked 1-D SRSA and mean Fourier spectrum
ax(5) = nexttile(5);

p1 = plot(Bins12,NS(mean(MeanSRSp12(stackRows,:),1)));
hold on;
p2 = plot(FFT9{2}(:,1),NS(mean(FFT9{2}(:,2:end),2)));
hold off;

xlim([2.95,3.09]);
ylim([0,1.2]);
xline([3.000,3.040],'m--');

legend([p1,p2], ...
    'Normalized 1-D SRSA', ...
    'Normalized Mean Fourier Spectrum', ...
    'Box','off');

text(0.32,0.10, ...
    sprintf('Stacked $Q_L \\in [%d,%d]$',stackQL(1),stackQL(2)), ...
    'Units','normalized', ...
    'Color',[0,0.4470,0.7410], ...
    'Interpreter','latex');

% y13: stacked 1-D SRSA and mean Fourier spectrum
ax(6) = nexttile(6);

p1 = plot(Bins13,NS(mean(MeanSRSp13(stackRows,:),1)));
hold on;
p2 = plot(FFT9{3}(:,1),NS(mean(FFT9{3}(:,2:end),2)));
hold off;

xlim([2.95,3.09]);
ylim([0,1.2]);
xline([3.000,3.040],'m--');

legend([p1,p2], ...
    'Normalized 1-D SRSA', ...
    'Normalized Mean Fourier Spectrum', ...
    'Box','off');

text(0.32,0.10, ...
    sprintf('Stacked $Q_L \\in [%d,%d]$',stackQL(1),stackQL(2)), ...
    'Units','normalized', ...
    'Color',[0,0.4470,0.7410], ...
    'Interpreter','latex');

% % Row 2: df = 0.0200

% y21: local range
ax(7) = nexttile(7);

pcolor(Bins21,Q_L21(:,2),log(MeanSRSp21+1));
xlim([2.95,3.07]);
ylim([2500,3500]);
shading interp;
colorbar('northoutside','Color','k');

ylabel('$Q_L$', ...
    'Interpreter','latex', ...
    'FontName','Times New Roman', ...
    'FontSize',12);

% y22: local range
ax(8) = nexttile(8);

pcolor(Bins22,Q_L22(:,2),log(MeanSRSp22+1));
xlim([2.95,3.07]);
ylim([2500,3500]);
shading interp;
colorbar('northoutside','Color','k');

% y23: local range
ax(9) = nexttile(9);

pcolor(Bins23,Q_L23(:,2),log(MeanSRSp23+1));
xlim([2.95,3.07]);
ylim([2500,3500]);
shading interp;
colorbar('northoutside','Color','k');

% y21: stacked 1-D SRSA and mean Fourier spectrum
ax(10) = nexttile(10);

p1 = plot(Bins21,NS(mean(MeanSRSp21(stackRows,:),1)));
hold on;
p2 = plot(FFT9{4}(:,1),NS(mean(FFT9{4}(:,2:end),2)));
hold off;

xlim([2.95,3.07]);
ylim([0,1.2]);
xline([3.000,3.020],'m--');

ylabel('Normalized Spectra','FontName','Times New Roman','FontSize',12);

legend([p1,p2], ...
    'Normalized 1-D SRSA', ...
    'Normalized Mean Fourier Spectrum', ...
    'Box','off');

text(0.32,0.10, ...
    sprintf('Stacked $Q_L \\in [%d,%d]$',stackQL(1),stackQL(2)), ...
    'Units','normalized', ...
    'Color',[0,0.4470,0.7410], ...
    'Interpreter','latex');

% y22: stacked 1-D SRSA and mean Fourier spectrum
ax(11) = nexttile(11);

p1 = plot(Bins22,NS(mean(MeanSRSp22(stackRows,:),1)));
hold on;
p2 = plot(FFT9{5}(:,1),NS(mean(FFT9{5}(:,2:end),2)));
hold off;

xlim([2.95,3.07]);
ylim([0,1.2]);
xline([3.000,3.020],'m--');

legend([p1,p2], ...
    'Normalized 1-D SRSA', ...
    'Normalized Mean Fourier Spectrum', ...
    'Box','off');

text(0.32,0.10, ...
    sprintf('Stacked $Q_L \\in [%d,%d]$',stackQL(1),stackQL(2)), ...
    'Units','normalized', ...
    'Color',[0,0.4470,0.7410], ...
    'Interpreter','latex');

% y23: stacked 1-D SRSA and mean Fourier spectrum
ax(12) = nexttile(12);

p1 = plot(Bins23,NS(mean(MeanSRSp23(stackRows,:),1)));
hold on;
p2 = plot(FFT9{6}(:,1),NS(mean(FFT9{6}(:,2:end),2)));
hold off;

xlim([2.95,3.07]);
ylim([0,1.2]);
xline([3.000,3.020],'m--');

legend([p1,p2], ...
    'Normalized 1-D SRSA', ...
    'Normalized Mean Fourier Spectrum', ...
    'Box','off');

text(0.32,0.10, ...
    sprintf('Stacked $Q_L \\in [%d,%d]$',stackQL(1),stackQL(2)), ...
    'Units','normalized', ...
    'Color',[0,0.4470,0.7410], ...
    'Interpreter','latex');

% % Row 3: df = 0.0133

% y31: local range
ax(13) = nexttile(13);

pcolor(Bins31,Q_L31(:,2),log(MeanSRSp31+1));
xlim([2.95,3.0633]);
ylim([2500,3500]);
shading interp;
colorbar('northoutside','Color','k');

xlabel('Frequency (cpy)','FontName','Times New Roman','FontSize',12);

ylabel('$Q_L$', ...
    'Interpreter','latex', ...
    'FontName','Times New Roman', ...
    'FontSize',12);

% y32: local range
ax(14) = nexttile(14);

pcolor(Bins32,Q_L32(:,2),log(MeanSRSp32+1));
xlim([2.95,3.0633]);
ylim([2500,3500]);
shading interp;
colorbar('northoutside','Color','k');

xlabel('Frequency (cpy)','FontName','Times New Roman','FontSize',12);

% y33: local range
ax(15) = nexttile(15);

pcolor(Bins33,Q_L33(:,2),log(MeanSRSp33+1));
xlim([2.95,3.0633]);
ylim([2500,3500]);
shading interp;
colorbar('northoutside','Color','k');

xlabel('Frequency (cpy)','FontName','Times New Roman','FontSize',12);

% y31: stacked 1-D SRSA and mean Fourier spectrum
ax(16) = nexttile(16);

p1 = plot(Bins31,NS(mean(MeanSRSp31(stackRows,:),1)));
hold on;
p2 = plot(FFT9{7}(:,1),NS(mean(FFT9{7}(:,2:end),2)));
hold off;

xlim([2.95,3.0633]);
ylim([0,1.2]);
xline([3.000,3.0133],'m--');

xlabel('Frequency (cpy)','FontName','Times New Roman','FontSize',12);

ylabel('Normalized Spectra','FontName','Times New Roman','FontSize',12);

legend([p1,p2], ...
    'Normalized 1-D SRSA', ...
    'Normalized Mean Fourier Spectrum', ...
    'Box','off');

text(0.32,0.10, ...
    sprintf('Stacked $Q_L \\in [%d,%d]$',stackQL(1),stackQL(2)), ...
    'Units','normalized', ...
    'Color',[0,0.4470,0.7410], ...
    'Interpreter','latex');

% y32: stacked 1-D SRSA and mean Fourier spectrum
ax(17) = nexttile(17);

p1 = plot(Bins32,NS(mean(MeanSRSp32(stackRows,:),1)));
hold on;
p2 = plot(FFT9{8}(:,1),NS(mean(FFT9{8}(:,2:end),2)));
hold off;

xlim([2.95,3.0633]);
ylim([0,1.2]);
xline([3.000,3.0133],'m--');

xlabel('Frequency (cpy)','FontName','Times New Roman','FontSize',12);

legend([p1,p2], ...
    'Normalized 1-D SRSA', ...
    'Normalized Mean Fourier Spectrum', ...
    'Box','off');

text(0.32,0.10, ...
    sprintf('Stacked $Q_L \\in [%d,%d]$',stackQL(1),stackQL(2)), ...
    'Units','normalized', ...
    'Color',[0,0.4470,0.7410], ...
    'Interpreter','latex');

% y33: stacked 1-D SRSA and mean Fourier spectrum
ax(18) = nexttile(18);

p1 = plot(Bins33,NS(mean(MeanSRSp33(stackRows,:),1)));
hold on;
p2 = plot(FFT9{9}(:,1),NS(mean(FFT9{9}(:,2:end),2)));
hold off;

xlim([2.95,3.0633]);
ylim([0,1.2]);
xline([3.000,3.0133],'m--');

xlabel('Frequency (cpy)','FontName','Times New Roman','FontSize',12);

legend([p1,p2], ...
    'Normalized 1-D SRSA', ...
    'Normalized Mean Fourier Spectrum', ...
    'Box','off');

text(0.32,0.10, ...
    sprintf('Stacked $Q_L \\in [%d,%d]$',stackQL(1),stackQL(2)), ...
    'Units','normalized', ...
    'Color',[0,0.4470,0.7410], ...
    'Interpreter','latex');

% % Panel labels (a)-(r)

pcolorTiles = [1:3,7:9,13:15];

for k = 1:18

    if ismember(k,pcolorTiles)
        labelColor = 'w';
    else
        labelColor = 'k';
    end

    text(ax(k),0.025,0.97,sprintf('(%c)','a'+k-1), ...
        'Units','normalized', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'FontName','Times New Roman', ...
        'FontSize',11, ...
        'FontWeight','bold', ...
        'Color',labelColor);

end
