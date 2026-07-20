% *[Overview]
%  This script uses the M1 superconducting-gravimeter record following the 2011 event to demonstrate
%  the practical SRSA workflow and its data-driven application to real observations according to
%  the internal SRSA criteria. Within the predefined target frequency band of the 0S2 mode, neither
%  the number nor the locations of the spectral lines to be resolved are prescribed using theoretical
%  free-oscillation frequencies or previously observed spectral peaks; reported frequencies are used
%  only as references when interpreting the results.
%
%  To reduce the computational cost of the Hankel-SVD and SRSA calculations, the original record y_1,
%  sampled at 60 s, is band-pass filtered and downsampled to 300 s to obtain y_2. Comparison of their
%  Fourier spectra confirms that the principal spectral structure within the 0S2 frequency band is
%  retained after filtering and downsampling. Three adjacent model orders, K_used = 50, 55, and 60,
%  are then selected from the elbow region of the singular-value spectrum of y_2. The consistency of
%  the corresponding SRSA results indicates that moderate variations in the model order do not
%  materially alter the principal spectral conclusions.
%
%  The full 2-D SRS reveals five persistent high-coherence spectral ridges and constitutes the primary
%  result in the frequency-Q_L domain. On this basis, a stable trade-off interval of Q_L is selected,
%  and the spectral responses within this interval are averaged to obtain the 1-D SRS. The 1-D SRS is
%  therefore not an indiscriminate global average of the full 2-D spectrum, but a compact representation
%  of the stable spectral structure identified in the frequency-Q_L plane.
%
% *[Caption for Figure (1)]
%  Figure (1):
%       Panel  (a) y_1 and y_2; (b) singular-value spectrum of y_2;
%       Panels (c,g,m) full 2-D SRS for K_used = 50, 55, and 60;
%       Panels (d,h,n) selected stable-subwindow SRS;
%       Panels (e,i,o) corresponding mean-stacked 1-D SRS;
%       Panel  (f) power spectra of y_1 and y_2;
%       Panels (j,k) normalized power spectra in the 0S2 band;
%       Panel  (l) zoom of the singular-value-spectrum elbow.
%
% *[Computational Environment]
%  The calculations were performed using MATLAB R2024b Update 3 on a workstation equipped with an
%  Intel Core i9-14900K processor (24 cores and 32 threads), 128 GB of RAM, and an NVIDIA GeForce
%  RTX 4090 D GPU with 24 GB of memory. The process-based parallel pool used eight workers. The
%  present implementation does not invoke GPU computing; therefore, the GPU did not contribute to
%  the reported runtimes.
%
% *[Representative Runtime]
%  The total wall-clock time for this demo was approximately 1 min on the stated workstation,
%  including data loading, numerical computation, parallel-pool initialization, and figure
%  generation.
%----------------------------------------------------------------------------------------------------
% Author       :   Zhifeng Chen
% E-mail       :   zfchen@whu.edu.cn
% Created      :   2026-01-20
% Last updated :   2026-06-14
%====================================================================================================
clear;
clc;

%% Data and common settings

load m1_2011.mat

y_1 = series_m1{1}(:);                  % Sampling interval: 60 s
y_2 = series_m1{2}(:);                  % Sampling interval: 300 s

fontName = 'Times New Roman';
fontSize = 12;

modeFreq = [0.299965;0.304536;0.309193;0.313843;0.318433]; % mHz
freqLim  = [0.29,0.33];                                    % mHz
refColor = [0.7,0.7,0.7];

Kused    = [50,55,60];
Kcolor   = {'b','r','g'};
qlRange  = [641,2321];

%% Fourier power spectra

win1  = hann(numel(y_1),'periodic');
win2  = hann(numel(y_2),'periodic');
nfft1 = 2*numel(y_1);
nfft2 = 2*numel(y_2);

[power1,freq1] = periodogram(y_1,win1,nfft1,1/60,'power');
[power2,freq2] = periodogram(y_2,win2,nfft2,1/300,'power');

band1 = freq1>=freqLim(1)*1e-3 & freq1<=freqLim(2)*1e-3;
band2 = freq2>=freqLim(1)*1e-3 & freq2<=freqLim(2)*1e-3;

%% Singular-value spectrum of y_2

N = numel(y_2);
P = fix(N/2);                            % Number of Hankel columns
Q = N-P+1;                              % Number of Hankel rows
H = hankel(y_2(1:Q),y_2(Q:N));

singularV  = svd(H,'econ');
singularSN = (1:numel(singularV)).';

clear H

%% SRSA spectra

DeltaF = 1/(300*N*10);                  % Tenfold frequency-grid oversampling
bins   = 0.270e-3:DeltaF:0.350e-3;

SRSf = cell(1,3);
LVec = cell(1,3);
SRSp = cell(1,3);
stableRows = cell(1,3);

for k = 1:3
    [SRSf{k},LVec{k},SRSp{k}] = DoSRSA( ...
        y_2,300,bins,Kused(k),'AlgorithmMode','svd');

    stableRows{k} = LVec{k}(:,2)>=qlRange(1) & ...
                    LVec{k}(:,2)<=qlRange(2);
end

%% Figure (1)

fig = figure(1);
clf(fig);
tl = tiledlayout(fig,6,6,'TileSpacing','compact','Padding','compact');
colormap(fig,'pink');

% Panel (a): time series
ax = nexttile(tl,1,[2,2]);
plot(ax,1:numel(y_1),y_1,'k-');
hold(ax,'on');
plot(ax,downsample(1:numel(y_1),5),y_2,'r-');
hold(ax,'off');
xlabel(ax,'Time (min)');
ylabel(ax,'Amplitude');
legend(ax, ...
    'Preprocessed and model-corrected record $y_1$', ...
    'Band-pass-filtered and downsampled record $y_2$', ...
    'Interpreter','latex','Box','off');
addPanelLabel(ax,'(a)','k');

% Panel (f): full power spectra
ax = nexttile(tl,13,[2,2]);
plot(ax,freq1*1e3,power1,'k-');
hold(ax,'on');
plot(ax,freq2*1e3,power2,'r--');
hold(ax,'off');
xlim(ax,[0.1,1]);
xlabel(ax,'Frequency (mHz)');
ylabel(ax,'Power (nm/s^2)^2');
legend(ax,'Power spectrum of $y_1$','Power spectrum of $y_2$', ...
    'Interpreter','latex','Box','off');
addPanelLabel(ax,'(f)','k');

% Panels (j) and (k): normalized power spectra in the 0S2 band
bandAxes = [nexttile(tl,25,[2,1]),nexttile(tl,26,[2,1])];
bandFreq = {freq1(band1),freq2(band2)};
bandPower = {power1(band1),power2(band2)};
bandColor = {'k','r'};
bandName = {'$y_1$','$y_2$'};
bandTag = {'(j)','(k)'};

for k = 1:2
    ax = bandAxes(k);
    plot(ax,bandFreq{k}*1e3,rescale(bandPower{k},0,1), ...
        'Color',bandColor{k});
    addModeLines(ax,modeFreq,refColor);
    xlim(ax,freqLim);
    ylim(ax,[0,1.2]);
    xlabel(ax,'Frequency (mHz)');
    ylabel(ax,'Normalized Power');
    legend(ax,string(bandName{k})+" power spectrum in the ${}_0S_2$ band"+newline+ ...
        "$-\,-\,-$ Ref. Ding \& Shen (2013)", ...
        'Interpreter','latex','Box','off');
    addPanelLabel(ax,bandTag{k},'k');
end

% Panel (b): full singular-value spectrum
ax = nexttile(tl,3,[3,1]);
plot(ax,singularSN,singularV,'d','MarkerEdgeColor','k');
ylabel(ax,'Singular Value');
legend(ax,'Singular-value spectrum of $y_2$','Interpreter','latex','Box','off');
addPanelLabel(ax,'(b)','k',[0.06,0.95]);

% Panel (l): singular-value-spectrum elbow
ax = nexttile(tl,21,[3,1]);
p0 = plot(ax,singularSN,singularV,'d','MarkerEdgeColor','k');
hold(ax,'on');
pk = gobjects(1,3);
for k = 1:3
    pk(k) = plot(ax,Kused(k),singularV(Kused(k)),'d', ...
        'MarkerFaceColor',Kcolor{k},'MarkerEdgeColor',Kcolor{k}, ...
        'MarkerSize',8);
end
hold(ax,'off');
xlim(ax,[30,150]);
ylim(ax,[0,100]);
xlabel(ax,'Singular-Value Index');
ylabel(ax,'Singular Value');
legend(ax,[p0,pk], ...
    "Zoom of the elbow region in panel (b)", ...
    '$K_{used}=50$','$K_{used}=55$','$K_{used}=60$', ...
    'Interpreter','latex','Box','off');
addPanelLabel(ax,'(l)','k');

% Full 2-D, selected 2-D, and mean-stacked 1-D SRS panels
fullTiles     = [4,16,28];
selectedTiles = [5,17,29];
stackTiles    = [6,18,30];
fullTags      = {'(c)','(g)','(m)'};
selectedTags  = {'(d)','(h)','(n)'};
stackTags     = {'(e)','(i)','(o)'};

for k = 1:3
    % Full 2-D SRS
    ax = nexttile(tl,fullTiles(k),[2,1]);
    pcolor(ax,SRSf{k}*1e3,LVec{k}(:,2),SRSp{k});
    shading(ax,'interp');
    xlim(ax,freqLim);
    addQLRectangle(ax,qlRange,'data');
    ylabel(ax,'Subwindows $Q_L$','Interpreter','latex');
    if k==3, xlabel(ax,'Frequency (mHz)'); end
    addPanelLabel(ax,fullTags{k},'w');
    text(ax,0.03,0.30, ...
        sprintf('Selected trade-off subwindows $Q_L \\in [%d,%d]$',qlRange), ...
        'Units','normalized','Color','c','Interpreter','latex');
    text(ax,0.40,0.95,sprintf('$K_{used}=%d$',Kused(k)), ...
        'Units','normalized','Color',Kcolor{k},'Interpreter','latex');

    % Selected stable-subwindow SRS
    ax = nexttile(tl,selectedTiles(k),[2,1]);
    pcolor(ax,SRSf{k}*1e3,LVec{k}(stableRows{k},2), ...
        SRSp{k}(stableRows{k},:));
    shading(ax,'interp');
    xlim(ax,freqLim);
    addQLRectangle(ax,qlRange,'axes');
    ylabel(ax,'Subwindows $Q_L$','Interpreter','latex');
    if k==3, xlabel(ax,'Frequency (mHz)'); end
    addPanelLabel(ax,selectedTags{k},'w');
    text(ax,0.40,0.95,sprintf('$K_{used}=%d$',Kused(k)), ...
        'Units','normalized','Color',Kcolor{k},'Interpreter','latex');

    % Mean-stacked 1-D SRS
    ax = nexttile(tl,stackTiles(k),[2,1]);
    plot(ax,SRSf{k}*1e3,mean(SRSp{k}(stableRows{k},:),1), ...
        'Color',Kcolor{k});
    addModeLines(ax,modeFreq,refColor);
    xlim(ax,freqLim);
    ylabel(ax,'Normalized 1-D SRS');
    if k==3, xlabel(ax,'Frequency (mHz)'); end
    legend(ax,string(sprintf('Mean stacking spectrum ($K_{used}=%d$)',Kused(k)))+ ...
        newline+"$-\,-\,-$ Ref. Ding \& Shen (2013)",'Interpreter','latex','Box','off');
    addPanelLabel(ax,stackTags{k},'k');

end

% Apply a common font to all text-bearing graphics objects in one statement.
set(findall(fig,'-property','FontName'), ...
    'FontName',fontName,'FontSize',fontSize);



% % Local plotting utilities

function addPanelLabel(ax,label,color,position)
if nargin<4, position = [0.06,0.92]; end
text(ax,position(1),position(2),label, ...
    'Units','normalized','FontWeight','bold','Color',color);
end


function addModeLines(ax,modeFreq,color)
xline(ax,modeFreq,'--','Color',color,'HandleVisibility','off');
end


function addQLRectangle(ax,qlRange,mode)
xl = xlim(ax);

if strcmp(mode,'data')
    position = [xl(1),qlRange(1),diff(xl),diff(qlRange)];
    clipping = 'on';
else
    yl = ylim(ax);
    position = [xl(1),yl(1),diff(xl),diff(yl)];
    clipping = 'off';
end

hold(ax,'on');
rectangle(ax,'Position',position, ...
    'EdgeColor','c','LineWidth',1.5,'LineStyle','-', ...
    'HandleVisibility','off','Clipping',clipping);
hold(ax,'off');
end

