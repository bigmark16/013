function [] = cantenna_rti_v5(wavFile)
% function [] = cantenna_rti_v4(wavFile)
% 
% Produces an RTI (range x time intensity) image of the
% cantenna recording. Also applies a simple two-pulse 
% canceller to filter out clutter and CFAR normalization
% to improve visual detection. See inside the function
% for additional parameters.
%
%    wavFile = the filename of the .WAV file to process
%    
% MIT IAP 2016 Laptop Radar Course
% (c) 2016 Massachusetts Institute of Technology

% ---- setup constants and parameters ----
c = 299e6; % (m/s) speed of light
Tp = 10e-3; % (s) minimum pulse length
numPad = 64; % number of samples to pad for bandlimited interpolation & shifting
ovsTrg = 16; % oversampling factor applied when interpolating the trigger signal
ovsRng = 2; % oversampling factor applied when interpolating the range data
fStart = 2400e6; % (Hz) LFM start frequency
fStop = 2480e6; % (Hz) LFM stop frequency
nPulseCancel = 2; % number of pulses to use for canceller 
maxRange = 300; % (m) maximum range to display
doTracking = true; % attempt to detect and track targets (can be slow) - Tracker tuning parameters can be set in RadarTracker.m
% ----- end constants and parameters -----

% use a default filename if none is given
if ~exist('wavFile','var')
    wavFile = 'running_outside_20ms.wav';
end

% read the raw wave data
fprintf('Using %g MHz bandwidth\n', (fStop-fStart)*1e-6);
fprintf('Loading WAV file...\n');
[Y,Fs] = audioread(wavFile,'native');

% derived parameters
Np = round(Tp * Fs); % # of samples per pulse
BW = fStop - fStart; % (Hz) transmit bandwidth
delta_r = c/(2*BW); % (m) range resolution

% change sign of the input because it appears inverted in practice
trig = -Y(:,2); % the trigger signal is in the first channel
s = -Y(:,1); % the mixer output is in the second channel
clear Y;


% Clean up the sync signal.  
% Run sync through a schmitt comparator with levels of +threshold and
% -threshold.  Output is signal with +1,-1,0 levels.  Fiddle with threshold
% if needed 
cleantrig = trig;
locked = 0;
threshold = 1000;
fprintf('Preprocess Sync...\n');
for ii = 1:length(trig)
   if (trig(ii) > threshold)
       cleantrig(ii) = 1; locked = 1; 
   elseif (trig(ii) < -threshold)
       cleantrig(ii) = -1; locked = 1;
   else 
       if(locked)
           cleantrig(ii) = cleantrig(ii-1);
       else
           cleantrig(ii) = 0;
       end
   end
    
end
trig = cleantrig;

% parse the trigger signal (look for threshold crossings)
fprintf('Parsing the recording...\n');
[pulseStarts,a] = process_sync(trig);

% refine using measured parameters: the pulse width
Np = round(min(diff(pulseStarts))/2);
Tp = Np / Fs;
fprintf('Measured pulse width of %g ms \n', Tp*1e3);

% pre-compute some windows and other vectors 
Nrange = floor(ovsRng*Np/2); % number of output range samples
dataRange = (0:Nrange-1).' * (delta_r/ovsRng); % ranges of each bin (m)
dataRange = dataRange(dataRange <= maxRange); % apply range limits
Nrange_keep = numel(dataRange); % number of range bins to keep
%
rngWin = hann_window(Np); % the window applied to reduce range sidelobes
padWin = sin( (1:numPad).'/(numPad+1) * pi/2) .^2; % the window applied to the padded data
trgWin = hann_window(numPad*2+1); % the window applied to the trigger data

% get parameters about the data, like how many pulses there are
nSamples = numel(s);
pulseStarts = pulseStarts(pulseStarts+Np+numPad <= nSamples);
numPulses = numel(pulseStarts);
fprintf('Found %d pulses\n',numPulses);

% process pulses into a data matrix
sif = zeros(Nrange_keep,numPulses);
fprintf('Processing pulse data...\n');
for pIdx = 1:numPulses
    % bandlimited interpolate the trigger signal
    tmp = double(trig(pulseStarts(pIdx) + (-numPad:numPad))) .* trgWin;
    interpTmp = fft_interp(tmp,ovsTrg);
    interpTmp = interpTmp( (numPad*ovsTrg + 1) + (-ovsTrg:ovsTrg) );
    interpOffs = (-ovsTrg:ovsTrg)/ovsTrg;
    myIdx = find(diff(sign(interpTmp))==2)+1;
    if isempty(myIdx)
        disp(['Trigger edge not found... skipping pulse ' num2str(pIdx)]);
    end
    tmp2 = interpTmp( myIdx + (-1:0) );
    % linear interpolate to find the zero crossing
    fracOffset = -(interpOffs(myIdx) - tmp2(2)/(tmp2(2)-tmp2(1)) / ovsTrg);
    
    % time-align the data to the trigger event (the zero crossing)
    cInds = pulseStarts(pIdx) + (-numPad:(Np+numPad-1));
    tmp = double(s(cInds));
    tmp(1:numPad) = tmp(1:numPad) .* padWin;
    tmp(end:-1:(end-numPad+1)) = tmp(end:-1:(end-numPad+1)) .* padWin;
    % time delay applied in the frequency domain below
    tmp = fft(tmp);
    tmp = tmp .* exp( -1j*(0:(Np+2*numPad-1)).'/(Np+2*numPad)*2*pi*fracOffset );
    tmp = ifft(tmp,'symmetric');
    
    % compute & scale range data from the time-aligned mixer output
    tmp = ifft(tmp(numPad + (1:Np)) .* rngWin, 2*Nrange);
    sif(:,pIdx) = tmp(1:Nrange_keep);
end
%
clear s trig;
%
sif = sif.';

% display the RTI
figure;
imagesc(dataRange,(0:numPulses-1)*Tp*2,20*log10(abs(sif)));
ylabel('Time (s)');
xlabel('Range (m)');
title('RTI without clutter rejection');
colormap(jet(256));
colorbar;

% apply the N-pulse canceller
mti_filter = -ones(nPulseCancel,1)/nPulseCancel;
midIdx = round((nPulseCancel+1)/2);
mti_filter(midIdx) = mti_filter(midIdx) + 1;
sif = convn(sif,mti_filter,'same');

% display the MTI results
figure;
imagesc(dataRange,(1:numPulses)*Tp*2,20*log10(abs(sif)));
ylabel('Time (s)');
xlabel('Range (m)');
title('RTI with MTI clutter rejection');
colormap(jet(256));
colorbar;

% apply the median CFAR normalization
sif_dB = 20*log10(abs(sif));
sif_dB = sif_dB - repmat(median(sif_dB,1),[size(sif,1) 1]); % over time
sif_dB = sif_dB - repmat(median(sif_dB,2),[1 size(sif,2)]); % over range

% plot the CFAR normalized results
figure;
imagesc(dataRange,(1:numPulses)*Tp*2,sif_dB);
ylabel('Time (s)');
xlabel('Range (m)');
title('RTI with MTI+CFAR');
colormap(jet(256));
caxis([0 40]-3);
colorbar;

% apply tracker and detector
if doTracking
    trackOptions.buffer_size=numPulses;
    trackS = RadarTracker(trackOptions);
    for iPulse = 1:numPulses
        trackS.track_targets(sif_dB(iPulse,:).', iPulse*2*Tp, dataRange);
    end
    
    % plot the CFAR normalized results with detections and tracks overlayed
    figure;
    imagesc(dataRange,(1:numPulses)*Tp*2,sif_dB);
    ylabel('Time (s)');
    xlabel('Range (m)');
    title('RTI Detections and Tracks');
    colormap(jet(256));
    caxis([0 40]-3);
    colorbar;
    hold on
    trackS.plotAllTracksAndDets(gca,(1:numPulses)*Tp*2, false);
end

% ---- standard DSP helper functions below ----

function [y] = fft_interp(x,M)
% perform approximate bandlimited interpolation of x by a factor of M
L = 4;
winInds = (-L*M : L*M).'/M * pi;

% get the ideal antialiasing filter's impulse response of length 2*M + 1 
winInds(L*M + 1) = 1;
myWin = sin(winInds) ./ winInds;
myWin(L*M + 1) = 1;

% use the window method; apply a hann window
myWin = myWin .* hann_window(2*L*M + 1);

% insert zeros in data and apply antialias filter via FFT
nFFT = numel(x) * M;
if isreal(x)
    y = ifft( fft(myWin,nFFT) .* repmat(fft(x),[M 1]), 'symmetric');
else
    y = ifft( fft(myWin,nFFT) .* repmat(fft(x),[M 1]) );
end
y = y([L*M+1:end 1:L*M]);

function [w] = hann_window(N)
% create a hann (cosine squared) window
w = .5 + .5*cos(2*pi*((1:N).'/(N+1) - .5));

function [rise,fall,group_num]=process_sync(x,param)
% identify starts and stops of pulses and pulse groups

% create default param structure if not input
if nargin<2
    display('process_sync using default parameters')
    param.samples_per_sec=44100;
    param.height_thresh=.5;
    param.t_min_pulse=15e-3; % minimum duration of a pulse
    param.t_max_pulse=25e-3; % maximum duration of a pulse
    param.t_min_break=0.5; % minimum duration of gap between pulse groups
    param.n_min_pulses_per_group=25; % minimum number of pulses per group
    param.n_discard_first=1; % drop this many pulses at the beginning of a group
    param.n_discard_last=1; % drop this many pulses at the end of a group
end

% get parameters from the param structure
samples_per_sec=param.samples_per_sec;
height_thresh=param.height_thresh;
t_min_pulse=param.t_min_pulse;
t_max_pulse=param.t_max_pulse;
t_min_break=param.t_min_break;
n_min_pulses_per_group=param.n_min_pulses_per_group;
n_discard_first=param.n_discard_first;
n_discard_last=param.n_discard_last;
    
cleanup=true; % save memory
be_verbose=true;

pulsewidth_thresh=floor(t_min_pulse*samples_per_sec);
breakwidth_thresh=floor(t_min_break*samples_per_sec);

% find when signal has been above a threshold for some duration
x1=x>=height_thresh;
x2=conv(double(x1),ones(pulsewidth_thresh,1),'same');
if cleanup, clear x1, end
x3=x2>=pulsewidth_thresh;
if cleanup, clear x2, end
x4=[0;diff(x3)];
if cleanup, clear x3, end

% initial rise and fall estimates 
rise0=find(x4==1);
fall0=find(x4==-1);
if cleanup, clear x4, end
assert(all(fall0>rise0))

% zero crossings of original signal
x5=x>0;

% this will assume if the first sample is high, it's a rise, and if the
% last sample is high, then the following sample is a fall
x6=diff([false;x5;false]); 
if cleanup, clear x5, end

% candidate rise and fall times (zero crossings)
rise_candidates=find(x6==1);
fall_candidates=find(x6==-1);
if cleanup, clear x6, end
assert(all(fall_candidates>rise_candidates))

% for each rise and fall estimate, find nearest candidate
% (this crashes if there's only one rise or fall candidate)
rise=interp1(rise_candidates,rise_candidates,rise0,'nearest','extrap');
fall=interp1(fall_candidates,fall_candidates,fall0,'nearest','extrap');

assert(~any(isnan(rise)|isnan(fall)))

% this can happen with small positive signals because rise0 is based on
% (non-zero) threshold crossings, and rise_candiates is based on zero
% crossings
ind_discard=(rise>rise0)|(fall<fall0);
rise(ind_discard)=[];
fall(ind_discard)=[];

assert(all(fall>rise))

if be_verbose,display(['found ',num2str(length(rise)),' raw pulses!']),end

% filter on pulsewidth
delta=fall-rise;
ind_discard=(delta<t_min_pulse*samples_per_sec)|(delta>t_max_pulse*samples_per_sec);
rise(ind_discard)=[];
fall(ind_discard)=[];
if be_verbose,display(['discarded ',num2str(nnz(ind_discard)),' pulses based on pulsewidth']),end

% find which group each pulse is in
if length(rise)>1
    is_first_of_group=[true;diff(rise)>=breakwidth_thresh];
else
    is_first_of_group=true(length(rise),1);
end

% filter on number of pulses per group here
first_pulse_ind=find(is_first_of_group);
last_pulse_ind=[first_pulse_ind(2:end)-1;length(is_first_of_group)];
pulses_per_group=1+last_pulse_ind-first_pulse_ind;
ind_discard_group=(pulses_per_group<n_min_pulses_per_group);
ind_discard=false(size(rise));
for k=1:length(pulses_per_group)
    if ind_discard_group(k)
        ind_discard(first_pulse_ind(k):last_pulse_ind(k))=true;
    end
end
rise(ind_discard)=[];
fall(ind_discard)=[];
is_first_of_group(ind_discard)=[];
if be_verbose,display(['discarded ',num2str(nnz(ind_discard)),' pulses based on pulses per group']),end

% find which group each pulse belongs to
group_num=cumsum(is_first_of_group);

% throw away first and last n pulses in each group
temp=ones(n_discard_first+n_discard_last,1);
ind_discard=logical(conv(double(is_first_of_group),temp));
ind_discard(1:n_discard_last)=[]; % first n_discard_last pts
ind_discard((end-n_discard_first+2):end)=[]; % last n_discard_first-1 pts
ind_discard=logical(ind_discard);

rise(ind_discard)=[];
fall(ind_discard)=[];
group_num(ind_discard)=[];
if be_verbose,display(['discarded ',num2str(nnz(ind_discard)),' pulses at beginnings and ends of groups']),end

if be_verbose,display(['final answer: ',num2str(length(rise)),' pulses in ',num2str(group_num(end)),' groups!']),end
return

