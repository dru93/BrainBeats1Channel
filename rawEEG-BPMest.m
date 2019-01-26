addpath('MATLAB-Tempogram-Toolbox_1.0');
clear
close all

%% parameters
tempoWindow = 12; % window length for tempogram computation
smooth_len = 0.5; % window for local average substraction in novelty curve

sID = '14';       % stimulus    ID

%% configuration: directories, filenames ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dirData = 'data/raw/';
dirAnnoData = 'data/';
filenameWav = [ sID '.wav'];
filenameAnn = [ sID '.txt'];

%% load beat annotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
annoData = csvread( [dirAnnoData filenameAnn] );

% participant IDs
pIDs = [ 01 04 06 07 09 11 12 13 14 ];

% column: participant, row: eeg channel
BPMguess = zeros(length(pIDs), 64);

% iterate through all participants
for pIdx = 1:1%length( pIDs )
    % tranform to '0x' format
    pID = num2str( pIDs(pIdx), '%02d');
    disp('Participant:')
    disp(pID)
    
    filenamePart = ['P' pID '-eegdata.mat'];
    %% load EEG data 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dataEEG = load( [dirData filenamePart]);
    dataEEG = dataEEG.data;
    
    % iterate through all channels
    disp('eeg channels')
    for channelIdx = 3:3%64
        % tranform to '0x' format
        cID = num2str( channelIdx, '%02d');
        dataEEGchannel = dataEEG(channelIdx,:);
        
        disp(cID)

        featureRateEEG = 512; 

        parameterSmooth = [];
        parameterSmooth.smooth_len = smooth_len;

        % Consider accumulatedEEG as novelty curve and subtract local average 
        % -> results in highpass filtering
        parameterSmooth = [];
        parameterSmooth.smooth_len = smooth_len;
        [ noveltyCurveEEG, local_average ] = novelty_smoothedSubtraction_EEG( dataEEGchannel,parameterSmooth );
        
        %Visualize EEG data
        parameterVis = [];
        parameterVis.featureRate = featureRateEEG;
        parameterVis.plotAnn = annoData;

        % visualization of channel with local average curve
        [time_axis] = visualize_noveltyCurve(dataEEGchannel,parameterVis);
        title(['EEG Channel ' num2str(channelIdx), ' ID: ', num2str(sID) ])
        axis tight
        hold on;
        plot(time_axis, local_average,'r','LineWidth',2)

        % visualization of highpass filtered version of channel
        visualize_noveltyCurve(noveltyCurveEEG,parameterVis);
        axis tight
        title(['EEG Channel ' num2str(channelIdx), ' ID: ', num2str(sID) ])

        %% tempogram_fourier (EEG)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        featureRate = featureRateEEG;                  
        noveltyCurve = noveltyCurveEEG;

        parameterTempogram = [];
        parameterTempogram.featureRate = featureRate;
        parameterTempogram.tempoWindow = tempoWindow;   % window length in sec
        parameterTempogram.BPM = 30:1:240;              % tempo values

        [tempogram_fourier, T, BPM] = noveltyCurve_to_tempogram_via_DFT(noveltyCurve,parameterTempogram);
        tempogram_fourier = normalizeFeature(tempogram_fourier,2, 0.0001);

        visualize_tempogram(tempogram_fourier,T,BPM)
        title(['Tempogram (Fourier) on EEG data' pID '-' cID])

        % mean on tempogram
        figure;
        mean_tempogram = mean(abs(tempogram_fourier),2);
        stem(BPM, mean_tempogram)
        xlabel('Tempo (BPM)')
        title(['Mean tempogram (EEG)' pID '-' cID])

        % save figure to myplots folder
        name = ['myplots/tempogramHist-P' pID '-C' cID '.png'];
        saveas(gcf, name)

        [m, i] = max(mean(abs(tempogram_fourier),2));
        BPMguess(pIdx, channelIdx) = BPM(i);
        % close all;
    end
end

csvwrite(['stim-' sID '.csv'], BPMguess)

% BPMguess = load(['stim-47.csv']);
scores = mean(BPMguess, 1);
errors = abs(scores - 160);

[sortedErrors, sortedChannels] = sort(errors,'ascend');

BPMestimation = scores(sortedChannels(1));

plot(sortedErrors)
title('BPM error by EEG channel')
ylabel('BPM error')
xlabel('EEG channels (sorted by BPM error)')