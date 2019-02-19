%path must be in brain_data folder (where your raw data folders are)
%currently written to handle dyads
%required inpaint_nans, homer2 scripts, and huppertt (nirs-toolbox) scripts

%modified by MCW for Couples study; remodified for Jordan

%%%%%%%%%%%%%%%%%%%%%% USER INPUTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
raw_dir='Jordan_fNIRS_Data/To_Preprocess';
dataprefix='JOR';
probeInfoFile=0;
dyads=0;
Jordan=1; %change to 0 if using for other projects

%dosubs=[];

%inputs: raw_dir: folder path with all single-subject or dyad-level raw data
%                 files (path must be relative to current matlab path).
%       dataprefix: string. Prefix of every folder name that should be considered a
%       data folder. E.g., MIN for MIN_101, MIN_102, etc. 
%       probeInfoFile: 0 or 1. If 0, will look for it in the first data folder.
%                 If 1, will ask you to provide a probeInfoFile before
%                 running.  
%       dyads: 0 or 1. 1 if hyperscanning, 0 if single subject.
%       dosubs: if all, leave at []. if you want to do a specific subset,
%       specify based on folder order, where 1 refers to first folder in data, etc.
%               -- e.g. [1:5, 9:15]
%
%outputs: preprocessed and .nirs files in same directory level as raw_dir

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DEBUGGING TIPS:
%Note that this function is, well, funcitonal, but not 100% optimized for
%everything that could go wrong. A couple problems I have noticed so far,
%and how to fix them if it comes up for you:
%   - If you get an "Index exceeds matrix dimensions" error in
%   hmrMotionArtifact for a subject that's not the first file: 
%       Check the SD Mask structure in the .hdr of that subject to see if 
%       it matches the SD Mask structure of subject 1. If the wrong montage
%       was selected in recording, will cause this error. Simply copy-paste
%       the correct SD Mask and ChannelDistance list into the .hdr file.
%   - Delete all false start files from the data directory, or will cause
%   script to error out. 

addpath(genpath('preprocessing_tools'));
addpath(genpath('functions'));

currdir=dir(strcat(raw_dir,'/',dataprefix,'*'));
%if ~isempty(dosubs)
 %   currdir=currdir(dosubs);
%end
       
probecheck = 0;
if probeInfoFile
    [probefile,probepath] = uigetfile('*_probeInfo.mat');
    load(fullfile(probepath,probefile));
    probecheck=1;
end  


    %all again but no dyad stuff
    fprintf('\n\t Preprocessing ...\n')
    reverseStr = '';
    for i=1:length(currdir)
        if Jordan
            life_duration=324;
            choice_duration=281;
        end

        subj=currdir(i).name;
        subjdir=dir(strcat(raw_dir,'/',subj,'/',dataprefix,'*'));
        %subjdir=dir(strcat(raw_dir,'/',subj);
        msg = sprintf('\n\t subject number %d/%d ...',i,length(currdir));
        fprintf([reverseStr,msg]);
        reverseStr = repmat(sprintf('\b'),1,length(msg));
        for j=1:length(subjdir)
            event_ind=[];
            no_events=99;
            scanname = subjdir(j).name;
            subjfolder = strcat(raw_dir,'/',subj,'/',scanname);
            sprintf('\n\t doing scan %s',scanname)
            while ~probecheck
                probefile = dir(strcat(subjfolder,'/','*_probeInfo.mat'));
                probefilename = probefile(1).name;
                load(strcat(subjfolder,'/',probefilename));
                probecheck=1;
            end
            if ~probecheck
                error('ERROR: Cannot find probeInfo file in first subject folder');
            end
            
            outpath = strcat('PreProcessedFiles_cut/',subj,'/',scanname);
            if ~exist(outpath,'dir')
                %1) extract data values
                clear d sd_ind samprate wavelength s cutd;
                [d, sd_ind, samprate, wavelengths, s] = extractData(subjfolder);
                
                %This is for Jordan video files. It allows us to cut the
                %files (all to the same length, being the length of the video)
                %before they get preprocessed
                
                %if Jordan
                    if contains(scanname,'choice')
                        choice_duration_full=choice_duration*samprate;
                        event_ind=find(s==1); %where was the trigger?
                        if size(event_ind,1)>1 % if more than one occurred, use the last one
                            event_ind=event_ind(length(event_ind));
                        end
                        if strcmp(subj,'JOR_038')==1
                            event_ind=find(s==1);
                            event_ind=event_ind(1);
                        end
                        if size(event_ind,1)==0 % take note if no events
                            no_events=1;
                            d=d;
                        else
                            cutd=d(event_ind:(event_ind+(round(choice_duration_full)-1)),:);
                            d=cutd;
                            no_events=0;
                        end
                    
                    elseif contains(scanname,'life')
                        life_duration_full=life_duration*samprate;
                        event_ind=find(s==1); %where was the trigger?
                        if size(event_ind,1)>1 % if more than one occurred, use the last one
                            event_ind=event_ind(length(event_ind));
                        end
                        
                        if isempty(event_ind)
                            no_events=1;
                            %d=d;
                        else
                            cutd=d(event_ind:(event_ind+(round(life_duration_full)-1)),:);
                            no_events=0;
                            d=cutd;
                        end
                    end
               % end
    
                %2) identify and remove bad channels
                %bad channel defined as any where detector saturation occurs for >2sec, 
                %or where power spectrum variation is too high. 
                %Feel free to change these parameters if you have a good reason to do so
                %
                %reasoning for default choices:
                %- if saturation occurs, data will be 'NaN'. But if this only lasts a
                %short amount of time (e.g. <8 points=<2 seconds at 4Hz), we can fill in what 
                %those data points would have likely been with reasonable confidence.
                %
                %- power spectrum of the signal shows how many sine waves at each freq.
                %make up the raw signal. Good signal should have a large peak at lower
                %frequencies. Pure noise will have random numbers of all freqencies. 
                %We will use a modified version of the quartile coefficient of
                %dispersion
                %(https://en.wikipedia.org/wiki/Quartile_coefficient_of_dispersion)
                %to automatically decide which channels have good or bad
                %signal. Essentially, it sums the frequency amplitudes in the
                %first and third quartiles of the frequency range, and then
                %compares them via (Q1-Q3)/(Q1+Q3). Larger QCoD is cleaner
                %signal. Default threshold is set to 0.1. Change this to <0.1 
                %to allow for greater noise in the signal, or change to >0.1 
                %for more stringency. 
    
                satlength = 2; %in seconds
                QCoDthresh = 0.1;
                [channelmask, d] = removeBadChannels(d, samprate, satlength, QCoDthresh);
    
                %3) convert to .nirs format
                [SD, aux, t] = getRemainingNirsVars(d, sd_ind, samprate, wavelengths, probeInfo, channelmask);
            
                %4) motion filter, convert to hemodynamic changes
                d=d;
                SD=SD;
                
                 %see hmrMotionArtifact in Homer2 documentation for parameter description
    numchannels = size(d,2)/2;
    tInc = hmrMotionArtifact(d, samprate, SD, ones(length(d)), 0.5, 2, 10, 5); %gives us a list of all timepoints with 1 for data looked good across channels and 0 if there seemed to be a major motion artifact across all channels
    %see hmrMotionCorrectPCA in Homer2 documentation for parameter description
    
    % filter the good channels that are left
    
    mlAct = SD.MeasListAct; %what are the good channels?
    
    
%tInc = procResult.tIncAuto;        % identify motion (vector of 1-no motion, and 0-motion)

%dod = procResult.dod;  % delta OD
    nSV = .9;

    lstNoInc = find(tInc==0); %find wherever there was a motion artifact across channels
    lstAct = find(mlAct==1); %index the good channels

    if isempty(lstNoInc)
        dN = d;
        svs = [];
        nSV = 0;
        %return;
    else

    %
    % Do PCA
    %
    y = d(lstNoInc,lstAct); % find where the "good" channels have a bit of bad motion
    yc = y;
    yo = y;

    c = y.' * y;
    [V,St,foo] = svd(c);
    svs = diag(St) / sum(diag(St));

    svsc = svs;
    for idx = 2:size(svs,1)
        svsc(idx) = svsc(idx-1) + svs(idx);
    end
    if nSV<1 & nSV>0 % find number of SV to get variance up to nSV
        ev = diag(svsc<nSV);
        nSV = find(diag(ev)==0,1)-1;
    end

    ev = zeros(size(svs,1),1);
    ev(1:nSV) = 1;
    ev = diag(ev);

    yc = yo - y*V*ev*V';


    %
    % splice the segments of data together
    %
    lstMs = find(diff(tInc)==-1);%find starts 
    lstMf = find(diff(tInc)==1);% and ends of bad sections
    if isempty(lstMf) 
        lstMf = length(tInc);
    end
    if isempty(lstMs)
        lstMs = 1;
    end
    if lstMs(1)>lstMf(1)
        lstMs = [1;lstMs];
    end
    if lstMs(end)>lstMf(end)
        lstMf(end+1,1) = length(tInc);
    end
    lstMb = lstMf-lstMs;
    for ii=2:length(lstMb)
        lstMb(ii) = lstMb(ii-1) + lstMb(ii);
    end

    dN = d;
    
    for ii=1:length(lstAct)

        jj = lstAct(ii);

        lst = (lstMs(1)):(lstMf(1)-1);
        if lstMs(1)>1
            dN(lst,jj) = yc(1:lstMb(1),ii) - yc(1,ii) + dN(lst(1),jj);
        else
            dN(lst,jj) = yc(1:lstMb(1),ii) - yc(lstMb(1),ii) + dN(lst(end),jj);
        end

        for kk=1:(length(lstMf)-1)
            lst = (lstMf(kk)-1):lstMs(kk+1);
            dN(lst,jj) = d(lst,jj) - d(lst(1),jj) + dN(lst(1),jj);

            lst = (lstMs(kk+1)):(lstMf(kk+1)-1);
            dN(lst,jj) = yc((lstMb(kk)+1):lstMb(kk+1),ii) - yc(lstMb(kk)+1,ii) + dN(lst(1),jj);
        end

        if lstMf(end)<length(d)
            lst = (lstMf(end)-1):length(d);
            dN(lst,jj) = d(lst,jj) - d(lst(1),jj) + dN(lst(1),jj);        
        end

    end
    end


    
    
    dfiltered=dN;
    %[dfiltered,~,~] = hmrMotionCorrectPCA(SD, d, tInc, 2); 
    %see hmrIntensity2Conc in Homer2 documentation for parameter description
    [dconverted, ~] = hmrIntensity2Conc(dfiltered, SD, samprate, 0.005, 0.5, [6, 6]);
    dnormed = zscore(dconverted);
    dnormed(:,:,(SD.MeasListAct')==0)=NaN;
    dconverted(:,:,(SD.MeasListAct')==0)=NaN;
    oxy1 = zeros(length(dconverted), numchannels);
    deoxy1 = zeros(length(dconverted), numchannels);
    totaloxy1 = zeros(length(dconverted), numchannels);
    z_oxy1 = zeros(length(dnormed), numchannels);
    z_deoxy1 = zeros(length(dnormed), numchannels);
    z_totaloxy1 = zeros(length(dnormed), numchannels);
    for c = 1:numchannels
        oxy1(:,c) = dconverted(:,1,c);
        deoxy1(:,c) = dconverted(:,2,c);
        totaloxy1(:,c) = dconverted(:,3,c);
        z_oxy1(:,c) = dnormed(:,1,c);
        z_deoxy1(:,c) = dnormed(:,2,c);
        z_totaloxy1(:,c) = dnormed(:,3,c);
    end  
                
                
                if no_events==1
                    outpath=strcat(outpath,'_NO_EVENTS');
                end
                mkdir(outpath) 
                save(strcat(outpath,'/',scanname,'_preprocessed.mat'),'oxy1', 'deoxy1', 'totaloxy1','z_oxy1', 'z_deoxy1', 'z_totaloxy1');
                save(strcat(outpath,'/',scanname,'.nirs'),'aux','d','s','SD','t');
            end
        end
    end

