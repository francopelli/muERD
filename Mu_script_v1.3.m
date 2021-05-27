%fran copelli march 2020
%/Volumes/UNTITLED/paper_ripping/0_raw

% Pipeline
% ------------------------------------
% EEG lab 19.1
% MATLAB R2019b

% 1. Create EEG structure
% 2. label/remove unused channels (EXG7&8)
% 3. add channel locations reference (note: this will change depending on matlab folder location)
% 4. filter (high-pass: 1)
% 5. filter (low-pass: 60hz)
% 6. reject bad channels
% 7. interpolate removed channels
% 8. rereference (to average)
% 9. epoch, remove baseline
% 10. ICA (runica)
% 11. delete noisy epochs
% 12. run ICA again
% 13. DIPFIT3.3
% 14. IC label for IC rejection
% 15. create study file
% 16. study design and preselect components for clustering
% 17. compute component measures (PCA)
% 18. compute and visualize clusters (kmeans). on average, choose one component per participant per cluster.
% 19. export data
%%%NEED TO UPDATE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PART ONE
close all;
clear all;
eeglab;

for subject = [1001:1021]
    sub=int2str(subject);
    file1=[sub,'.bdf'];
    inpath=['/Volumes/UNTITLED/paper_ripping/2020/0_raw/'];
    input1=[inpath,file1];
    outpath=['/Volumes/UNTITLED/paper_ripping/2020/1_preprocess/'];
    output=[sub,'_processed.set'];
    splitpath=['/Volumes/UNTITLED/paper_ripping/2020/2_analysis_source/'];
    checkfile=exist(inpath);
    
    %insert iff statement forr other EEG types
    %user input prompt forr how many channels total, incorporate into line 48
        EEG = pop_biosig(input1,'ref',[65,66]); 
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
        EEG = eeg_checkset( EEG );
        [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); 
        EEG = eeg_checkset( EEG );
        
        %names the file
        EEG.setname=sub;
        EEG = eeg_checkset( EEG );
        EEG.subject = sub;
        
        %remove unused channels
        %dependent on EEG system type
        EEG = pop_select( EEG,'nochannel',{'EXG7' 'EXG8'}); 
        
        % label EXG channels
        EEG = pop_chanedit(EEG,'changefield',{65 'labels' 'M1' });
        EEG = pop_chanedit(EEG,'changefield',{66 'labels' 'M2' });
        EEG = pop_chanedit(EEG,'changefield',{67 'labels' 'LO1'});
        EEG = pop_chanedit(EEG,'changefield',{68 'labels' 'LO2'});
        EEG = pop_chanedit(EEG,'changefield',{69 'labels' 'IO1'});
        EEG = pop_chanedit(EEG,'changefield',{70 'labels' 'IO2'});

        %save
        EEG = pop_saveset( EEG,'filename',[sub,'_Aremoveexg78.set'],'filepath',outpath);
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        
        %channel locations
        %user prompt to select EEGLAB directory needed
        EEG=pop_chanedit(EEG, 'lookup','/Users/smartlab/Documents/MATLAB/eeglab2019_1/plugins/dipfit-3.3/standard_BEM/elec/standard_1005.elc');
        EEG = eeg_checkset( EEG );
        
        %highpass filter
        EEG = pop_eegfiltnew(EEG, 'locutoff',1);
        
        %lowpass filter
        EEG = pop_eegfiltnew(EEG, 'hicutoff',50);
        
        %save
        EEG = pop_saveset( EEG,'filename',[sub,'_Bfilter.set'],'filepath',outpath);
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        
        %interpolate noisy scalp channels
        %depends on EEG type
        EEGtemp = EEG;
        [EEGtemp, badchans] = pop_rejchan(EEG, 'elec',[1:64],'threshold',5,'norm','on','measure','kurt');
        EEG = pop_interp(EEG, badchans);
       
        %reject noisy EXG channels
        [EEG, badEXGchans] = pop_rejchan(EEG, 'elec',[65:70] ,'threshold',5,'norm','on','measure','kurt');
        %noisy scalp channels are interpolated. noisy EXg channels are
        %rejected.
        
        %save
        EEG = pop_saveset( EEG,'filename',[sub,'_Cinterp.set'],'filepath',outpath);
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        
        %save interpolated channel list
        badEXGchans = badEXGchans + 64;
        chanrej = numel([badchans badEXGchans]);
        chanrejlist = [chanrej,badchans,badEXGchans];
        fileID = fopen([sub,'_interpolatedchannels.txt'],'w');
        fprintf(fileID, '%i\n', chanrejlist);
        fclose(fileID);
        %note that the first value of the txt file is the total amount of
        %rejected or interpolated channels. if channels 1-64 appear in the list, they would have been
        %interpolated, and if channels 65-70 appear in the list, they would have been rejected.
     
        %rereference to average
        EEG = pop_reref( EEG, [],'exclude',[65:68]);
        
        %save
        EEG = pop_saveset( EEG,'filename',[sub,'_Drereftoavg.set'],'filepath',outpath);
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        
        %epoch
        %user prompt forr epoch length
        EEG = pop_epoch( EEG, {}, [-1.5  8], 'newname', 'BDF file epochs', 'epochinfo', 'yes');
        EEG = eeg_checkset( EEG );
        
        %remove baseline
        EEG = pop_rmbase( EEG, [-1500 0]);
        EEG = eeg_checkset( EEG );
        
        %save
        EEG = pop_saveset( EEG,'filename',[sub,'_Eepoch.set'],'filepath',outpath);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        
        %%decomp with ICA
        EEG = eeg_checkset( EEG );
        ICs = 69-chanrej; % PCA dimension must be n-1, then reduced to exclude bad/interpolated channels. this will be your resulting ICs
        EEG = pop_runica(EEG, 'extended',1,'interupt','on','pca',ICs);
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        
        %save
        EEG = pop_saveset( EEG,'filename',[sub,'_FICA.set'],'filepath',outpath);
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        
        %automatic epoch rejection
        [EEG, rmepochs] = pop_autorej(EEG, 'nogui','on','eegplot','off',‘icacomps’, [1:size(EEG.icaact,1)]);
        EEG = pop_rejepoch( EEG, rmepochs ,0);
        
        %save
        EEG = pop_saveset( EEG,'filename',[sub,'_Gepochrej.set'],'filepath',outpath);
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        
        %save list of rejected epochs
        epochrej = numel(rmepochs);
        epochrejlist = [epochrej,rmepochs];
        fileID = fopen([sub,'_rejectedepochs.txt'],'w');
        fprintf(fileID, '%i\n', epochrejlist);
        fclose(fileID);
        %note that the first value of the txt file is the total amount of
        %rejected epochs. all the following values are the individual
        %epoch numbers.
        
        %%decomp with ICA  round 2
        EEG = eeg_checkset( EEG );
        EEG = pop_runica(EEG, 'extended',1,'interupt','on','pca',ICs);
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        
        %save
        EEG = pop_saveset( EEG,'filename',[sub,'_HICA2.set'],'filepath',outpath);
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        
        % dipfit single dipole
        EEG = pop_dipfit_settings( EEG, 'hdmfile','/Users/smartlab/Documents/MATLAB/eeglab2019_1/plugins/dipfit-3.3/standard_BEM/standard_vol.mat','coordformat','MNI','mrifile','/Users/smartlab/Documents/MATLAB/eeglab2019_1/plugins/dipfit-3.3/standard_BEM/standard_mri.mat','chanfile','/Users/smartlab/Documents/MATLAB/eeglab2019_1/plugins/dipfit-3.3/standard_BEM/elec/standard_1005.elc','chansel',[1:64] );
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);
        EEG = pop_multifit(EEG, [1:ICs] ,'threshold',100,'rmout','on');
        %fit all available independent components
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        
        %save
        EEG = pop_saveset( EEG,'filename',[sub,'_Idipfit.set'],'filepath',outpath);
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        
        %IC labeling
        EEG = pop_iclabel(EEG, 'default');
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        
        %muscle reject values
        muscle = EEG.etc.ic_classification.ICLabel.classifications(:,2);
        musclereject = find(muscle >= 0.9);
        musclenum = numel(musclereject);
        ICrejlist = [musclenum,musclereject'];
        fileID = fopen([sub,'_musclereject.txt'],'w');
        fprintf(fileID, '%i\n', ICrejlist);
        fclose(fileID);
        %note that the first value of the txt file is the total amount of
        %rejected ICs. all the following values are the individual
        %IC numbers.
        
        eye = EEG.etc.ic_classification.ICLabel.classifications(:,3);
        eyereject = find(eye >= 0.9);
        eyenum = numel(eyereject);
        ICrejlist = [eyenum,eyereject'];
        fileID = fopen([sub,'_eyereject.txt'],'w');
        fprintf(fileID, '%i\n', ICrejlist);
        fclose(fileID);
        
        heart = EEG.etc.ic_classification.ICLabel.classifications(:,4);
        heartreject = find(heart >= 0.9);
        heartnum = numel(heartreject);
        ICrejlist = [heartnum,heartreject'];
        fileID = fopen([sub,'_heartreject.txt'],'w');
        fprintf(fileID, '%i\n', ICrejlist);
        fclose(fileID);
        
        line = EEG.etc.ic_classification.ICLabel.classifications(:,5);
        linereject = find(line >= 0.9);
        linenum = numel(linereject);
        ICrejlist = [linenum,linereject'];
        fileID = fopen([sub,'_linereject.txt'],'w');
        fprintf(fileID, '%i\n', ICrejlist);
        fclose(fileID);
        
        channel = EEG.etc.ic_classification.ICLabel.classifications(:,6);
        channelreject = find(channel >= 0.9);
        channelnum = numel(channelreject);
        ICrejlist = [channelnum,channelreject'];
        fileID = fopen([sub,'_channelreject.txt'],'w');
        fprintf(fileID, '%i\n', ICrejlist);
        fclose(fileID);
        
        other = EEG.etc.ic_classification.ICLabel.classifications(:,7);
        otherreject = find(other >= 0.9);
        othernum = numel(otherreject);
        ICrejlist = [othernum,otherreject'];
        fileID = fopen([sub,'_otherreject.txt'],'w');
        fprintf(fileID, '%i\n', ICrejlist);
        fclose(fileID);
        
        %IC rejection
        [EEG] = pop_icflag(EEG, [NaN NaN;0.9 1;0.9 1;0.9 1;0.9 1;0.9 1;0.9 1]);
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        
        %save
        EEG = pop_saveset( EEG,'filename',[sub,'_JICrejection.set'],'filepath',outpath);
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
        
        %relabel epochs
        
        EEG = pop_selectevent( EEG, 'type',[11:10:41] ,'renametype','A','deleteevents','off','deleteepochs','off','invertepochs','off');
        EEG = pop_selectevent( EEG, 'type',[12:10:42] ,'renametype','CA','deleteevents','off','deleteepochs','off','invertepochs','off');
        EEG = pop_selectevent( EEG, 'type',[13:10:43] ,'renametype','CV','deleteevents','off','deleteepochs','off','invertepochs','off');
        EEG = pop_selectevent( EEG, 'type',[14:10:44] ,'renametype','V','deleteevents','off','deleteepochs','off','invertepochs','off');
        EEG = pop_selectevent( EEG, 'type',[15:10:45] ,'renametype','CAV','deleteevents','off','deleteepochs','off','invertepochs','off');
        EEG = pop_selectevent( EEG, 'type',[16:10:46] ,'renametype','AV','deleteevents','off','deleteepochs','off','invertepochs','off');
        EEG = pop_selectevent( EEG, 'type',[17:10:47] ,'renametype','ACV','deleteevents','off','deleteepochs','off','invertepochs','off');
        EEG = pop_selectevent( EEG, 'type',[18:10:48] ,'renametype','VCA','deleteevents','off','deleteepochs','off','invertepochs','off');
        EEG = pop_selectevent( EEG, 'type',[19:10:49] ,'renametype','ACT','deleteevents','off','deleteepochs','off','invertepochs','off');

        EEG = pop_saveset( EEG,'filename',[sub,'_Kprocessed.set'],'filepath',outpath)
        
close all;
clear all;
eeglab;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PART TWO
%create study file

close all;
clear all;
eeglab;

%create study file
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
[STUDY ALLEEG] = std_editset( STUDY, [], 'name','Paper Ripping Experiment','commands',{{'index' 1 'load' '/Volumes/UNTITLED/paper_ripping/2020/1_preprocess/1001_Kprocessed.set'} {'index' 2 'load' '/Volumes/UNTITLED/paper_ripping/2020/1_preprocess/1002_Kprocessed.set'} {'index' 3 'load' '/Volumes/UNTITLED/paper_ripping/2020/1_preprocess/1003_Kprocessed.set'} {'index' 4 'load' '/Volumes/UNTITLED/paper_ripping/2020/1_preprocess/1004_Kprocessed.set'} {'index' 5 'load' '/Volumes/UNTITLED/paper_ripping/2020/1_preprocess/1005_Kprocessed.set'} {'index' 6 'load' '/Volumes/UNTITLED/paper_ripping/2020/1_preprocess/1006_Kprocessed.set'} {'index' 7 'load' '/Volumes/UNTITLED/paper_ripping/2020/1_preprocess/1007_Kprocessed.set'} {'index' 8 'load' '/Volumes/UNTITLED/paper_ripping/2020/1_preprocess/1008_Kprocessed.set'} {'index' 9 'load' '/Volumes/UNTITLED/paper_ripping/2020/1_preprocess/1009_Kprocessed.set'} {'index' 10 'load' '/Volumes/UNTITLED/paper_ripping/2020/1_preprocess/1010_Kprocessed.set'} {'index' 11 'load' '/Volumes/UNTITLED/paper_ripping/2020/1_preprocess/1011_Kprocessed.set'} {'index' 12 'load' '/Volumes/UNTITLED/paper_ripping/2020/1_preprocess/1012_Kprocessed.set'} {'index' 13 'load' '/Volumes/UNTITLED/paper_ripping/2020/1_preprocess/1013_Kprocessed.set'} {'index' 14 'load' '/Volumes/UNTITLED/paper_ripping/2020/1_preprocess/1014_Kprocessed.set'} {'index' 15 'load' '/Volumes/UNTITLED/paper_ripping/2020/1_preprocess/1015_Kprocessed.set'} {'index' 16 'load' '/Volumes/UNTITLED/paper_ripping/2020/1_preprocess/1016_Kprocessed.set'} {'index' 17 'load' '/Volumes/UNTITLED/paper_ripping/2020/1_preprocess/1017_Kprocessed.set'} {'index' 18 'load' '/Volumes/UNTITLED/paper_ripping/2020/1_preprocess/1018_Kprocessed.set'} {'index' 19 'load' '/Volumes/UNTITLED/paper_ripping/2020/1_preprocess/1019_Kprocessed.set'} {'index' 20 'load' '/Volumes/UNTITLED/paper_ripping/2020/1_preprocess/1020_Kprocessed.set'} {'index' 21 'load' '/Volumes/UNTITLED/paper_ripping/2020/1_preprocess/1021_Kprocessed.set'} {'inbrain' 'on' 'dipselect' 0.15}},'updatedat','on','rmclust','on' );
[STUDY ALLEEG] = std_checkset(STUDY, ALLEEG);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
%edit study design
STUDY = std_makedesign(STUDY, ALLEEG, 1, 'name','All conditions','delfiles','off','defaultdesign','off','variable1','type','values1',{'A' 'ACT' 'ACV' 'AV' 'CA' 'CAV' 'CV' 'V' 'VCA'},'vartype1','categorical','subjselect',{'1001' '1002' '1003' '1004' '1005' '1006' '1007' '1008' '1009' '1010' '1011' '1012' '1013' '1014' '1015' '1016' '1017' '1018' '1019' '1020' '1021'});
[STUDY EEG] = pop_savestudy( STUDY, EEG, 'filename','Paper_ripping_2020.study','filepath','/Volumes/UNTITLED/paper_ripping/2020/1_preprocess/');
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
%precompute component measures
[STUDY ALLEEG] = std_precomp(STUDY, ALLEEG, 'components','savetrials','on','recompute','on','scalp','on','spec','on','specparams',{'specmode' 'fft' 'logtrials' 'off'},'ersp','on','erspparams',{'cycles' [3 0.8]  'nfreqs' 100 'ntimesout' 200},'itc','on');
eeglab('redraw');
%preclustering array
[STUDY ALLEEG] = std_preclust(STUDY, ALLEEG, 1,{'dipoles' 'weight' 1},{'ersp' 'npca' 10 'freqrange' [8 13]  'timewindow' [] 'weight' 1});
%cluster
[STUDY] = pop_clust(STUDY, ALLEEG, 'algorithm','kmeans','clus_num',  round(length(STUDY.cluster.comps)/length(STUDY.subject)) , 'outliers',  3 );
%number of clusters is the average amount of components per subject

%save
[STUDY EEG] = pop_savestudy( STUDY, EEG, 'filename','Paper_ripping_2020.study','filepath','/Volumes/UNTITLED/paper_ripping/2020/1_preprocess/');
%ready to visualize!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PART THREE
%joseph rovetti march 2020
[STUDY erspdata ersptimes erspfreqs pgroup pcond pinter] = std_erspplot(STUDY,ALLEEG,'clusters', 10) % generate ersp for given cluster (specify cluster number)
% generate ersp for given cluster (specify cluster number)
[STUDY specdata specfreqs] = std_specplot(STUDY,ALLEEG,'clusters', 10) % generate cluster spectra for given cluster (specify cluster number)  

cd '/Users/smartlab/Documents/paper_ripping/2020/' % sets the directory of where to export the data at the end
ersp_collapsed = zeros(size(erspdata{1},3),size(erspdata,1)); % makes an empty matrix to put the data in
for aa = 1:length(erspdata) % loops through each of the nine conditions
    data = erspdata{aa}; % isolates the conditions one at a time
    data = data(:,24:end,:); % takes just the post-zero (i.e., non-baseline) data
    data = mean(data,[1,2]); % averages across frequency and time
    data = permute(data,[1,3,2]); % turns the data, still 3D, into a 1D column
    ersp_collapsed(:,aa) = data; % puts that column of data into the the data matrix we made before
end
ersp_means = mean(ersp_collapsed); % finds the mean of each condition (across all components)
table = array2table(ersp_collapsed,... % puts the data into a table to be exported
    'VariableNames',{'A','ACT','ACV','AV','CA','CAV','CV','V','VCA'}); % labels all conditions
writetable(table,'ersp_collapsed.xlsx') % exports the data table to the directory chosen before