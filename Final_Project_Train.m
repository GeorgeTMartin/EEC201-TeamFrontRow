% EEC 201 Final Project - UC Davis, Winter 2022
% Authors : George Martin & Victor Madrid
% 1. Load Files                 
% 2. Normalize and Cutoff Tails 
% 3. MFCC
%     a. Frame Blocking         
%     b. Windowing 
%     c. FFT
%     d. Log Scale
%     e. Mel Frequency
%     f. Cepstrum
% 4. Cluster & Create Centroids of Training Files
% 5. Refine Codebook
% 6. Testing 
%% 
clc;clear;close all;
FilePath = which('FinalProject_Train');
[FileParts,~,~] = fileparts(FilePath);
TrainingPath = fullfile(FileParts,'Data','Training_Data');
addpath(TrainingPath);
d = dir([TrainingPath, '\*.wav']);
num_data = length(d);
train_set = cell(num_data,2);
%%
for i = 1:num_data
    filename = ['train',num2str(i)];
    train_set{i,1} = join(['train',num2str(i)]);
    train_set{i,2} = join(['Fs',num2str(i)]);
    [train_set{i,1},train_set{i,2}] = audioread([filename,'.wav']);
    train_set{i,1} = (train_set{i,1}-mean(train_set{i,1}))./max(abs(train_set{i,1})); %Normalization to 1
end
for i=1:(num_data) 
    [n,m] = size(train_set{i,1});
    if m == 2
         train_set{i,1}(:,2) = [];
    end
end
%% 
%Generate Ceptsrum Coefficients -------------------------------------------
N=512; %Length of Frame Block  %used to be 256                                                                     %STD 256
M=256; %Shift Between Each  Frame Block  %used to be 100                                                             %STD 64
Window = hamming(N);
%notchSpecs  = fdesign.notch('N,F0,Q,Ap',2,0.01,.5,2);
%notchSpecs = fdesign.bandstop('N,Fp1,Fst1,Fst2,Fp2,C',60,3.6e3,4.8e3,7.4e3,8.6e3,20500)
%notchfilt = design(notchSpecs);
%fvtool(notchfilt)
for i = 1:num_data
    n1 = 1:N-M:length(train_set{i,1})-N;
    n2 = n1 + N;
    train_set_MELFB{i} = melfb(20,N,train_set{i,2});
        for j = 1:1:length(n1)
            %FRAMEBLOCK
            train_set_slice = train_set{i,1};
            train_set_frames{i,j} = train_set_slice(n1(j):(n2(j)-1));
            %WINDOW
            train_set_windowed{i,j} = Window .* train_set_frames{i,j};
            train_set_FFT{i,j} = fft(train_set_windowed{i,j});
           % train_set_FFT{i,j} = filter(notchfilt,train_set_FFT{i,j});            %%COMMENT OUT IF NOT REMOVING SIGNAL BAND
            %MEL FREQ WRAPPING
            N2 = 1 + floor(N/2);
            train_set_MEL{i,j} = train_set_MELFB{i} * abs(train_set_FFT{i,j}(1:N2)).^2;
            %CEPSTRUM
            train_set_CEP{i,j} = dct(log(train_set_MEL{i,j}));
            train_set_CEP{i,j} = train_set_CEP{i,j}(2:end);
        end
end
%% 
%CENTROIDS------------------------------------------------------------------
k=1;
num_entries = size(train_set_CEP);
for i = 1:num_entries(1)
    ts_cep_matrix=[];
    for j = 1:num_entries(2)
        if(isempty(train_set_CEP{i,j}))
            continue
        else
            ts_cep_matrix(:,k) = train_set_CEP{i,j};
            k=k+1;
        end
    end
    codebook{i,1}=mean(ts_cep_matrix,2);
    k=1;
end
%% 
epsilon = .01; %Splitting Factor                                                           %STD .01 [doesnt matter]
while size(codebook,2) < 8                                                                 %STD 16
  for i = 1:num_entries(1) %DOUBLE CODEBOOK SIZE
     for l = 1:size(codebook,2)
     codebook_new{i,2*l-1} = (1+epsilon).*codebook{i,l};
     codebook_new{i,2*l} = (1-epsilon).*codebook{i,l};
     end
  end
  for i = 1:num_entries(1) %OPTIMIZE CURRENT CENTROIDS
            for p = 1:100
             centroid_bin=cell(size(codebook_new,2),1);
             for k = 1:num_entries(2)
                 for m = 1:size(codebook_new,2)
                     if(isempty(train_set_CEP{i,k}))
                        continue
                     else
                         dist_centroid(m) = norm(codebook_new{i,m}-train_set_CEP{i,k},2);
                     end
                 end
                [min_val,min_index] = min(dist_centroid);
                if isempty(centroid_bin{min_index})
                    centroid_bin{min_index} = train_set_CEP{i,k};
                else
                    centroid_bin{min_index} = [centroid_bin{min_index} train_set_CEP{i,k}];
                end
             end
             for m = 1:size(codebook_new,2)
                 if isempty(centroid_bin{m})
                     continue
                 else
                    codebook_new{i,m} = mean(centroid_bin{m},2); %BIG ISSUE
                 end
             end            
           end
  end
 codebook = codebook_new;
end
writecell(codebook,'Codebook')

 for i = 1:size(codebook,2) %CODEBOOK CENTROIDS TRAINSET 1
     d1(i) = codebook{10,i}(1);
     d2(i) = codebook{10,i}(2);
     d3(i) = codebook{5,i}(1);
     d4(i) = codebook{5,i}(2);
 end
 for i = 1:size(train_set_CEP,2)
       if(isempty(train_set_CEP{1,i}))
                    continue
       else
            s1(i) = train_set_CEP{10,i}(1);
            s2(i) = train_set_CEP{10,i}(2);
            s3(i) = train_set_CEP{5,i}(1);
            s4(i) = train_set_CEP{5,i}(2);
       end
 end
 figure
 scatter(d1,d2,'b')
 hold on
 scatter(s1,s2,'g')
 scatter(d3,d4,'black','*')
 scatter(s3,s4,'r','*')
 legend('Codebook 1','Data 1','Codebook 2','Data 2')
 xlabel('Dimension 1')
 ylabel('Dimension 2')
 title('Two-Dimensional Display of Centroid and Clustering Behavior')

%% STFT Stuff 
W = 512;
for i = 1:length(train_set)
    figure()
    stft(train_set{i,1},train_set{i,2},'Window',hamming(W),'OverlapLength',round(W/3),'FFTLength',512)
end

%% STFT Stuff 
W = 512;
for i = 1:length(train_set)
    figure()
    stft(train_set{i,1},train_set{i,2},'Window',hamming(W),'OverlapLength',round(W/3),'FFTLength',512)
end

 MELFBFULL = melfb(20,length(train_set{1,1}),train_set{1,2});
 TRAINSETFFT = fft(train_set{1,1});
 subplot(2,1,1)
 spectrogram(TRAINSETFFT)
 TRAINSETMEL = MELFBFULL * abs(TRAINSETFFT(1:6529)).^2;
 subplot(2,1,2)
 spectrogram(TRAINSETMEL)


 