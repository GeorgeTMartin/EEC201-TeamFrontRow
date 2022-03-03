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
clc;clear;close;
num_data = 11;
train_set = cell(num_data,2);
for i = 1:num_data
    filename = ['train',num2str(i)];
    train_set{i,1} = join(['train',num2str(i)]);
    train_set{i,2} = join(['Fs',num2str(i)]);
    [train_set{i,1},train_set{i,2}] = audioread([filename,'.wav']);
    train_set{i,1} = (train_set{i,1}-mean(train_set{i,1}))./max(abs(train_set{i,1})); %Normalization to 1
end
for i=9:11
    train_set{i,1}(:,2) = [];
end
%Generate Ceptsrum Coefficients -------------------------------------------
N=256; %Length of Frame Block
M=100; %Shift Between Each  Frame Block
Window=.54-.46*cos(2*pi().*[1:N]/(N-1));
for i = 1:num_data
    n1 = 1:N-M:length(train_set{i,1})-N;
    n2 = n1 + N;
    train_set_MELFB{i} = melfb(20,N,train_set{i,2});
        for j = 1:1:length(n1)
            %FRAMEBLOCK
            train_set_slice = train_set{i,1};
            train_set_frames{i,j} = train_set_slice(n1(j):(n2(j)-1));
            %WINDOW
            train_set_windowed{i,j} = Window' .* train_set_frames{i,j};
            train_set_FFT{i,j} = fft(train_set_windowed{i,j});
            %MEL FREQ WRAPPING
            N2 = 1 + floor(N/2);
            train_set_MEL{i,j} = train_set_MELFB{i} * abs(train_set_FFT{i,j}(1:N2)).^2;
            %CEPSTRUM
            train_set_CEP{i,j} = dct(log(train_set_MEL{i,j}));
            train_set_CEP{i,j} = train_set_CEP{i,j}(2:end);
        end
end
%CLUSTER ------------------------------------------------------------------
%Centroids
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
    codebook{i,1}=mean(ts_cep_matrix');
    k=1;
end

epsilon = .001; %Splitting Factor
error = 100;
prev_error = 0; 
l=1;
while abs(error)> 4
        for i = 1:num_entries(1)
            codebook_new{i,2*l-1} = (1+epsilon)*codebook{l};
            codebook_new{i,2*l} = (1-epsilon)*codebook{l};
            k=1;
            for j = 1:num_entries(2)
                if(isempty(train_set_CEP{i,j}))
                    continue
                else
                    dist_centroid_up(j) = norm(codebook_new{i,2*l-1}-train_set_CEP{i,j},2);
                    dist_centroid_down(j) = norm(codebook_new{i,2*l}-train_set_CEP{i,j},2);
                        if dist_centroid_up(j) > dist_centroid_down(j)
                         centroid_up_new(:,k) = train_set_CEP{i,j};
                        else
                         centroid_down_new(:,k) = train_set_CEP{i,j};
                        end
                    k=k+1;
                end
            end
            codebook_new{i,2*l-1}=mean(centroid_up_new');
            codebook_new{i,2*l}=mean(centroid_down_new');
            k=1;
        end
    codebook = codebook_new;
    error = mean(dist_centroid_down') + mean(dist_centroid_up') - prev_error;
    prev_error = mean(dist_centroid_down') + mean(dist_centroid_up');
    l=l+1;
end

% %PLOT CLUSTER PROOF
% for i = 1:size(codebook,1)
%     d1(i) = codebook(1,)
% end