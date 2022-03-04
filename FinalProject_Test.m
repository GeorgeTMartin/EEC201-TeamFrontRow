% EEC 201 Final Project - UC Davis, Winter 2022
% Authors : George Martin & Victor Madrid
clc;clearvars -EXCEPT codebook;close;
addpath(genpath(pwd))
num_data = 11;
train_set = cell(num_data,2);
for i = 1:num_data
    filename = ['s',num2str(i)];
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
%Window=.54-.46*cos(2*pi().*[1:N]/(N-1));
Window = hamming(N);
for i = 1:num_data
    n1 = 1:N-M:length(train_set{i,1})-N;
    n2 = n1 + N;
    test_set_MELFB{i} = melfb(20,N,train_set{i,2});
        for j = 1:1:length(n1)
            %FRAMEBLOCK
            test_set_slice = train_set{i,1};
            test_set_frames{i,j} = test_set_slice(n1(j):(n2(j)-1));
            %WINDOW
            test_set_windowed{i,j} = Window .* test_set_frames{i,j};
            test_set_FFT{i,j} = fft(test_set_windowed{i,j});
            %MEL FREQ WRAPPING
            N2 = 1 + floor(N/2);
            test_set_MEL{i,j} = test_set_MELFB{i} * abs(test_set_FFT{i,j}(1:N2)).^2;
            %CEPSTRUM
            test_set_CEP{i,j} = dct(log(test_set_MEL{i,j}));
            test_set_CEP{i,j} = test_set_CEP{i,j}(2:end);
        end
end
%Total_Distortion = zeros(size(codebook,1),size(codebook,2));
for i = 1:size(test_set_CEP,1)
    %Distortion_sum = [];
    Nearest_Ref = zeros(1,size(test_set_CEP,1));
    for j = 1:size(test_set_CEP,2)
        for k = 1:size(codebook,1)
            for m = 1:size(codebook,2)
                if isempty(test_set_CEP{i,j})
                   continue
                else
                   centroid_dist(k,m) = norm(codebook{k,m}(:)-test_set_CEP{i,j}(:),2);
                end
            end
        end
%         Total_Distortion = Total_Distortion + centroid_dist;
          [minval, minidx] = min(centroid_dist(:));
          [p,q] = ind2sub(size(centroid_dist),minidx);
          Nearest_Ref(p) = Nearest_Ref(p)+1;
    end
[temp, out_ind] = max(Nearest_Ref);
out(i) = out_ind;
% Distortion_sum = sum(Total_Distortion,2);
% [mini, ind] = min(Distortion_sum);
% out(i)=ind;

%  ATTEMPT 1
%  [minval, minidx] = min(Total_Distortion(:));
%  [p,q] = ind2sub(size(Total_Distortion),minidx);
%  out(i)=p;
end
out


% min_count_bin = zeros(1,size(test_set_CEP,1));
% for i = 1:size(test_set_CEP,1)
%     for j = 1:size(test_set_CEP,2)
%         for k = 1:size(codebook,1)
%           for m = 1:size(codebook,2)
%               if isempty(test_set_CEP{i,j})
%                   continue
%               else
%                   centroid_dist(k,m) = norm(codebook{k,m}(:)-test_set_CEP{i,j}(:),2);
%               end
%           end
%         end 
%       [minval, minidx] = min(centroid_dist(:));
%       [p,q] = ind2sub(size(centroid_dist),minidx);
%        min_count_bin(p) = min_count_bin(p) + 1;
%     end
%     [max_count max_ind] = max(min_count_bin);
%     output(i)= max_ind;
% end
% 
% output'
