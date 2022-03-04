%% Grab data
num_data = 11;
test_set = cell(num_data,2);
for i = 1:num_data
    filename = ['s',num2str(i)];
    test_set{i,1} = join(['s',num2str(i)]);
    test_set{i,2} = join(['Fs',num2str(i)]);
    [test_set{i,1},test_set{i,2}] = audioread([filename,'.wav']);
    test_set{i,1} = (test_set{i,1}-mean(test_set{i,1}))./max(abs(test_set{i,1})); %Normalization to 1
end
%% Remove stero
for i=9:11
    test_set{i,1}(:,2) = [];
end

%% 
codebook = readmatrix("Codebook_final.txt");

%% %% Generate Ceptsrum Coefficients 
N=256; %Length of Frame Block
M=100; %Shift Between Each  Frame Block
Window=.54-.46*cos(2*pi().*[1:N]/(N-1));
for i = 1:num_data
    n1 = 1:N-M:length(test_set{i,1})-N;
    n2 = n1 + N;
    test_set_MELFB{i} = melfb(20,N,test_set{i,2});
        for j = 1:1:length(n1)
            %FRAMEBLOCK
            test_set_slice = test_set{i,1};
            test_set_frames{i,j} = test_set_slice(n1(j):(n2(j)-1));
            %WINDOW
            test_set_windowed{i,j} = Window' .* test_set_frames{i,j};
            test_set_FFT{i,j} = fft(test_set_windowed{i,j});
            %MEL FREQ WRAPPING
            N2 = 1 + floor(N/2);
            test_set_MEL{i,j} = test_set_MELFB{i} * abs(test_set_FFT{i,j}(1:N2)).^2;
            %CEPSTRUM
            test_set_CEP{i,j} = dct(log(test_set_MEL{i,j}));
            guest_speaker{i,j} = test_set_CEP{i,j}(2:end);
        end
end
%% format the code book 
for i=1:num_data
    p2 = 19;
    p1 =1;
    for j=1:((length(codebook)/19))
        code_cell{i,j} = codebook(i,p1:p2)';
        p1 = p2+1;
        p2 = p1 + 18;
    end
end
%% Find how many elements are in each row
size_ = (cellfun('size',guest_speaker,1));
length_row =[];
for i=1:num_data
    col = nnz(size_(i,:));
    length_row = [length_row,col];
end

%% get the euclidean and closest centroid
% total_centr = [];
% per_speaker_dist = [];
% euc_dist = [];
% euc_dist2 = [];
% temp = [];
% for i=1:num_data
%     for j=1:(length_row(i))
%         temp =[];
%         voice_num = [];
%         speaker_num = [];
%         for k =1:num_data
%             reg = [];
%             for L=1:8 % being the number of centroids
%                 euc_dist(L) = disteu(code_cell{k,L},test_set_CEP{i,j});
%                 reg = [reg,euc_dist(L)];% 8 centroids - 1 point
%             end
%             [temp,voice] = min(reg);
%             voice_num = [voice_num,voice];
%             reg2=[];
%             for q = 1:length(voice_num)
%                 euc_dist2(L) = disteu(code_cell{q,voice_num(q)},test_set_CEP{i,j});
%                 reg2 = [reg2,euc_dist(L)];
%             end
%             [temp, speaker] = min(reg2);
%             speaker_num = [speaker_num, speaker];
%         end
%         
%         
%     end
% end

%% find out the length of the code book
size_ = (cellfun('size',code_cell,1));
len_code_cell =[];
for i=1:num_data
    col = nnz(size_(i,:));
    len_code_cell = [len_code_cell,col];
end
len_code_cell = len_code_cell(1);
%% get the euc distance and match speakers
euc_dist = [];
reg =[];
centr_avg =[];
guest_avg =[];
index_guess = [];
for i=1:num_data%number of known speakers
    guest_avg =[];
    for F=1:11%number of guest speakers
        centr_avg =[];
        for L=1:len_code_cell%number of centroids
            reg =[];
            for j=1:(length_row(F))%number of data points per guest speaker
                euc_dist = disteu(code_cell{i,L}, guest_speaker{F,j});
                reg = [reg,euc_dist];
            end
            prev_avg = min(reg); %USED TO BE MEAN NOW ITS MIN
            centr_avg = [centr_avg,prev_avg];
        end
        prev_avg = mean(centr_avg);
        guest_avg = [guest_avg,prev_avg];
    end
    [M,I] = min(guest_avg);
    index = I;
    index_guess = [index_guess,index];
end
disp(index_guess)
    



 