%% Grab data
num_test = 8;
num_train = 11;
test_set = cell(num_test,2);
for i = 1:num_test
    filename = ['s',num2str(i)];
    test_set{i,1} = join(['s',num2str(i)]);
    test_set{i,2} = join(['Fs',num2str(i)]);
    [test_set{i,1},test_set{i,2}] = audioread([filename,'.wav']);
    test_set{i,1} = (test_set{i,1}-mean(test_set{i,1}))./max(abs(test_set{i,1})); %Normalization to 1
end
%% Remove stero
% for i=9:num_data
%     test_set{i,1}(:,2) = [];
% end

%% 
codebook = readmatrix("Codebook_final.txt");

%% %% Generate Ceptsrum Coefficients 
N=256; %Length of Frame Block
M=100; %Shift Between Each  Frame Block
Window=.54-.46*cos(2*pi().*[1:N]/(N-1));
for i = 1:num_test
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
for i=1:num_train
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
for i=1:num_test
    col = nnz(size_(i,:));
    length_row = [length_row,col];
end


%% find out the length of the code book
size_ = (cellfun('size',code_cell,1));
len_code_cell =[];
for i=1:num_train
    col = nnz(size_(i,:));
    len_code_cell = [len_code_cell,col];
end
len_code_cell = len_code_cell(1);

%% get the euc dist and match speakers
guess =[];
for i=1:num_test
    tot_min_avg =[];
    for f=1:num_train
        dist_mat =[];
        for l = 1:len_code_cell
            reg =[];
            for j = 1:length_row(i)
                euc_dist = disteu(code_cell{f,l}, guest_speaker{i,j});
                reg = [reg,euc_dist];     
            end 
            dist_mat(l,:) = reg;
        end 
        min_dist = min(dist_mat);
        min_dist = mean(min_dist);
        tot_min_avg = [tot_min_avg,min_dist];
    end
    [M,I] = min(tot_min_avg);
    guess = [guess,I];
end
for i=1:length(guess)
    X = ['Guest Speaker ' num2str(i),' is known speaker  ',...
        num2str(guess(i)),'.'];
    disp(X)
end




 