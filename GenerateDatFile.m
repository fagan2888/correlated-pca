%%basic script to convert a given matrix into a format that TikZ can use to
%%generate a 3d surf plot. 

clear
clc

%load the desired mat file
%load('data/phase_trans_vs_n_mc50')


%%X is data matrix of size (numi, numj)

X = PhaseTrans;
[numi, numj] = size(X);

%co-ordinate axis labels
% rrange = [1 : 10];
% nrange = unique(ceil(linspace(50, 1000, 10)));


FileName = 'PhaseTransvsn_gaussian_rv_n.dat';
fileID = fopen(Filename, 'w'); %can vary read/write/overwrite if needed

%note here that all entries are integers
for ii = 1 : numi
    for jj = 1 : numj
        fprintf(fileID, '%d %d %d\n', nrange(ii), AlRange(jj), X(ii, jj));
    end
    fprintf(fileID, '\n');
end
fclose(fileID);


%%Part --2:

%%basic script to obtain a table to generate TikZ line-plot

clear
clc

%load desired mat file
%load('data/SE_bnd.mat');

SE_temp = SE_theory(3, :);
temp1 = all_errors{3};
max_SE_temp = max(temp1, [], 1);
mean_SE_temp = mean(temp1, 1);

FileName = 'bound_SE_n100r10.dat';
fileID = fopen(FileName, 'w');

%note here that entries may be floats -- choose precision
for ii = 1 : length(AlRange)
        fprintf(fileID, '%d, %.4f, %.4f,%.4f\n', AlRange(ii), ...
            mean_SE_temp(ii),  max_SE_temp(ii), SE_temp(ii));
    fprintf(fileID, '\n');
end

fclose(fileID);
