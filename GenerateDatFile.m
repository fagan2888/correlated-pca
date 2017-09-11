%%script to convert matrix into dat file for tikz phase transition plots


clear
clc
load('data/phase_trans_vs_n_mc50')
%%X is data matrix of size (numi, numj)
X = PhaseTrans;
%rrange = [1 : 10];
% nrange = unique(ceil(linspace(50, 1000, 10)));

nrange = unique(ceil(linspace(50, 1000, 10)));
[numi, numj] = size(X);
fileID = fopen('PhaseTransvsr_gaussian_Bt.dat', 'w');

for ii = 1 : numi
    for jj = 1 : numj
        fprintf(fileID, '%d %d %d\n', nrange(ii), AlRange(jj), X(ii, jj));
    end
    fprintf(fileID, '\n');
end

fclose(fileID);


%%script to convert simple plot into data file

clear
clc

SE_temp = SE_theory(8, :);
max_SE_temp = max(temp1, [], 1);
mean_SE_temp = mean(temp1, 1);

fileID = fopen('bound_SE.dat', 'w');

for ii = 1 : length(AlRange)
        fprintf(fileID, '%d, %.4f, %.4f,%.4f\n', AlRange(ii), mean_SE_temp(ii),  max_SE_temp(ii), SE_temp(ii));
    fprintf(fileID, '\n');
end

fclose(fileID);
