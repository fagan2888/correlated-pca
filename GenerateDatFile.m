%%script to convert matrix into dat file for tikz phase transition plots


clear
clc
%load('data/phase_trans_vs_n_mc50')
%%X is data matrix of size (numi, numj)
X = PhaseTrans;
[numi, numj] = size(X);
fileID = fopen('PhaseTransvsr.dat', 'w');

for ii = 1 : numi
    for jj = 1 : numj
        fprintf(fileID, '%d %d %d\n', ii, AlRange(jj), X(ii, jj));
    end
    fprintf(fileID, '\n');
end

fclose(fileID);