function yt = QAM_mod(ut, M,K)
% u: vector of input bits
% M: Size of signal constellation, e.g., 2, 4, 16, ...
% Note: number of elements in u must be an integer multiple of log2(M)
persistent QAMmod
if isempty(QAMmod)
    QAMmod = comm.RectangularQAMModulator('ModulationOrder', M, ...
                 'BitInput', true, ...
                 'SymbolMapping', 'Gray', ...
                 'NormalizationMethod', 'Average power');
end

for ik = 1:K
    u = ut(:,ik);
y = step(QAMmod, u);
yt(:,ik) = y;
end