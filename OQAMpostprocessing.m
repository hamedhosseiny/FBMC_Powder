function S2 = OQAMpostprocessing(S, nAct, nSymb)
% OQAM processing

S = reshape(S, nAct, 2*nSymb);

SR = real(S(:,1:2:end));
SI = real(S(:,2:2:end));
S2 = SR + 1i*SI;
S2 = S2(:);