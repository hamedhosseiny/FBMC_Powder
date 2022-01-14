function [Xt, st1] = SMT_PPN_mod_pilot_mu(St, pulse_shape, nSymb, N, act_subc, K, pilots, power_cf, max_power)
% SMT Transmitter
st1 = zeros(64,nSymb*2,K);
for ik = 1:K
S = St(:,ik) * power_cf(ik);

nAct = numel(act_subc);

S = reshape(S, nAct, nSymb); % Serial to Parallel
SR = real(S);
SI = imag(S);

% OQAM processing
sIn = zeros(N, 2*nSymb); % sIn is the input matrix to the IDFT block

sIn(act_subc,1:2:2*nSymb) = SR; 
sIn(act_subc,2:2:2*nSymb) = SI; 
sIn(1:N) = pilots(:,ik) * max_power;
% phase adjustment
% phase_mat = 1i*ones(size(sIn));
% phase_mat(1:2:end,1:2:end) = 1;
% phase_mat(2:2:end,2:2:end) = 1;


phase_mat = 1i*ones(size(sIn));
for ii =1:size(phase_mat,1)
    for jj=1:size(phase_mat,2)
        phase_mat(ii, jj) = phase_mat(ii, jj)^(ii+jj-2);
    end
end
sIn = phase_mat .* sIn; 
st1(:,:,ik) = sIn;
X = SFB_SMT(sIn, N, pulse_shape);
X = X(:);

Xt(:,ik) = X;
end