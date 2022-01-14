function [S,received_pilot] = SMT_PPN_demod2(Y, pulse_shape, nSymb, N, K, act_subc)

if size(Y,1) ~= 1
    Y = Y.'; % Transform Y to a row vector
end

% IFFT output matrix 
sOut = AFB_SMT(Y, N, pulse_shape);
received_pilot = sOut(:); 
S    = sOut(act_subc,:);
S    = S(:);

% % Phase adjustment
% phase_mat = -1i*ones(size(sOutEq));
% phase_mat(1:2:end,1:2:end) = 1;
% phase_mat(2:2:end,2:2:end) = 1;
% sOutEq = sOutEq .* phase_mat;
% 
% % Remove the transient
% delay = 2*(K-1);
% S = sOutEq(act_subc, delay + (1:2*nSymb)); 
% S = S(:);