function out = SFB_SMT(in, M, pulse_shape)
% Description:
%   This function implements the synthesis filter bank that consists 
%   of simple multipliers, IFFT, polyphase filtering, and parallel to serial
%   converter (upsamplers by M/2 and delay chain).
% 
% Usage: 
%   out = SynthesisFB(M, K, opt, in)
%
% Input variables:
%     M: number of subchannels
%    in: data matrix of dimension Mx2SymLen 
%
% Output variables:
%   out: output vector
% 
% Written by: Markku Renfors
% Modified by: Amir Aminjavaheri

K = round(length(pulse_shape)/M);

persistent PRS_InitFlag; % Helper variable, only used for initialisation of OTHER persistent variables
persistent PRS_hp; % prototype filter 
persistent PRS_G;  % polyphase filters
persistent PRS_mul; % extra multipliers

if isempty(PRS_InitFlag) % After first "persistent PRS_InitFlag", the variable has the value "[]".
    PRS_InitFlag = 0;
    PRS_hp = pulse_shape; % generate prototype filter
    PRS_G = reshape(PRS_hp,M,K); % MxK matrix 
    PRS_mul = exp(1i*2*pi*(0:M-1)'/M*(K*M)/2);        
end 

filtermemory = zeros(M,2*(K-1)); % polyphase filter memory
PtoSmemory = zeros(M/2,1); % P/S memory

% Add extra zeros so that whole data go through the polyphase filters
extra = 2*(K-1)+1;
in = [in zeros(M,extra)]; 

Dlen = size(in,2); % data length
out  = zeros(1,Dlen*M/2); % output vector
for k = 1:Dlen % block by block processing
    
    % Multiplication by simple multipliers
    OutMul = in(:,k).*PRS_mul;

    % IFFT
    OutIFFT = ifft(OutMul)*M;

    % Polyphase filtering
    FilterStructure = [OutIFFT filtermemory];
    tmp = FilterStructure(:,1:2:end).*PRS_G;
    OutPF = conj((sum(conj(tmp')))');
    filtermemory = FilterStructure(:,1:end-1);

    % Parallel to serial conversion
    out((1:M/2)+(k-1)*M/2) = reshape(OutPF(1:M/2),1,M/2) + reshape(PtoSmemory,1,M/2);
    PtoSmemory = OutPF(M/2+1:end);       

end 
