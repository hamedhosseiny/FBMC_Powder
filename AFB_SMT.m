function out = AFB_SMT(in, M, pulse_shape)
%   This function implements the analysis filter bank that consists of
%   serial to parallel converter (delay chain and downsamplers by M/2), 
%   polyphase filtering, FFTs, and imple multipliers.
% 
% Usage: 
%   out = AFB_SMT(M, K, opt, in)
%
% Input variables:
%     M: number of subchannels
%     K: overlapping factor, i.e., L=K*M
%
% Output variables:
%   out: output matrix
% 
% Written by: Markku Renfors
% Modified by: Amir Aminjavaheri

K = round(length(pulse_shape)/M);

persistent PRS_InitFlag; % Helper variable, only used for initialisation of OTHER persistent variables
persistent PRS_hp; % prototype filter
persistent PRS_Q; % polyphase filters
persistent PRS_mul; % extra multipliers

if isempty(PRS_InitFlag) 
    PRS_InitFlag = 0;
    PRS_hp = pulse_shape; 
    PRS_Q = fliplr(reshape(PRS_hp,M,K)); % MxK matrix
    PRS_mul = exp(-1i*2*pi*(0:M-1)'/M*(K*M)/2);        
end 

filtermemory = zeros(M,2*(K-1)); % polyphase filter memory

% This makes sure that data length is a multiple of M/2
len = length(in); 
extra = ceil(len/(M/2))*M/2-len;

% Add extra zeros so that whole data go through the polyphase filters
extra2 = (2*(K-1)+1)/2*M;

in = [in zeros(1,extra+extra2)]; 
Dlen = length(in)/(M/2)-1;

out = zeros(M,Dlen); % output matrix

for k = 1:Dlen % block by block processing

    % Serial to parallel conversion
    in2 = in((1:M)+(k-1)*M/2);
    OutStoP = in2(:);

    % Polyphase filtering
    FilterStructure = [OutStoP filtermemory];
    tmp = FilterStructure(:,1:2:end).*PRS_Q;
    OutPF = conj((sum(conj(tmp')))');
    filtermemory = FilterStructure(:,1:end-1);

    % FFT
    OutFFT = fft(OutPF);

    % Multiplication by simple multipliers
    out(:,k) = OutFFT.*PRS_mul;

end 

