%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	Author(s): C. Nicolas Barati nicobarati@rice.edu 
%		Rahman Doost-Mohamamdy: doost@rice.edu
%
% Multiple iterations of a single-shot transmissions from one client or UE
% to one base station radio (UE stands for User Equipment).
% The script explores Bit Error Rate (BER) as a function of Signal-to-Noise
% Ratio (SNR) and therefore iterates over different SNR values (sim_SNR_db
% variable). Within each iteration, only a single frame transmission takes
% place.
%
% We define two modes: OTA (Over-the-air) and SIM_MOD (simulation).
% In simulation mode we simply use a Rayleigh channel whereas the OTA mode
% relies on the Iris hardware for transmission and reception.
% In both cases the client transmits an OFDM signal that resembles a
% typical 802.11 WLAN waveform. If the transmission is OTA, then the user
% specifies a schedule that tells the client when to transmit its frame
% The base station initiates the schedule by sending a beacon signal that
% synchronizes the client. After that, the client will simply transmit its
% frame.
%
%---------------------------------------------------------------------
% Original code copyright Mango Communications, Inc.
% Distributed under the WARP License http://warpproject.org/license
% Copyright (c) 2018-2019, Rice University
% RENEW OPEN SOURCE LICENSE: http://renew-wireless.org/license
% ---------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;

% [version, executable, isloaded] = pyversion;
% if ~isloaded
%     pyversion /usr/bin/python
%     py.print() %weird bug where py isn't loaded in an external script
% end

% Params:
WRITE_PNG_FILES         = 0;           % Enable writing plots to PNG
CHANNEL                 = 11;          % Channel to tune Tx and Rx radios

SIM_MOD                 = 1;

if SIM_MOD
    chan_type               = "rayleigh";
    nt                      = 1;
    sim_SNR_db              = 30;1:20;
    nsnr                    = length(sim_SNR_db);
    snr_plot                = 20;
    TX_SCALE                = 1;         % Scale for Tx waveform ([0:1])
else
    nt                      = 1;
    nsnr                    = 1;
    TX_SCALE                = 1;         % Scale for Tx waveform ([0:1])
    chan_type               = "iris";
end
ber_SIM = zeros(nt,nsnr);           % BER
berr_th = zeros(nsnr,1);            % Theoretical BER
fprintf("Channel type: %s \n",chan_type);

%Iris params:
N_BS_NODE 		= 1;
N_UE 			= 1;
TX_FRQ                  = 3.6e9;
RX_FRQ                  = TX_FRQ;
TX_GN                   = 75;
RX_GN                   = 65;
SMPL_RT                 = 5e6;
N_FRM                   = 10;

if TX_GN > 81
    display('WARNING: MAXIMUM TX GAIN IS 81!');
    TX_GN = 81;
end

bs_ids = string.empty();
bs_sched = string.empty();
ue_ids = string.empty();
ue_sched = string.empty();


% Waveform params
N_OFDM_SYM              = 2;         % Number of OFDM symbols for burst, it needs to be less than 47
MOD_ORDER               = 16;           % Modulation order (2/4/16/64 = BSPK/QPSK/16-QAM/64-QAM)

% OFDM params
SC_IND_PILOTS           = [8 22 44 58];                           % Pilot subcarrier indices
SC_IND_DATA             = [2:7 9:21 23:27 39:43 45:57 59:64];     % Data subcarrier indices
N_SC                    = 64;                                     % Number of subcarriers
CP_LEN                  = 16;                                     % Cyclic prefix length
N_DATA_SYMS             = N_OFDM_SYM * length(SC_IND_DATA);       % Number of data symbols (one per data-bearing subcarrier per OFDM symbol)
N_LTS_SYM               = 2;                                      % Number of 
N_SYM_SAMP              = N_SC + CP_LEN;                          % Number of samples that will go over the air
N_ZPAD_PRE              = 90;                                     % Zero-padding prefix for Iris
N_ZPAD_POST             = N_ZPAD_PRE - 14;                         % Zero-padding postfix for Iris

% Rx processing params
FFT_OFFSET                    = 16;          % Number of CP samples to use in FFT (on average)
DO_APPLY_PHASE_ERR_CORRECTION = 1;           % Enable Residual CFO estimation/correction

%% Define the preamble

% LTS for fine CFO and channel estimation
% lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 ...
%     1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
lts_f = [-1 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 ...
    1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
lts_t = ifft(lts_f, 64); %time domain

preamble = [lts_t(33:64) lts_t lts_t];

%% Generate a payload of random integers
tx_data = randi(MOD_ORDER, 1, N_DATA_SYMS) - 1;

tx_syms = mod_sym(tx_data, MOD_ORDER);

% Reshape the symbol vector to a matrix with one column per OFDM symbol
tx_syms_mat = reshape(tx_syms, length(SC_IND_DATA), N_OFDM_SYM);

% Define the pilot tone values as BPSK symbols
pilots = [1 1 -1 1].';

% Repeat the pilots across all OFDM symbols
pilots_mat = repmat(pilots, 1, N_OFDM_SYM);

%% Time Domain

% Do yourselves: construct the TD input matrix
fdd_mat = zeros(N_SC, N_OFDM_SYM);

% Insert the data and pilot values; other subcarriers will remain at 0
fdd_mat(SC_IND_DATA, :)   = tx_syms_mat;
fdd_mat(SC_IND_PILOTS, :) = pilots_mat;

% Do yourselves: get TD samples:
tdd_tx_payload_mat = ifft(fdd_mat, N_SC, 1);

% Insert the cyclic prefix
if(CP_LEN > 0)
    % Do yourselves: Insert CP
    tx_cp = tdd_tx_payload_mat((end-CP_LEN+1 : end), :);
    tdd_tx_payload_mat = [tx_cp; tdd_tx_payload_mat];
end

% Reshape to a vector
tx_payload_vec = reshape(tdd_tx_payload_mat, 1, numel(tdd_tx_payload_mat));

load ('test1.mat');
  L=4;
 At = TD_est_kongmimoV2(pilot_pos, pilot_set, h2, M, K, ones(K,1)*L,(1:M).' );   %Generate estimation matrix

% Construct the full time-domain OFDM waveform
tx_vec = [zeros(1,N_ZPAD_PRE) preamble tx_payload_vec zeros(1,N_ZPAD_POST) X1.'/10 zeros(1,N_ZPAD_POST)];

% Leftover from zero padding:
tx_vec_iris = tx_vec.';
% Scale the Tx vector to +/- 1
tx_vec_iris = TX_SCALE .* tx_vec_iris ./ max(abs(tx_vec_iris));

for isnr = 1:nsnr
    for it = 1:nt
if (SIM_MOD)

    % Iris nodes' parameters
    bs_ids = ones(1, N_BS_NODE);
    ue_ids = ones(1, N_UE);
    n_samp = length(tx_vec_iris);
    bs_sdr_params = struct(...
        'id', bs_ids, ...
        'n_sdrs', N_BS_NODE, ...        % number of nodes chained together
        'txfreq', [], ...
        'rxfreq', [], ...
        'txgain', [], ...
        'rxgain', [], ...
        'sample_rate', [], ...
        'n_samp', n_samp, ...          % number of samples per frame time.
        'n_frame', [], ...
        'tdd_sched', [], ...     % number of zero-paddes samples
        'n_zpad_samp', N_ZPAD_PRE ...
        );

    ue_sdr_params = bs_sdr_params;
    ue_sdr_params.id =  ue_ids;
    ue_sdr_params.n_sdrs = N_UE;
    ue_sdr_params.txgain = [];

    rx_vec_iris = getRxVec(tx_vec_iris, N_BS_NODE, N_UE, chan_type, sim_SNR_db(isnr), bs_sdr_params, ue_sdr_params, []);
    %rx_vec_iris = getRxVec(tx_vec_iris, N_BS_NODE, N_UE, chan_type, sim_SNR_db(isnr));
    rx_vec_iris = rx_vec_iris.'; % just to agree with what the hardware spits out.
    
else
%% Init Iris nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the Iris experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % Create two Iris node objects:
    bs_ids = ["RF3E000246"];
    ue_ids = ["RF3E000119"];
    
    bs_sched = ["BGGGGGRG"];           % BS schedule
    ue_sched = ["GGGGGGPG"];               % UE schedule
    
    n_samp = length(tx_vec_iris);    
    
    % Iris nodes' parameters
    bs_sdr_params = struct(...
        'id', bs_ids, ...
        'n_sdrs', N_BS_NODE, ...
        'txfreq', TX_FRQ, ...
        'rxfreq', RX_FRQ, ...
        'txgain', TX_GN, ...
        'rxgain', RX_GN, ...
        'sample_rate', SMPL_RT, ...
        'n_samp', n_samp, ...          % number of samples per frame time.
        'n_frame', N_FRM, ...
        'tdd_sched', bs_sched, ...     % number of zero-paddes samples
        'n_zpad_samp', N_ZPAD_PRE ...
        );
    
    ue_sdr_params = bs_sdr_params;
    ue_sdr_params.id =  ue_ids;
    ue_sdr_params.n_sdrs = N_UE;
    ue_sdr_params.tdd_sched = ue_sched;
    
    rx_vec_iris = getRxVec(tx_vec_iris, N_BS_NODE, N_UE, chan_type, [], bs_sdr_params, ue_sdr_params, []);

end
    rx_vec_iris = rx_vec_iris.';
    l_rx_dec=length(rx_vec_iris);


%% Correlate for LTS

% Complex cross correlation of Rx waveform with time-domain LTS
a = 1;
unos = ones(size(preamble.'))';
v0 = filter(flipud(preamble'),a,rx_vec_iris);
v1 = filter(unos,a,abs(rx_vec_iris).^2);
m_filt = (abs(v0).^2)./v1; % normalized correlation
lts_corr = m_filt;
[rho_max, ipos] = max(lts_corr);

payload_ind = ipos +1;
lts_ind = payload_ind - N_LTS_SYM*(N_SC + CP_LEN);

% Re-extract LTS for channel estimate
rx_lts = rx_vec_iris(lts_ind : lts_ind+159);
rx_lts1 = rx_lts(-64+-FFT_OFFSET + [97:160]);
rx_lts2 = rx_lts(-FFT_OFFSET + [97:160]);

% Received LTSs
% Do yourselves: take the FD pilots:
rx_lts1_f = fft(rx_lts1);
rx_lts2_f = fft(rx_lts2);

% Do yourselves: Calculate channel estimate from average of 2 training symbols: 
rx_H_est = mean([rx_lts1_f./lts_f.'   rx_lts2_f./ lts_f.'], 2); 

%% Rx payload processing
% Extract the payload samples (integer number of OFDM symbols following preamble)
if( (length(rx_vec_iris) - payload_ind ) > (N_SYM_SAMP * N_OFDM_SYM) )
    payload_vec = rx_vec_iris(payload_ind : payload_ind + (N_SYM_SAMP * N_OFDM_SYM),:);
else
    payload_vec = rx_vec_iris(payload_ind : end,:);
end
missed_samps = (N_SC+CP_LEN) * N_OFDM_SYM - length(payload_vec); %sometimes it's below 0.

if (missed_samps > 0) 
    payload_vec = [payload_vec; zeros(missed_samps,1)];
elseif (missed_samps < 0)
    payload_vec = payload_vec(1:end+missed_samps,:);
end

% % % % save('Y.mat');%,'rx_H_est','payload_vec');
% % % % 
% % % % 
% % % %  load('Y.mat');
% load('test1.mat');

 ber=1;
Y = rx_vec_iris( payload_ind + (N_SYM_SAMP * N_OFDM_SYM)+N_ZPAD_POST:payload_ind + (N_SYM_SAMP * N_OFDM_SYM)+N_ZPAD_POST+length(X1),:);
%   L=10;
%  At = TD_est_kongmimoV2(pilot_pos, pilot_set, h2, M, K, ones(K,1)*L,(1:M).' );   %Generate estimation matrix
  CSI=0;
    MSE_arr = zeros(length(N_arr),K);
    ber_arr = zeros(length(N_arr),K);
  
    
    N   = N_arr(end);  % Number of BS antennas
    
    %        beta = cell_free_setup(N, K, cell_side, 'Grid'); %Pathloss coefficients based on realistic scenario
     %%% Exponential channel PDP
% %     alpha_vec = .05*5*(1:10)';    %alpha based on a linear increase for users
% %     pdp_arr   = zeros(L,K,N);
% %     for ik = 1:K
% %         for jk = 1:N
% %              
% % %             alpha_vec(ik) =  rand/2;
% %             pdp_arr(:,ik,jk) = exp_pdp(L,alpha_vec(ik));
% % % %             pdp_arr(:,ik,jk) = exp_pdp(L,alpha_vec(ik))+[zeros(L/2,1); exp_pdp(L/2,alpha_vec(ik)*4)]/16  ;  % * beta(ik,jk); %rand*ik; %
% % %             pdp_arr(:,ik,jk) = [exp_pdp(L/2,alpha_vec(ik)); zeros(L/2,1)]+[zeros(L/2,1); exp_pdp(L/2,alpha_vec(ik)*2)]/1.3 ;
% % %             [maxx, indxx] = max(pdp_arr(:,ik,jk));
% % %             pdp_arr(:,ik,jk) = circshift(pdp_arr(:,ik,jk),indxx-1);
% %              pdp_arr(:,ik,jk) = pdp_arr(:,ik,jk)/sum(pdp_arr(:,ik,jk),'all');
% %         end
% %     end
% %     pdp_arr = repmat(pdp_arr,1,1,N_ap);
% %       h = rayleigh_chan_gen(pdp_arr,N*N_ap,K);
    
%         %%% LTE channel models

pdp_arr   = zeros(L,K,N);
h   = zeros(L,N,K);

%  for ik = 1:K
%          for jk = 1:N
%         fs = 1.92 * 10^6;
%         fs = 3.84 * 10^6*4;
% %         fs = 3.84 * 10^6;
% %         [h, pdp]= lte_channel_gen(fs,'ETU',N_arr(end),K);
% % % % % % % %                     [h, pdp] = rand_channel_gen(N_arr(end),K,L);
%      [h(:,:,ik),  pdp_arr(:,ik,:)] = lte_channel_gen(fs,'TDL-C',N_arr(end),1);
%          end
%  end
% %          fs = 3.84 * 10^6;
% %          [h, pdp]= lte_channel_gen(fs,'ETU',N_arr(end),K);
% %         pdp_arr = squeeze(pdp(:,1,:));
% %             pdp_arr = repmat(pdp_arr,1,1,N_ap*N_arr(end));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %h = ones(size(h));
%     %power control
%     %power_cf = frac_power_control(beta,K,N,max_power,nu); %Detemine power coefficents
%     X= X1; % .* (power_cf).';
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Simulation using data
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Y = zeros(size(X,1)+L-1,N*N_ap);
%     Y = Y + mean(std(X))*noise_std/sqrt(2) * (randn(size(Y)) + 1i*randn(size(Y))); % AWGN
% %                  Y = Y + noise_std/sqrt(2) * (randn(size(Y)) + 1i*randn(size(Y))); % AWGN
%     for ik = 1:K
%         for inn=1:N*N_ap
%             Y(:,inn) = Y(:,inn) + conv(X(:,ik) , h(:,inn,ik));
%         end
%     end
%     Y = Y(1:size(X,1),:);
    %Y = awgn(Y,SNR_db,0);
    
    
    %clear X;
    
    % Demodulation
    % Here, we only perform the analysis filter bank and the decimation for
    % each antenna chain. We DO NOT apply any equalization as well as the
    % phase compensation here.
    %    mod_pulse = mod_pulse_estimator(pulse_shape,pdp_arr,[]);
%      pdp_est =  pdp_estimate(h(:,1:N,ik),ik,'average'); 
    
%     mo_pulse = modified_prototype(pulse_shape,pdp_arr(:,1,1));
%     R1= zeros(7668,N*N_ap);
% %              R1= zeros(3408,N*N_ap);
%     RP = zeros(6816*2,N*N_ap);
    for inn = 1:N*N_ap
        [R1(:,inn), RP(:,inn)] = SMT_PPN_demod2(Y(:,inn), pulse_shape, nSymb, M, OL, act_subc);
    end
    %clear Y;
    % Peform MSE calculation for each case of antenna numbers in N_arr
    %%% channel estimation goes here
    %        pilot_set_eff = repmat(pilot_set,1,K) .* power_cf.';
    
    %At = TD_est_kongmimoV2(pilot_pos, pilot_set_eff, h2, M, K, ones(K,1)*L,act_subc );   %Generate estimation matrix
    h_est = zeros(size(h));
    if CSI ==0
        for inn=1:N*N_ap
            h_est(:,inn,:) = channel_estiation(RP(1:(fl_delay+2)*M,inn), pilot_pos, At, M, K, L, OL, (1:M).',Vtt);
        end
        h = ifft(rx_H_est);
        h_est = [h_est; zeros(M-L,1)];
        [max_i, ind_i] = max(abs(h));
%         h_est = circshift(h_est, ind_i-1)*10;
        MSE = var(squeeze(h)-squeeze(h_est))/var(squeeze(h)) ;
% %         for ikk=1:K
% %             for innn=1:N
% %                 MSE = var(squeeze(h(:,innn,ikk))-squeeze(h_est(:,innn,ikk)))/var(squeeze(h(:,innn,ikk))) + MSE;
% %                 if ikk==1
% %                 error_vec(:,iter) = squeeze(h(:,innn,ikk))-squeeze(h_est(:,innn,ikk));
% %                 else
% %                 error_vec2(:,iter) = squeeze(h(:,innn,ikk))-squeeze(h_est(:,innn,ikk));
% %                 end
% % 
% %             end
% %         end
    end
    % %                 for ikk=1:K
    % %                     for inn=1:N*N_ap
    % %                         h_est(:,inn,ikk) = h(:,inn,ikk) + mean(std(h(:,inn,ikk)))* ...
    % %                             est_noise/sqrt(2) * (randn(size(h(:,inn,ikk))) + 1i*randn(size(h(:,inn,ikk)))); % AWGN
    % %                     end
    % %                 end
    if CSI == 1
        hfreq = fft(h,M,1);
    else if CSI== 0
            hfreq = fft(h_est,M,1);
        end
    end
    hfreq2 = fft(h,M,1);
    hfreq2 = hfreq2(act_subc,:,:);
    hfreq2 = permute(hfreq2,[2 3 1]);
    
    hfreq = hfreq(act_subc,:,:);
    hfreq = permute(hfreq,[2 3 1]);
    
 %%   
    for in = 1:length(N_arr)
        
        N   = N_arr(in);  % Number of BS antennas
%                    N=1000;
        % MIMO equalization using channel frequency gains (MRC, ZF, MMSE)
        R2 = eq_mimo(R1(:,1:N*N_ap), hfreq(1:N*N_ap,:,:), SNR_db, eq_mode);
        %R2 = squeeze(R1(:,1,1));
        % Consider terminal k
        for ik = 1:K
            R2_k = R2(:,ik);
            MSET = ones(20,1)*10;
            bert = ones(20,1)*10;
            for delay =  [Lphi]; %[9, 11, 13, 15 , 17, 19, 21, 23] %:2:17
                if EnEQ
                    %pdp_est = pdp_arr(:,ik)/sum(pdp_arr(:,ik));
                    if pdp_CSI ==1
                        if N*N_ap>0
                            %[p2, sort_i] = sort(squeeze(sum(pdp_arr(:,ik,1:N),1)));
                            %pdp_arr(:,ik,1:N) = pdp_arr(:,ik,sort_i);
                            pdp_est = sum(pdp_arr(:,ik,1:N*N_ap),3)/ sum(sum(pdp_arr(:,ik,1:N*N_ap),3));  %PerfectCSI Average PDP
%                              [hhhh,  pdp_est] = lte_channel_gen_2(fs,'TDL-C',1,1);
                        else
                            pdp_est = sum(pdp_arr(:,ik,1:7),3)/sum(sum(pdp_arr(:,ik,1:7),3));  %PerfectCSI Average PDP
                        end
                    else if pdp_CSI ==0
                            if CSI ==0
                                %                                         pdp_est = squeeze(pdp_arr(:,ik,1));
                                
                                pdp_est =  pdp_estimate(h_est(:,1:N,ik),ik,pdp_est_type);
                                
                            else if CSI ==1
                                      pdp_est =  pdp_estimate(h(:,1:N,ik),ik,pdp_est_type);
                                end
                            end
                        end
                            end
                              
%                               pdp_est = exp_pdp(L,.1);
%                     pdp_est = pdp_est + randn(size(pdp_est)) * 0.01;
                    
                    
% %                     pdp_est = conv(pdp_arr(:,1),conv(pulse_shape, pulse_shape));
  switch post_eq_mode 
    case 'eq_channel'
                                        if CSI == 0  %, err(iter, in, ik)
                                         [R3_k] = persubcarrier_eq_mu_est(R2_k,M,squeeze(h_est(:,1:N*N_ap,:)),...
                                             act_subc,pulse_shape, hfreq(1:N*N_ap,:,:),K,ik,delay,...
                                             noise_std2,eq_mode,noise_std2^2  / 64,squeeze(h(:,1:N*N_ap,:)), compensation); %./ power_cf(ik); error:MSE/L/K
%                                         iii=iii+1;5

                                        else
                                            [R3_k]= persubcarrier_eq_mu(R2_k,M,squeeze(h(:,1:N*N_ap,:)),...
                                             act_subc,pulse_shape, hfreq(1:N*N_ap,:,:),K,ik,delay,...
                                             noise_std2,eq_mode,noise_std2^2  / 64,pdp_est,pdp_arr(:,ik,1)); %./ power_cf(ik);%error: MSE/L/K
                                        end
    case 'pdp_MMSE'

pdp_est = pdp_est ./ sum(pdp_est, 'all');
                                R3_k = persubcarrier_eq_pdp(R2_k,M,squeeze(h_est(:,1:N*N_ap,:)),...
                                             act_subc,pulse_shape,ik,delay,...
                                            noise_std2,noise_std2^2  / 64,pdp_est,compensation); %./ power_cf(ik);
     end

%                         pdp_est = pdp_arr(:,ik,1);
%  pdp_est = pdp_est/sum(pdp_est);
%                            R3_k = pdp_eq(R2_k,M,pdp_est,act_subc, noise_std2^2/64*L ); %./ power_cf(ik);
                    %R3_k = actual_eq_mu(R2_k,M,squeeze(h(:,1:N,:)),act_subc,pdp_est,pulse_shape, hfreq(1:N,:,:),K,ik);
                    %R3_k = symbol_spaced_eq(R2_k,M,squeeze(h(:,1:N,:)),act_subc,pdp_est,pulse_shape, hfreq(1:N,:,:),K,ik);
                    %R3_k = persubcarrier_eq(R2_k,M,squeeze(h(:,1:N,:)),act_subc,pdp_est,pulse_shape, hfreq(1:N,:,:),K,ik,delay);
                    
                    
                else
                    R3_k = R2_k; %./ power_cf(ik);
                end
                R4_k = phase_adj(R3_k,M,nSymb,OL,act_subc);
                R5_k = OQAMpostprocessing(R4_k, Nact, nSymb);
                %MSE_arr(in,ik) = mean(abs((R5_k) - (S(:,ik))).^2); %(1:Nact:end)
                MSET(delay-8) = mean(abs((R5_k((Ng+1)*M+1:end)) - (S((Ng+1)*M+1:end,ik))).^2);
                
                
                if ber == 1
             
                %%%%BER
                            dataSymbolsOut_n = qamdemod(R5_k(Nact*Ng+1:end),nQAM,'bin');
                            dataOutMatrix_n = de2bi(dataSymbolsOut_n,log2(nQAM));
                            dataOut_n = dataOutMatrix_n(:);                   % Return data in column vector
                
                            dataSymbolsOut_n = qamdemod(S(Nact*Ng+1:end,ik),nQAM,'bin');
                            dataOutMatrix_n = de2bi(dataSymbolsOut_n,log2(nQAM));
                            dataIn_n = dataOutMatrix_n(:);                   % Return data in column vector
                
                            [numErrors_n,berrate_n] = biterr(dataIn_n,dataOut_n)
% %                             fprintf('\nEstimation bit error rate = %5.2e, based on %d errors\n', ...
% %                                 berrate_n,numErrors_n);
                            bert(delay-8) = berrate_n;
                end
            end
            [MSE_arr(in,ik), idxx] = min(MSET);
             ber_arr(in,ik) = bert(idxx);
        end
        
    end
    %clear S R1 R2 hfreq ;
    
    SINR_sim = 1./MSE_arr;
    %SINR_arr = cat(3,SINR_arr,SINR_sim);
    SINR_arr(:,:,it) = SINR_sim;
    bert_arr(:,:,it) = ber_arr;
    %%%%%%%%%%%%%%%%%%%%%%
    % Plot the results
    %%%%%%%%%%%%%%%%%%%%%%
    

    
    end
end



        if rem(iter,5) == 0
    
    %%% user of intereset SINR
    % %         user_idx = 1;
    % %         SINR_arr1 = squeeze(SINR_arr(:,user_idx,:));
    % %         pdp1 = pdp_arr(:,user_idx);
    
            SINR_arr1 = squeeze(sum(SINR_arr(:,:,:),2)/K);
    %        pdp1 = pdp_arr(:,user_idx);
    
            sinr_mean_sim = mean(pow2db(SINR_arr1),2);
            sinr_std_sim  = std(pow2db(SINR_arr1),1,2);
    
            figure(aaa); clf;
            hold all;
            %errorbar(N_arr,sinr_mean_sim,sinr_std_sim,'-','LineWidth',2);
            plot(N_arr/N_wrap,sinr_mean_sim,'-','LineWidth',2);
    
    
            %         sinr_analytic = pow2db(FBMC_SINR_saturation_level(M,pdp1,pulse_shape));
            %         plot(N_arr,sinr_analytic*ones(numel(N_arr),1),':k','LineWidth',2);
    
            %         sinr_analytic = pow2db(FBMC_SINR_MRC(M,N_arr,SNR_db,pdp_arr,pulse_shape,user_idx));
            %         plot(N_arr,sinr_analytic,'--x','LineWidth',2);
            %
            %         sinr_analytic = pow2db(FBMC_SINR_ZF(M,N_arr,SNR_db,pdp_arr,pulse_shape,user_idx));
            %         plot(N_arr,sinr_analytic,'--*','LineWidth',2);
    
            set(gca,'XScale','log');
            xlim([1 N_arr(end)]); grid on;
            xlabel('N'); ylabel('SINR (dB)');
            legend('Sim','Sat','MRC','ZF','Location','SouthEast');
            grid on; box on; hold off; drawnow;
    
        end

scatterplot(R5_k)


fprintf(['\n' datestr(datetime('now')) '\n']);
fprintf('Simulation time %.2f hours\n', toc / 3600);