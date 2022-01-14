%%CEll_free after new scenario
%Modified for SIR last time
%clc; 
clear all;
aaa = 8;
% Paraetersrr
M         = 64;  % Number of subcarriers
K         = 2;   % Number 1of users
eq_mode   = 'ZF';  % 'MF', 'ZF', 'MMSE', 'MFDC'
nSymb     = 100;  % Number of OFDM symbols to transmit at each iteration
nIter     = 1;  % Averaging the results over nIter channel realizations
N_wrap    = 1;
N_arr     = 10;[K+1 2.^(3:1:9)+1  1000]; % Run the sims for an array of N (number of antennas)
N_ap      = 1;    %Number of antennas at each access point
SNR_db    = 10;
OL        = 4;    % FBMC Overlapping factor
EnEQ      = 0;    % FBMC PDF equalizer
ber       = 0;    % Bit error rate
max_power = 1;    %maximum transmit power in watts
nu = 0;           %Power control exponent
CSI =0; %CSI =1 when perfect CSI and CSI = 0 when imperfect
pdp_CSI = 1;
pdp_est_type = 'average'; %Can take alpha, alpha_average, average
post_eq_mode = 'pdp_MMSE'; %Can take pdp_MMSE, eq_channel
compensation = 0 ; %zero or one
fl_delay  = 2*(OL-1);    %transient time
Lphi=11;
nQAM     = 4; % QAM size
Nact     = (floor(round(75/128*M)/4)*4); % Number of active subcarriers
L        = 16; % Nact/K;   %round(M/6);
act_subc = [1:ceil(Nact/2) (ceil(-Nact/2):-1)+M]' + 1; % Active subcarriers   %
% % act_subc = (1:M);
% pulse_shape = GeneratePulse(M, OL, 'SRRC', 0.9);
% h2 = pulse_shape; %modify pulse_shaping filter for channel estimation proccess
% pulse_shape(end) = [];
pulse_shape = GeneratePulse(M, OL, 'PHYDYAS', 1);
h2 = [pulse_shape; pulse_shape(end)];
% % % Standard deviation of AWGN
SNR_db2 = SNR_db + pow2db(Nact/M);
noise_std = 1 / sqrt(db2pow(SNR_db2));
noise_std2 = 1 / sqrt(db2pow(SNR_db));
%noise_std = sqrt(290 * 20*10^6 * 10^(0.9) * 1.3 * 10^-23 );
est_noise = 1 / sqrt(db2pow(SNR_db+5));
% for saving the results
filename = ['SINR_' eq_mode '_EQ' num2str(EnEQ)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iteration over different channel realizations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
minsnr=0;
maxsnr=0;
SINR_arr = [];
power=0;
%%%%%%%%%%%
% Pilot Asignment
Ng = 2;   %Number of guard symbols
Guards = zeros(Ng *Nact,1);   %generate guard symbol matrix
pilot_set = 2 * (randi([0,1], M/K, K)-0.5); %Assign Random pilot set
pilot_set = pilot_set *sqrt(2); % * sqrt(M/K);
% % pilot_set(1) = 1;
% % pilot_set(3) = 1;
% % pilot_set(2) = 1;
%pilot_set = ones(size(pilot_set));
[pilot_pos, pilots] = pilot_assign(K, M,pilot_set,act_subc);

% % pilot_pos = [1; 10; 19; 28];
% % pilots = zeros(Nact,K);
% % pilots(1,1) = pilot_set;
% % pilots(10,2) = pilot_set;
% % pilots(19,3) = pilot_set;
% % pilots(28,4) = pilot_set;
if CSI == 0
    At = TD_est_kongmimoV2(pilot_pos, pilot_set, h2, M, K, ones(K,1)*L,(1:M).' );   %Generate estimation matrix
end
SINR_arr = zeros(length(N_arr),K,nIter);
bert_arr = zeros(length(N_arr),K,nIter);

% Parallel pool midification
%    parpool(8); delete(gcp);
%distcomp.feature( 'LocalUseMpiexec', false )
%%
S = zeros(nSymb*Nact,K);
bIn = zeros(Nact * nSymb * log2(nQAM),K);
%%%%%%
%%%%Noise correaltion matrix
Vtt = eye(M);
temp = pilot_pos.';
pilot_post =temp(:);
imag_int = [  1i*.239] * 2;

for ii = 1:M
    for jj =1:M
        if pilot_post(jj) -pilot_post(ii) ==-1
            Vtt(ii,jj) = imag_int(1);
        end
        if pilot_post(jj) -pilot_post(ii) == 1
            Vtt(ii,jj) = -1*imag_int(1);
        end
    end
end
Vtt = Vtt.' * noise_std2^2 * .78^2;


%%
MSE =0;
eqch_error = zeros(K,length(N_arr),nIter,L);
err=zeros(nIter,numel(N_arr),K);
iii=1;
error_vec = zeros(L,nIter);
error_vec2 = zeros(L,nIter);

% % %         %%% LTE channel models for constant channels

pdp_arr   = zeros(L,K,N_arr(end));
h   = zeros(L,N_arr(end),K);

%  for ik = 1:K
% %          for jk = 1:N
%         fs = 1.92 * 10^6;
%         fs = 3.84 * 10^6*4;
% %         fs = 3.84 * 10^6;
% % %         [h, pdp]= lte_channel_gen(fs,'ETU',N_arr(end),K);
% % % % % % % % %                     [h, pdp] = rand_channel_gen(N_arr(end),K,L);
%      [h(:,:,ik),  pdp_arr(:,ik,:)] = lte_channel_gen(fs,'TDL-C',N_arr(end),1);
% %          end
%  end
for iter = 1:nIter
    % Transmitter
    bIn     = randi([0,1], Nact * nSymb * log2(nQAM), K);
    S =  QAM_mod(bIn, nQAM,K);
    %S(1:Nact,:) = pilots;
    S(1:Nact*(Ng),:) = repmat(Guards,1,K);
    % %     pilots = zeros(M,K);
    X1 = SMT_PPN_mod_pilot_mu(S, pulse_shape, nSymb, M, act_subc,K,pilots,ones(1,K),1);
    
    
    
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

 for ik = 1:K
         for jk = 1:N
        fs = 1.92 * 10^6;
        fs = 3.84 * 10^6*4;
%         fs = 3.84 * 10^6;
%         [h, pdp]= lte_channel_gen(fs,'ETU',N_arr(end),K);
% % % % % % %                     [h, pdp] = rand_channel_gen(N_arr(end),K,L);
     [h(:,:,ik),  pdp_arr(:,ik,:)] = lte_channel_gen(fs,'TDL-C',N_arr(end),1);
         end
 end
% %          fs = 3.84 * 10^6;
% %          [h, pdp]= lte_channel_gen(fs,'ETU',N_arr(end),K);
% %         pdp_arr = squeeze(pdp(:,1,:));
% %             pdp_arr = repmat(pdp_arr,1,1,N_ap*N_arr(end));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %h = ones(size(h));
    %power control
    %power_cf = frac_power_control(beta,K,N,max_power,nu); %Detemine power coefficents
    X= X1; % .* (power_cf).';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulation using data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y = zeros(size(X,1)+L-1,N*N_ap);
    Y = Y + mean(std(X))*noise_std/sqrt(2) * (randn(size(Y)) + 1i*randn(size(Y))); % AWGN
%                  Y = Y + noise_std/sqrt(2) * (randn(size(Y)) + 1i*randn(size(Y))); % AWGN
    for ik = 1:K
        for inn=1:N*N_ap
            Y(:,inn) = Y(:,inn) + conv(X(:,ik) , h(:,inn,ik));
        end
    end
    Y = Y(1:size(X,1),:);
    %Y = awgn(Y,SNR_db,0);
    
    
    %clear X;
    
    % Demodulation
    % Here, we only perform the analysis filter bank and the decimation for
    % each antenna chain. We DO NOT apply any equalization as well as the
    % phase compensation here.
    %    mod_pulse = mod_pulse_estimator(pulse_shape,pdp_arr,[]);
%      pdp_est =  pdp_estimate(h(:,1:N,ik),ik,'average'); 
    
    mo_pulse = modified_prototype(pulse_shape,pdp_arr(:,1,1));
    R1= zeros(7668,N*N_ap);
%              R1= zeros(3408,N*N_ap);
    RP = zeros(6816*2,N*N_ap);
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
        
        for ikk=1:K
            for innn=1:N
                MSE = var(squeeze(h(:,innn,ikk))-squeeze(h_est(:,innn,ikk)))/var(squeeze(h(:,innn,ikk))) + MSE;
                if ikk==1
                error_vec(:,iter) = squeeze(h(:,innn,ikk))-squeeze(h_est(:,innn,ikk));
                else
                error_vec2(:,iter) = squeeze(h(:,innn,ikk))-squeeze(h_est(:,innn,ikk));
                end

            end
        end
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
                
                %%
                if ber == 1
             
                %%%%BER
                            dataSymbolsOut_n = qamdemod(R5_k(Nact*Ng+1:end),nQAM,'bin');
                            dataOutMatrix_n = de2bi(dataSymbolsOut_n,log2(nQAM));
                            dataOut_n = dataOutMatrix_n(:);                   % Return data in column vector
                
                            dataSymbolsOut_n = qamdemod(S(Nact*Ng+1:end,ik),nQAM,'bin');
                            dataOutMatrix_n = de2bi(dataSymbolsOut_n,log2(nQAM));
                            dataIn_n = dataOutMatrix_n(:);                   % Return data in column vector
                
                            [numErrors_n,berrate_n] = biterr(dataIn_n,dataOut_n);
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
    SINR_arr(:,:,iter) = SINR_sim;
    bert_arr(:,:,iter) = ber_arr;
    %%%%%%%%%%%%%%%%%%%%%%
    % Plot the results
    %%%%%%%%%%%%%%%%%%%%%%
    
    % % %     if rem(iter,5) == 0
    % % %
    % % % %%% user of intereset SINR
    % % % % %         user_idx = 1;
    % % % % %         SINR_arr1 = squeeze(SINR_arr(:,user_idx,:));
    % % % % %         pdp1 = pdp_arr(:,user_idx);
    % % %
    % % %         SINR_arr1 = squeeze(sum(SINR_arr(:,:,:),2)/K);
    % % % %        pdp1 = pdp_arr(:,user_idx);
    % % %
    % % %         sinr_mean_sim = mean(pow2db(SINR_arr1),2);
    % % %         sinr_std_sim  = std(pow2db(SINR_arr1),1,2);
    % % %
    % % %         figure(aaa); clf;
    % % %         hold all;
    % % %         %errorbar(N_arr,sinr_mean_sim,sinr_std_sim,'-','LineWidth',2);
    % % %         plot(N_arr/N_wrap,sinr_mean_sim,'-','LineWidth',2);
    % % %
    % % %
    % % %         %         sinr_analytic = pow2db(FBMC_SINR_saturation_level(M,pdp1,pulse_shape));
    % % %         %         plot(N_arr,sinr_analytic*ones(numel(N_arr),1),':k','LineWidth',2);
    % % %
    % % %         %         sinr_analytic = pow2db(FBMC_SINR_MRC(M,N_arr,SNR_db,pdp_arr,pulse_shape,user_idx));
    % % %         %         plot(N_arr,sinr_analytic,'--x','LineWidth',2);
    % % %         %
    % % %         %         sinr_analytic = pow2db(FBMC_SINR_ZF(M,N_arr,SNR_db,pdp_arr,pulse_shape,user_idx));
    % % %         %         plot(N_arr,sinr_analytic,'--*','LineWidth',2);
    % % %
    % % %         set(gca,'XScale','log');
    % % %         xlim([1 N_arr(end)]); grid on;
    % % %         xlabel('N'); ylabel('SINR (dB)');
    % % %         legend('Sim','Sat','MRC','ZF','Location','SouthEast');
    % % %         grid on; box on; hold off; drawnow;
    % % %
    % % %         %saveas(gcf, ['../fig/',filename], 'fig');
    % % %         %save(['../data/',filename], 'SINR_arr','M','K','SNR_db','iter','EnEQ');
    % % %
    % % %         fprintf(['\n' datestr(datetime('now')) '\n']);
    % % %         fprintf('Iter %d. Simulation time %.2f hours\n', iter, toc / 3600);
    % % %     end
end
fprintf(['\n' datestr(datetime('now')) '\n']);
fprintf('Simulation time %.2f hours\n', toc / 3600);
%{
%% plot
        
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
%}
% %         SINR_arr1 = squeeze(SINR_arr(:));
% %          [ff, xx] = ecdf(pow2db(SINR_arr1(:,1)));
% %          figure(2); hold on; plot(xx,ff,'LineWidth',1.5)
%%
SINR_arr1 = squeeze(sum(SINR_arr(:,:,:),2)/K);
sinr_mean_sim = mean(pow2db(SINR_arr1),2);
sinr_std_sim  = std(pow2db(SINR_arr1),1,2);
figure(aaa); hold on;

%errorbar(N_arr,sinr_mean_sim,sinr_std_sim,'-','LineWidth',2);
plot(N_arr/N_wrap,sinr_mean_sim,'-s','LineWidth',1.5);
set(gca,'XScale','log');
set(gca, 'FontName', 'Times')
box on;
xlim([N_arr(1) N_arr(end)]); grid on;
xlabel('N'); ylabel('SINR (dB)');
10*log10(MSE/K/sum(N_arr)/nIter);

if ber == 1
bert_arr1 = squeeze(sum(bert_arr(:,:,:),2)/K);
bert_mean_sim = mean(bert_arr1,2);
figure(aaa*10); hold on;

%errorbar(N_arr,sinr_mean_sim,sinr_std_sim,'-','LineWidth',2);
plot(N_arr/N_wrap,bert_mean_sim,'-s','LineWidth',1.5);
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca, 'FontName', 'Times')
box on;
xlim([N_arr(1) N_arr(end)]); grid on;
xlabel('SNR (dB)'); ylabel('BER');
end
% % 10*log10(MSE/K/sum(N_arr)/nIter)
tgprintf('Simulation Finished')