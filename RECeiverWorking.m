% close all;
clear all;

%% Demodulation

%Ayaz Test
 load('Parameters_test.mat');
 Y = read_complex_binary('received.dat',100000);

%Local Transmitted test
% load('receivedSignal');
% load('Parameters.mat');
% Y = read_complex_binary('transmit.dat',10000);
% Y = [zeros(1,1); Y];
% % %%Timing phase acqusition
% % for i =1:L-1
% %     p(i) = xBB(i:L:end)' * xBB(i:L:end);
% %     
% % end
% % [mx, idx] = max(p);
% % % figure, plot((1:L-1),p);xlabel('Sampling shift')
% % % ylabel('Power')
% % % fname = sprintf('timingphaseXRF3test.eps');
% % % saveas(gcf,fname,'eps')
% % 
% % xBBd = xBB(idx:L:end);
% % 
% % %%%%Preamble detection
% % itn=1;
% % Np=32;
% % for iii =2*Np+2:length(xBBd)
% %     ryy(itn) = sum (xBBd(iii-2:-1:iii-Np-1) .* conj ( xBBd(iii-Np-2:-1:iii-2*Np-1)));
% %     itn=itn+1;
% % end


%%%%OFDM SISO Method
%% Correlate for LTS
% Complex cross correlation of Rx waveform with time-domain LTS
a = 1;
SS = zeros(Nact*1,1);
preamble = SMT_PPN_mod_pilot_mu_rec(SS, pulse_shape, 1, M, act_subc,K,pilots,ones(1,K),1); %Transmit Signal
preamble = preamble(1:OL*M);
unos = ones(size(preamble.'))';
v0 = filter(flipud(preamble'),a,Y);
v1 = filter(unos,a,abs(Y).^2);
m_filt = (abs(v0).^2)./v1; % normalized correlation
lts_corr = m_filt;
[rho_max, ipos] = max(lts_corr);

lts_corr = abs(conv(Y,flipud(preamble)));
inds = (1:length(lts_corr)); ipos = mean(inds(lts_corr > 12));
ppos = 1 + ipos - length(preamble)/2 * 5
Y2=Y(27);
% L=6;

% At = TD_est_kongmimoV2(pilot_pos, pilot_set, h2, M, K, ones(K,1)*L,(1:M).' );   %Generate estimation matrix
% for ipos = 250+12+1:95000
% ipos = ipos - M * 
Y = Y2(1:end);
  Y= Y .* exp( 1i * 0.5);
%     Y= Y .* exp(- 1i * 0.5);
% %
% % payload_ind = ipos +1;
% % lts_ind = payload_ind - (Ng*2)*M;

% % % Re-extract LTS for channel estimate
% % rx_lts = rx_vec_iris(lts_ind : lts_ind+159);
% % rx_lts1 = rx_lts(-64+-FFT_OFFSET + [97:160]);
% % rx_lts2 = rx_lts(-FFT_OFFSET + [97:160]);


%% %%%%%%%%%%%


% Here, we only perform the analysis filter bank and the decimation for
% each antenna chain. We DO NOT apply any equalization as well as the
% phase compensation here.
%    mod_pulse = mod_pulse_estimator(pulse_shape,pdp_arr,[]);
%mo_pulse = modified_prototype(pulse_shape,pdp_arr(:,1));
% % R1= zeros(7668,N);
% % RP = zeros(64*213,N);
% nSymb = 100;
for inn = 1:N
    [R1(:,inn), RP(:,inn)] = SMT_PPN_demod2(Y(:,inn), pulse_shape, nSymb, M, OL, act_subc);
end
%clear Y;
%%%Time Domain channel estimation
h2 = [pulse_shape; pulse_shape(end)];
L=20;
At = TD_est_kongmimoV2(pilot_pos, pilot_set(:,1), h2, M, K, ones(K,1)*L,(1:M).' );   %Generate estimation matrix

h_est = zeros(L,N,K);
fl_delay  = 2*(OL-1);    %transient time
for inn=1:N
    h_est(:,inn,:) = channel_estiation(RP(1:(fl_delay+2)*M,inn), pilot_pos, At, M, K, L, OL, (1:M).');
end
% Channel MSE Calculation
MSE=0;
% % for ikk=1:K
% %     for innn=1:N
% %         MSE = var(squeeze(h(:,innn,ikk))-squeeze(h_est(:,innn,ikk)))/var(squeeze(h(:,innn,ikk))) + MSE;
% %     end
% % end

%    h_est = h;
hfreq = fft(h_est,M,1);
hfreq = hfreq(act_subc,:,:);
hfreq = permute(hfreq,[2 3 1]);
N_arr = N;
for in = 1:length(N_arr)
    
    N   = N_arr(in);  % Number of BS antennas
    % MIMO equalization using channel frequency gains (MRC, ZF, MMSE)
    SNR_db =10;
    eq_mode = 'MF';
    R2 = eq_mimo(R1(:,1:N), hfreq(1:N,:,:), SNR_db, eq_mode);
    %R2 = squeeze(R1(:,1,1));
    % Consider terminal k
    delay = 11;
    for ik = 1:K
        R2_k = R2(:,ik);        
        R3_k = R2_k; %./ power_cf(ik); 

        
                %%%%%%%%%%DDCR 
        s_hat = R3_k; %(1:Nact:end);
phi=zeros(size(s_hat)); s1=zeros(size(s_hat)); mu=.01;
for n=1:length(s_hat)-1
    s1(n)=s_hat(n)*exp(-1i*phi(n));
    s2=sign(real(s1(n)))+1i*sign(imag(s1(n)));
    s12=s1(n)*s2'; e=imag(s12)/real(s12);
    phi(n+1)=phi(n)+mu*e;
end

  phi = 0.5;
s_hat = s_hat .* exp(-1i * phi);
%              R3_k = s_hat;
                R3_k = persubcarrier_eq_mu_cell(R3_k,M,squeeze(h_est(:,1:N,:)),...
        act_subc,pulse_shape, hfreq(1:N,:,:),K,ik,delay,...
        noise_std*0,eq_mode,0*noise_std^2 / M,1,1);
        R4_k = phase_adj(R3_k,M,nSymb,OL,act_subc);
        R5_k = OQAMpostprocessing(R4_k, Nact, nSymb);
        
        

        %MSE_arr(in,ik) = mean(abs((R5_k) - (S(:,ik))).^2); %(1:Nact:end)
        MSE_arr(in,ik) = mean(abs(R5_k(Ng * Nact * Np+1:end) - S(Ng * Np* Nact+1:end,ik)).^2);
        
        %%
        %%%%BER
                    dataSymbolsOut_n = qamdemod(R5_k,nQAM,'bin');
                    dataOutMatrix_n = de2bi(dataSymbolsOut_n,log2(nQAM));
                    dataOut_n = dataOutMatrix_n(:);                   % Return data in column vector
                    bin2file(dataOut_n,'part_IV.txt');

        % %
                            dataSymbolsOut_n = qamdemod(R5_k(Nact*Ng+1:end),nQAM,'bin');
                            dataOutMatrix_n = de2bi(dataSymbolsOut_n,log2(nQAM));
                            dataOut_n = dataOutMatrix_n(:);                   % Return data in column vector
                
                            dataSymbolsOut_n = qamdemod(S(Nact*Ng+1:end,ik),nQAM,'bin');
                            dataOutMatrix_n = de2bi(dataSymbolsOut_n,log2(nQAM));
                            dataIn_n = dataOutMatrix_n(:);                   % Return data in column vector
                
                            [numErrors_n,berrate_n] = biterr(dataIn_n,dataOut_n);
% %                             fprintf('\nEstimation bit error rate = %5.2e, based on %d errors\n', ...
% %                                 berrate_n,numErrors_n);
                            bert(1) = berrate_n
        
        
    end
    
end
%clear S R1 R2 hfreq ;

SINR_sim = 1./MSE_arr;
%SINR_arr = cat(3,SINR_arr,SINR_sim);
SINR_arr(:,:) = SINR_sim;
%%%%%%%%%%%%%%%%%%%%%%
% Plot the results
%%%%%%%%%%%%%%%%%%%%%%

% % fprintf(['\n' datestr(datetime('now')) '\n']);
% % fprintf('Simulation time %.2f hours\n', toc / 3600);
SINR_arr1 = squeeze(sum(SINR_arr(:,:,:),2)/K);
sinr_mean_sim = mean(pow2db(SINR_arr1),2) %SINR at the Receiver
% % sinr_std_sim  = std(pow2db(SINR_arr1),1,2);
% % %figure(3); hold on;
% % %errorbar(N_arr,sinr_mean_sim,sinr_std_sim,'-','LineWidth',2);
% % % % plot(N_arr,sinr_mean_sim,'-s','LineWidth',1.5);
% % % % set(gca,'XScale','log');
% % % % set(gca, 'FontName', 'Times')
% % % % box on;
% % % % xlim([N_arr(1) N_arr(end)]); grid on;
% % % % xlabel('N'); ylabel('SINR (dB)');
% % % % 10*log10(MSE/K/sum(N_arr))
% % sinrd(ipos) = sinr_mean_sim;
% % clear R1 RP;
% end