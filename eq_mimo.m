function Y = eq_mimo(X, H, SNR_dB, method)
% MIMO equalization
% size of X: NFreq*NTime x NR
% size of H: NR x NT x NFreq
% size of Y: NFreq*NTime x NT

NR    = size(H, 1);
NT    = size(H, 2);
NFreq = size(H, 3);
NTime = size(X,1)/NFreq;

X = reshape(X, [NFreq, NTime, NR]);
X = permute(X, [3 2 1]);

Y = zeros(NT, NTime, NFreq);

method = upper(method);
switch method
    case 'MF'
        for iFreq = 1:NFreq
            H_temp = H(:, :, iFreq);
            X_temp = X(:, :, iFreq);
            D = diag(H_temp' * H_temp);
            D = diag(1 ./ D);
            Y(:, :, iFreq) = D * H_temp' * X_temp;
        end
        
    case 'ZF'
        for iFreq = 1:NFreq
            H_temp = H(:, :, iFreq);
            X_temp = X(:, :, iFreq);
            test_temp=(H_temp' * H_temp);
            Y(:, :, iFreq) = (H_temp' * H_temp) \ H_temp' * X_temp;
        end
        
    case 'MFDC'
        for iFreq = 1:NFreq
% % %             H_temp = H(:, :, iFreq) ;   %%%Works with Co-located
% % %             X_temp = X(:, :, iFreq);
% % %             Y(:, :, iFreq) = H_temp' * X_temp /NR ;
            H_temp = H(:, :, iFreq) ; 
            H_temp = H_temp ./(abs(H_temp).^2);
            X_temp = X(:, :, iFreq);
            Y(:, :, iFreq) = H_temp' * X_temp/NR ;

        end
    case 'MMSE'
        SNR = db2pow(SNR_dB);
        for iFreq = 1:NFreq
            H_temp = H(:, :, iFreq);
            X_temp = X(:, :, iFreq);
            Y(:, :, iFreq) = (H_temp' * H_temp + 1/SNR * eye(NT)) \ H_temp' * X_temp;
        end
            
            otherwise
                error('Unexpected input value (method)');
        end
        
        Y = permute(Y, [3 2 1]);
        Y = reshape(Y, [NFreq*NTime, NT]);
        
        
        
        
        
        