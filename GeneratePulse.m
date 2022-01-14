function P = GeneratePulse(L, K, Method, alpha)

% L: Number of samples in each symbol period
% K: Length of filter in terms of symbol period

Method = upper(Method);
switch Method
    
        case 'SRRC'        
         P = sr_Nyquist_p(L*K, L, alpha,1);

       % P(end) = [];    % Even Length! 
        % P = [P; zeros(L-1,1)]; % Even Length

    case 'RC'        

        P = sr_cos_p(L*K, L, alpha);
       % P(end) = [];    % Even Length! 
        % P = [P; zeros(L-1,1)]; % Even Length
        
    case 'PHYDYAS'
        switch K
            case 3
                C = [0.41143783 -0.91143783 1 -0.91143783 0.41143783]';
            case 4
                C = [-0.23514695 +0.70710678 -0.97195983 1 ...
                     -0.97195983 +0.70710678 -0.23514695]';
            case 5
                C = [+0.12747868 -0.50105361 +0.86541624 -0.99184131 +1 ...
                     -0.99184131 +0.86541624 -0.50105361 +0.12747868]';
            case 6
                C = [-0.06021021 +0.31711593 -0.70710678 +0.94838678 ...
                     -0.99818572 1 -0.99818572 +0.94838678 -0.70710678 ...
                     +0.31711593 -0.06021021]';
            case 7
                C = [+0.03518546 -0.20678881 0.53649931 -0.84390076 ...
                      0.97838560 -0.99938080 1 -0.99938080 0.97838560 ...
                     -0.84390076 0.53649931 -0.20678881 +0.03518546]';
            case 8
                C = [-0.03671221 +0.18871614 -0.44756522 +0.70710678 ...
                     -0.89425129 +0.98203168 -0.99932588 1 -0.99932588 ...
                     +0.98203168 -0.89425129 +0.70710678 -0.44756522 ...
                     +0.18871614 -0.03671221]';
            otherwise
                error('Undefined Variable')
        end;
        Lc = length(C);
        C(end+1:K*L) = 0;
        C = circshift(C,[(1-Lc)/2 0]);
        P = ifft(C);

    case 'IOTA' % Please check later
        P = EGF(L, L*K, alpha);

    otherwise
        error('Undefined Variable')
end;

P = P./sqrt(P.'*P); % Normalization