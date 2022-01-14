function y = QAM_demod(r, M, noiseVar)
% r: vector of received QAM symbols
% M: Size of signal constellation, e.g., 2, 4, 16, ...
% noiseVar: noise variance in case of soft decision decoding

persistent QAMdemod;

if nargin < 3 % Hard decision
    if isempty(QAMdemod)
    QAMdemod = comm.RectangularQAMDemodulator('ModulationOrder', M, ...
        'BitOutput', true, ...
        'SymbolMapping', 'Gray', ...
        'NormalizationMethod', 'Average power', ...
        'DecisionMethod', 'Hard decision');       
    end
    y = step(QAMdemod, r);
    
else % Soft decision
    if isempty(QAMdemod)
    QAMdemod = comm.RectangularQAMDemodulator('ModulationOrder', M, ...
        'BitOutput', true, ...
        'SymbolMapping', 'Gray', ...
        'NormalizationMethod', 'Average power', ...
        'DecisionMethod', 'Approximate log-likelihood ratio', ...
        'VarianceSource', 'Input port');   
    end
    y = step(QAMdemod, r, noiseVar);
end


