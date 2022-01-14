function h = channel_estiation(y,pilots_pos, At, M, Nt, L, OL, act_subc, Vtt2)

pilottot = zeros(Nt, numel(pilots_pos(1,:)));
s_I2 = phase_adj(y,M,1,OL,act_subc);
for nt =1:Nt
% %     y_I = reshape(y, numel(act_subc), []);
% %     y_I = y_I(:,1:2:end);   % Taking the Inphase (includes piltos)
% %     s_I=[];
% %     % Apply Phase Adjustment
% %     for kp = 1:numel(act_subc) %or M/2
% %         s_I = [s_I;((-1i)^(act_subc(kp)-1)*y_I(kp,:))];
% %     end
    
    % Arrange Pilots Vector
    pilottot(nt,:) =  (s_I2(pilots_pos(nt,:),1)); 
end
pilottot = pilottot.';
pilottot = pilottot(:);

%LS estimate
  correlatedestimate = (At'*At)\At'* pilottot;
%%%Noise correlation corrected
% %         Vtt2 = Vtt * sigma2;
%           correlatedestimate = (At'*inv(Vtt2)*At)\At'* inv(Vtt2) *pilottot;

h = reshape(correlatedestimate, L , Nt);
