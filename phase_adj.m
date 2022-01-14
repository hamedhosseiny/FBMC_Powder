function [S, res, res2] = phase_adj(X,M,nSymb,OL,act_subc)

Nact = numel(act_subc);
X    = reshape(X,Nact,[]);

% Phase adjustment
% % phase_mat = -1i*ones(M,size(X,2));
% % phase_mat(1:2:end,1:2:end) = 1;
% % phase_mat(2:2:end,2:2:end) = 1;
% % phase_mat = phase_mat(act_subc,:);

phase_mat = -1i*ones(M,size(X,2));
for ii =1:size(phase_mat,1)
    for jj=1:size(phase_mat,2)
        if rem(ii,4) == 1 || rem(ii,4) == 2
            phase_mat(ii, jj) = phase_mat(ii, jj)^(ii+jj)  ;
        else
            phase_mat(ii, jj) = phase_mat(ii, jj)^(ii+jj)   ;
        end
    end
end
phase_mat = phase_mat .* 1;
 phase_mat = phase_mat(act_subc,:);
S = X .* phase_mat;

% Remove the transient
delay = 2*(OL-1);
%res = S(2,1:30);
%res2 = S(4,1:30);
S = S(:, delay + (1:2*nSymb)); 
S = S(:);