
function At = TD_est_kongmimoV2(pilots_pos,pilots,h,n_sub,Nt, Lvector,act_subc)
At = [];
for u_index = 1:Nt
npilots = numel(pilots_pos(u_index,:));
zeta = zeros(npilots,sum(Lvector));

%%%%%%%%%Matrix of effect of the user
L=Lvector(u_index);
Ap=zeros(n_sub,L);
Apilots2 = zeros(1,n_sub);
Apilots2(act_subc(pilots_pos(u_index,:)).')=pilots(1:npilots,u_index);

for m=act_subc(pilots_pos(u_index,:)).'-1 
    for l=0:L-1
        for K=-(size(h,1)-1)/2:(size(h,1)-1)/2
            for p=act_subc(pilots_pos(u_index,:)).'-1
                if (K-l+1+(size(h,1)-1)/2)>0
                    phi = Apilots2(p+1)*h(K-l+1+(size(h,1)-1)/2)...
                        *h(K+1+(size(h,1)-1)/2)*exp(1i*2*pi*(p-m)*K/n_sub)* ...
                        exp(1i*pi*(p-m)/2)*exp(-1i*2*pi*p*l/n_sub); %%(1i^(mod(p-m,2))) phasea adjustment related part
                    Ap(m+1,l+1) = Ap(m+1,l+1)+phi;
                    
                end
            end
        end
    end
end
Ap1 = Ap(act_subc(pilots_pos(u_index,:)).', :);
zeta(:,sum(Lvector(1:u_index-1))+1:sum(Lvector(1:u_index-1))+L) = Ap(act_subc(pilots_pos(u_index,:)), :);
%%%%%%%%for other users in the cell
for nt = 1:Nt
    if (nt ~= u_index)
        npilotsuc = numel(pilots_pos(nt,:));
        L=Lvector(nt);
        pilot_pos2 = act_subc(pilots_pos(nt,:).');
        Apilots2 = zeros(1,n_sub);
        Apilots2(pilot_pos2)=pilots(1:npilotsuc,nt);
        Ap=zeros(n_sub,L);
        for m= 0:n_sub-1   %pilots_pos(u_index,:)-1
            for l=0:L-1
                for K=-(size(h,1)-1)/2:(size(h,1)-1)/2
                    for p= 0:n_sub-1  %pilot_pos2-1
                        if (K-l+1+(size(h,1)-1)/2)>0
                    phi = Apilots2(p+1)*h(K-l+1+(size(h,1)-1)/2)...
                        *h(K+1+(size(h,1)-1)/2)*exp(1i*2*pi*(p-m)*K/n_sub)* ...
                        exp(1i*pi*(p-m)/2)*exp(-1i*2*pi*p*l/n_sub);
                    Ap(m+1,l+1) = Ap(m+1,l+1)+phi;
                            
                        end
                    end
                end
            end
        end
        zeta(:,sum(Lvector(1:nt-1))+1:sum(Lvector(1:nt-1))+L) = ...
            (Ap(act_subc(pilots_pos(u_index,:)),:));
    end
end

%Noise Correlation
V=zeros(n_sub,n_sub);
At = [ At; zeta];
end

end