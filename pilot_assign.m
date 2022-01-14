function [pilot_pos, pilots] = pilot_assign(Nt, n_sub, pilot_set,act_subc)
pilots = zeros(n_sub, Nt);
%%%%%%%%%%Assign Pilots' position
switch Nt
    case 1
        pilot_pos = (1:1:n_sub);
    case 2
        pilot_pos = (1:2:n_sub);
        pilot_pos(2,:) = (2:2:n_sub);
    case 4
        pilot_pos = (1:4:n_sub);
        pilot_pos(2,:) = (2:4:n_sub);
        pilot_pos(3,:) = (3:4:n_sub);
        pilot_pos(4,:) = (4:4:n_sub);
    case 8
        pilot_pos = (1:8:n_sub);
        pilot_pos(2,:) = (2:8:n_sub);
        pilot_pos(3,:) = (3:8:n_sub);
        pilot_pos(4,:) = (4:8:n_sub);
        pilot_pos(5,:) = (5:8:n_sub);
        pilot_pos(6,:) = (6:8:n_sub);
        pilot_pos(7,:) = (7:8:n_sub);
        pilot_pos(8,:) = (8:8:n_sub);
    otherwise
        error('Undefined Variable')
end

for ik = 1:Nt
    pilots((pilot_pos(ik,:)), ik) = pilot_set(:,ik);
end
