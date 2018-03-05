% MU_PFS
alapha = 0.1;
R = [R1; R2; R3; R4];
uepair = [[uepair1 zeros(size(uepair1,1), 3)]; [uepair2 zeros(size(uepair2,1), 2)]; [uepair3 zeros(size(uepair3,1), 1)]; uepair4];
T = zeros(1,size(R,2)); %   comment

if drop == 1
    [~,ind] = max(sum(R,2));
    T = R(ind,:);
else
    [~,ind] = max(sum(R./(repmat(T, size(R,1), []) + 1e-10), 2));
    T = (1-alapha)*T + alapha*R(ind,:);
end



for i = 1:size(R,1)
    if i == 1
        [~,ind] = max(sum(R,2));
        T = R(ind,:);
    else
        
    end
end