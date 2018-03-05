% Proportionl Fair MATLAB CODE
i = 1;
snr = 10.^([5 10 20]/10);
R = 100e3*log2(1+snr);
ic = 10;
Tcur = zeros(1,3);
Rcur = zeros(1,3);
while i < 1e4
    R1(i) = R(floor(unifrnd(1,4)));
    R2(i) = R(floor(unifrnd(1,4)));
    R3(i) = R(floor(unifrnd(1,4)));
    Rcur(i,:) = [R1(i) R2(i) R3(i)];
    if i == 1
        [a,b] = max(Rcur(i,:));
        Tcur(i+1,b) = Rcur(i,b);
    else
        Rcur_chk = Rcur(i,:)./(Tcur(i,:)+1e-10);
        [a,b] = max(Rcur_chk);
        Tcur(i+1,:) = (1-1/ic)*Tcur(i,:);
        Tcur(i+1,b) = (1-1/ic)*Tcur(i,b) + (1/ic)*Rcur(i,b);
    end
    Rgot(i,b) = Rcur(i,b);
    i = i+1;
end
Rmean = [sum(Rgot(:,1))/i sum(Rgot(:,2))/i sum(Rgot(:,3))/i]
Rmeansys = sum(Rmean)
