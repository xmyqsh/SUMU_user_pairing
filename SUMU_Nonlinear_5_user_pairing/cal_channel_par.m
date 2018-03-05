% load channel_10s.mat
k = 1;
m=k+(0:RBnumPerSB*syspar.Nr-1);
Hq=cell(numUsers,1);
for u=1:numUsers
    Htmp(1:RBnumPerSB*syspar.Nr,:) = Hf(m,:,u);
    HRBtmp1 = reshape(Htmp.',[],RBnumPerSB);
    HRBtmp2 = reshape(HRBtmp1,syspar.Nt,syspar.Nr,[]);
    HRB{u} = permute(HRBtmp2,[2 1 3]);
end

size(HRB{u}(:,:,1))

for u=1:numUsers
    user_cov(u) = cond(HRB{u}(:,:,1));
    user_amp(u) = sqrt(sum(sum(abs(HRB{u}(:,:,1)).^2)))/4;
end

%%%%%%%%% Fig 4.1 %%%%%%%%%%%%%%%%%%
% cond1 = [1.9487	2.8936	2.012	8.6867	2.2653	2.266	15.0987	3.4174	12.3873	2.0274	4.2443	5.547	2.7607	3.0884	2.511	2.6113	4.2734	3.8443	3.2919	3.2634];
out = cdff(cond1);
plot(out.x, out.y, 'b-o'), grid on
xlabel('channel condition number')
ylabel('F(cond)')
title ('channel condition number CDF')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Fig 4.2 %%%%%%%%%%%%%%%%%%
gain1 = [6.3526	1.1125	0.9081	1.3354	8.7131	10.1584	2.0606	1.032	0.8662	1.0162	1.0833	1.0522	10.9846	4.0274	0.9431	4.1648	5.0146	0.8283	1.5045	0.8354];
out = cdff(gain1);
plot(out.x, out.y, 'b-o'), grid on
xlabel('channel average gain')
ylabel('F(gain)')
title ('channel average gain CDF')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Fig 4.3 %%%%%%%%%%%%%%%%%%
out = cdff(C1);
plot(out.x, out.y, 'b-o'), grid on
xlabel('channel capacity')
ylabel('F(capacity)')
title ('channel capacity CDF')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Fig 4.4 %%%%%%%%%%%%%%%%%%
out1 = cdff(C1);
out2 = cdff(C2);
out3 = cdff(C3);
out4 = cdff(C4);
out5 = cdff(C5);
out6 = cdff(C6);
out7 = cdff(C7);

plot(out1.x, out1.y, 'b-o'), grid on, hold on
plot(out2.x, out2.y, 'r-o')
plot(out3.x, out3.y, 'y-')
plot(out4.x, out4.y, 'k-')
plot(out5.x, out5.y, 'g-')
plot(out6.x, out6.y, 'b-')
plot(out7.x, out7.y, 'r-')
xlabel('paired channel capacity')
ylabel('F(capacity)')
title ('paired channel capacity CDF')
legend('1 user', '2 users', '3 users', '4 users', '5 users', '6 users', '7 users')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Fig 4.5 %%%%%%%%%%%%%%%%%%
out1 = cdff(C1);
out2 = cdff(C2./2);
out3 = cdff(C3./3);
out4 = cdff(C4./4);
out5 = cdff(C5./5);
out6 = cdff(C6./6);
out7 = cdff(C7./7);

plot(out1.x, out1.y, 'b-o'), grid on, hold on
plot(out2.x, out2.y, 'r-o')
plot(out3.x, out3.y, 'y-')
plot(out4.x, out4.y, 'k-')
plot(out5.x, out5.y, 'g-')
plot(out6.x, out6.y, 'b-')
plot(out7.x, out7.y, 'r-')
xlabel('average paired channel capacity per user')
ylabel('F(capacity)')
title ('average paired channel capacity per user CDF')
legend('1 user', '2 users', '3 users', '4 users', '5 users', '6 users', '7 users')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Fig 4.6 %%%%%%%%%%%%%%%%%%
out1 = cdff(sum(num_rank1,2));
out2 = cdff(sum(num_rank2,2));
out3 = cdff(sum(num_rank3,2));
out4 = cdff(sum(num_rank4,2));
out5 = cdff(sum(num_rank5,2));
out6 = cdff(sum(num_rank6,2));
out7 = cdff(sum(num_rank7,2));

plot(out1.x, out1.y, 'b-o'), grid on, hold on
plot(out2.x, out2.y, 'r-o')
plot(out3.x, out3.y, 'y-')
plot(out4.x, out4.y, 'k-')
plot(out5.x, out5.y, 'g-')
plot(out6.x, out6.y, 'b-')
plot(out7.x, out7.y, 'r-')
xlabel('paired channel rank')
ylabel('F(rank)')
title ('paired channel rank CDF')
legend('1 user', '2 users', '3 users', '4 users', '5 users', '6 users', '7 users')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Fig 4.7 %%%%%%%%%%%%%%%%%%
out1 = cdff(reshape(10*log10(SINR1), 1, []));
out2 = cdff(reshape(10*log10(SINR2), 1, []));
out3 = cdff(reshape(10*log10(SINR3), 1, []));
out4 = cdff(reshape(10*log10(SINR4), 1, []));
out5 = cdff(reshape(10*log10(SINR5), 1, []));
out6 = cdff(reshape(10*log10(SINR6), 1, []));
out7 = cdff(reshape(10*log10(SINR7), 1, []));

plot(out1.x, out1.y, 'b-o'), grid on, hold on
plot(out2.x, out2.y, 'r-o')
plot(out3.x, out3.y, 'y-')
plot(out4.x, out4.y, 'k-')
plot(out5.x, out5.y, 'g-')
plot(out6.x, out6.y, 'b-')
plot(out7.x, out7.y, 'r-')
xlabel('sub-channel SINR(dB)')
ylabel('F(SINR)')
title ('sub-channel SINR CDF')
legend('1 user', '2 users', '3 users', '4 users', '5 users', '6 users', '7 users')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% 20 users %%%%%%%%%%%%%%%%%%%%
user_cov_paired1 = cond(HRB{6}(:,:,1))
user_cov_paired2 = cond([HRB{6}(:,:,1); HRB{13}(:,:,1)])
user_cov_paired3 = cond([HRB{5}(:,:,1); HRB{6}(:,:,1); HRB{13}(:,:,1)])
user_cov_paired4 = cond([HRB{5}(:,:,1); HRB{6}(:,:,1); HRB{13}(:,:,1); HRB{17}(:,:,1)])
user_cov_paired5 = cond([HRB{1}(:,:,1); HRB{5}(:,:,1); HRB{6}(:,:,1); HRB{13}(:,:,1); HRB{17}(:,:,1)])
user_cov_paired6 = cond([HRB{1}(:,:,1); HRB{5}(:,:,1); HRB{6}(:,:,1); HRB{13}(:,:,1); HRB{17}(:,:,1); HRB{18}(:,:,1)])
user_cov_paired7 = cond([HRB{1}(:,:,1); HRB{5}(:,:,1); HRB{6}(:,:,1); HRB{13}(:,:,1); HRB{17}(:,:,1); HRB{18}(:,:,1); HRB{19}(:,:,1)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% 10 users %%%%%%%%%%%%%%%%%%%%%
user_cov_paired1 = cond(HRB{5}(:,:,1))
user_cov_paired2 = cond([HRB{1}(:,:,1); HRB{5}(:,:,1)])
user_cov_paired3 = cond([HRB{2}(:,:,1); HRB{3}(:,:,1); HRB{5}(:,:,1)])
user_cov_paired4 = cond([HRB{2}(:,:,1); HRB{3}(:,:,1); HRB{5}(:,:,1); HRB{6}(:,:,1)])
user_cov_paired5 = cond([HRB{2}(:,:,1); HRB{3}(:,:,1); HRB{5}(:,:,1); HRB{6}(:,:,1); HRB{8}(:,:,1)])
user_cov_paired6 = cond([HRB{2}(:,:,1); HRB{3}(:,:,1); HRB{5}(:,:,1); HRB{6}(:,:,1); HRB{8}(:,:,1); HRB{9}(:,:,1)])
user_cov_paired7 = cond([HRB{2}(:,:,1); HRB{3}(:,:,1); HRB{5}(:,:,1); HRB{6}(:,:,1); HRB{8}(:,:,1); HRB{9}(:,:,1); HRB{10}(:,:,1)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load result_10_users_AMC.mat
TargetSINR = [-7.5176,-5.6776,-3.5476,-1.5076,0.4924,2.3924,4.3972,6.2272,8.1972,10.1191,12.0191,13.9691,15.8191,17.7491,19.6391];
for k = 1:length(SINR1user)
    [~, MCS1user(k,1)] = min(abs(SINR1user(k,1) - TargetSINR));
    [~, MCS2user(k,1)] = min(abs(SINR2user(k,1) - TargetSINR));
    [~, MCS2user(k,2)] = min(abs(SINR2user(k,2) - TargetSINR));
    [~, MCS3user(k,1)] = min(abs(SINR3user(k,1) - TargetSINR));
    [~, MCS3user(k,2)] = min(abs(SINR3user(k,2) - TargetSINR));
    [~, MCS3user(k,3)] = min(abs(SINR3user(k,3) - TargetSINR));
end
MCS1user = MCS1user - 1;
MCS2user = MCS2user - 1;
MCS3user = MCS3user - 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m = 1:length(H3user)
    cond_H3user(m) = cond(H3user(:,:,m));
end
stem(cond_H3user)
out = cdff(cond_H3user);
figure, plot(out.x, out.y, 'r-o')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m = 1:length(H4user)
    cond_H4user(m) = cond(H4user(:,:,m));
end
% stem(cond_H4user)
% out = cdff(cond_H4user);
% figure, plot(out.x, out.y, 'r-o')
sel_ind = find(cond_H4user < 10)
uepair=combntns(1:30,4);
sel_cha = uepair(sel_ind,:);
unique(sel_cha)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for u = 1:30
    cal(u) = sum(sum(selected_uepair == u));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
