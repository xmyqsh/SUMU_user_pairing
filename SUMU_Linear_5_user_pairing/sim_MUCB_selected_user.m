function capacity=sim_MUCB_selected_user(ctrlpar)
global user_selected channel
capacity = [];

if nargin>1, error('Too many input arguments.');end
tic
addpath('SCME');

[syspar,scmpar,antpar]=argdeliv(ctrlpar);
cbpar=struct('AntennaNum',syspar.Nt,...
    'AntennaConfig',strcmp(syspar.ScmOptions,'polarized'),...% choice {0:ULA,1:X-polarized}
    'CorrFlag',strcmp(syspar.antconfig,'diversity-ant'),...% choice {0:closely spaced antennas,1:diversity antennas}
    'UEPolarizeSlant',45);

ebno_range=syspar.ebno_range;%仿真信噪比的范围'ebno_range',10.^(([-5:5:20])/10)
num_drops=syspar.ndrops;%循环次数ctrlpar.ndrops = 100;
numUsers=syspar.nusers;%用户数ctrlpar.nusers = 5;
maxlayer=syspar.maxlayer;%最大层数ctrlpar.maxlayer = 1;
FBgran=syspar.FBgran;%ctrlpar.FBgran = 'SB Level';????
capacity_vector=zeros(length(ebno_range),1);

phase_error = ctrlpar.phase_error;  %   degree
amp_error = ctrlpar.amp_error;    %   dB

%%%%%%%% PFS parameter %%%%
alapha = 0.1;
T = zeros(1,numUsers); %   comment
%%%%%%%%%%%%%%%%

% for drop = 1:num_drops
%     fprintf('drop: %d/%d\n', drop, num_drops);
%     Hf=zeros(syspar.ntimesamples*syspar.Nr*syspar.ntones,syspar.Nt,numUsers);%%%%%不懂？？？？？？？？
%     for ueIdx=1:numUsers
%         linkpar=linkparset;
%         if ueIdx==1
%             d=linkpar.MsBsDistance;
%         else
%             % users are randomly dropped with the same geometry, so they
%             % have same SNR
%             linkpar.MsBsDistance=d;
%         end
%         
%         if drop == 1
%             Frame_time{ueIdx} = 0;
%             [Hf(:,:,ueIdx) chan_par_backup{ueIdx}] =calc_channel(syspar,scmpar,linkpar,antpar,Frame_time{ueIdx});
%         else
%             Frame_time{ueIdx} = Frame_time{ueIdx} + 10;  %1e-3;
%             [Hf(:,:,ueIdx) chan_par_backup{ueIdx}] =calc_channel(syspar,scmpar,linkpar,antpar,Frame_time{ueIdx},chan_par_backup{ueIdx});
%         end
%         
%     end
%     channel(:,:,:,drop) = Hf;
% 
% end
%     save('channel.mat', 'channel')
    
%     load channel_50u_8t_2r.mat;
for drop = 1:num_drops   %[1 2 10 100 1000]  %[1]
    Hf = channel(:,:,:,drop);
    %% FBgran == 'SB Level'
    RBFreqSB = 1 ;
    TTIgran  = 1 ;
    FBdelay  = 0 ;
    RBnumPerSB = RBFreqSB*TTIgran;
    for k=1:(RBnumPerSB*syspar.Nr):size(Hf,1) % 重复执行的次数与RB的数目相同
        m=k+(0:RBnumPerSB*syspar.Nr-1);
        Hq=cell(numUsers,1);
        for u=1:numUsers
            Htmp(1:RBnumPerSB*syspar.Nr,:) = Hf(m,:,u);
            HRBtmp1 = reshape(Htmp.',[],RBnumPerSB);
            HRBtmp2 = reshape(HRBtmp1,syspar.Nt,syspar.Nr,[]);
            HRB{u} = permute(HRBtmp2,[2 1 3]);
        end
        
%         numUsers = 50;
        HRB_temp = HRB;
        clear HRB;
%         user_selected = [2 3 8 9 10 11 12 15 18 20];
%         user_selected = [1	2	3	5	10	12	13	15	16	17	19	20	21	25	26	27	29	30	31	32	37	38	39	40	41	42	43	44	49	50];
%         
        for u = 1:length(user_selected)
            HRB{u} = HRB_temp{user_selected(u)};
        end
        CMax=zeros(length(ebno_range),1);
        % MU-MIMO scheduling procedure
        % multiplexing of two users
        [C_err1, C_real1, R1, uepair1, SINR1, num_rank1] = user_paired_1(length(user_selected), RBnumPerSB, syspar, HRB, ebno_range, phase_error, amp_error);
        [C_err2, C_real2, R2, uepair2, SINR2, num_rank2] = user_paired_2(length(user_selected), RBnumPerSB, syspar, HRB, ebno_range, phase_error, amp_error);
        [C_err3, C_real3, R3, uepair3, SINR3, num_rank3] = user_paired_3(length(user_selected), RBnumPerSB, syspar, HRB, ebno_range, phase_error, amp_error);
        [C_err4, C_real4, R4, uepair4, SINR4, num_rank4] = user_paired_4(length(user_selected), RBnumPerSB, syspar, HRB, ebno_range, phase_error, amp_error);

        
% add PFS
        % MU_PFS
        R = [R1; R2; R3; R4];
        uepair = [[uepair1 zeros(size(uepair1,1), 3)]; [uepair2 zeros(size(uepair2,1), 2)]; [uepair3 zeros(size(uepair3,1), 1)]; uepair4];
        C_err = [C_err1; C_err2; C_err3; C_err4];
        C_real = [C_real1; C_real2; C_real3; C_real4];

        if drop == 1
            [t,ind] = max(sum(R,2));
            T = R(ind,:);
        else
            [t,ind] = max(sum(R./(repmat(T, size(R,1), []) + 1e-10), 2));
            T = (1-alapha)*T + alapha*R(ind,:);
        end
        R_got(drop,:) = R(ind,:);
        selected_uepair(drop,:) = uepair(ind,:);
        selected_C_err(drop) = C_err(ind);
        selected_C_real(drop) = C_real(ind);
        num_uepair(drop) = sum(selected_uepair(drop,:)~=0);
        
    end
end

eval(['save ','F:\efijlmj\zhongxing\results/result_linear_users_', num2str(length(user_selected)), '_', 'ebno', '_', num2str(10*log10(ebno_range)), 'dB.mat']);
toc
rmpath('SCME');
end

function [r, num_rank] = unique_r(r_1, num_user2paired)
    r_1(17:end) = [];
    r = [r_1(1)];
    num_rank = [];
    count = 1;
    for n = 2:length(r_1)
        if r(end) ~= r_1(n)
            r = [r r_1(n)];
            num_rank = [num_rank count];
            count = 0;
        end
        count = count + 1;
    end
    num_rank = [num_rank count];
    for n = length(r)+1:num_user2paired
        r(n) = 0;
        num_rank(n) = 0;
    end
    
end

function [H_eff, M_BD] = BD_channel(H, numUsers, spatial_freedom)
    num_UA = size(H,1)/numUsers;
    M_BD = [];
    H_eff = [];
    for u_ind = 1:numUsers
        H_temp = H;
        H_temp([(u_ind-1)*num_UA+1:u_ind*num_UA], :) = [];
        [u s v] = svd(H_temp);
        M_BD = [M_BD v(:,end-spatial_freedom(u_ind)+1:end)];
    end
    H_eff = H*M_BD;
end

function [sorted_H, sorted_user] = sort_users(H, numUsers)
    num_UA = size(H,1)/numUsers;
    sorted_H = [];
    sorted_user = [];
    num_user = 1:numUsers;
    while (numUsers-1)
        r = zeros(1,numUsers);
        det_H = det(H*H');
        for u = 1:numUsers
            %         H_u = H((u-1)*num_UA+1:u*num_UA,:);
            H_temp = H;
            H_temp([(u-1)*num_UA+1:u*num_UA], :) = [];
            r(u) = det_H/det(H_temp*H_temp'); % ( det(H*H')/det(H_temp*H_temp') ).^( 1/(2*num_UA) )
        end
        [t,ind] = max(abs(r));
        
        sorted_H = [H((ind-1)*num_UA+1:ind*num_UA,:); sorted_H];
        sorted_user = [num_user(ind), sorted_user];
        
        H([(ind-1)*num_UA+1:ind*num_UA], :) = [];
        num_user(ind) = [];
        numUsers = numUsers-1;
    end
    
    sorted_H = [H; sorted_H];
    sorted_user = [num_user, sorted_user];

end

function [C_err, C_real, R, uepair, SINR, num_rank] = user_paired_1(numUsers, RBnumPerSB, syspar, HRB, ebno_range, phase_error, amp_error)
        num_user2paired = 1;
        spatial_freedom = [8];
        Nr = syspar.Nr;
        R = [];
        C_real = [];
        SINR = [];
        num_rank = [];
        %%%%%%%%%%%%%
        uepair=combntns(1:numUsers,num_user2paired);% different pair of UEs   % uepair=combntns(1:numUsers,numPair);
        for m=1:size(uepair,1) % a specific UE pairing
            u1=uepair(m,1);
%             u2=uepair(m,2);
%             u3=uepair(m,3);
%             u4=uepair(m,4);
            %%%%%%%%%%%%
            H_u1 = HRB{u1};
%             H_u2 = HRB{u2};
%             H_u3 = HRB{u3};
%             H_u4 = HRB{u4};
            for k = 1:RBnumPerSB
%                 H{k} = [H_u1(:,:,k); H_u2(:,:,k); H_u3(:,:,k); H_u4(:,:,k)];
%                 H{k} = [H_u1(:,:,k); H_u2(:,:,k)];
                H{k} = [H_u1(:,:,k)];

                
                                %   add channel error
                real_H = H{k};
                error_phase = pi*5/180*( phase_error(1) + (phase_error(2)-phase_error(1)).*rand(size(real_H,1), size(real_H,2)) );
                error_amp = 10.^( 0.1.*( amp_error(1) + (amp_error(2)-amp_error(1)).*rand(size(real_H,1), size(real_H,2)) ) );
                error_H = real_H.*error_amp.*exp(1i*error_phase);
                delta_H = real_H - error_H;
                
                if cond(error_H) > 10
                    continue;
                end
                
                [H_eff, M_BD] = BD_channel(error_H, num_user2paired, spatial_freedom);
                
                U = [];
                V = [];
                S = [];
                shape = [];
                start = 1;
                for u_ind = 1:num_user2paired
%                     [u s v] = svd( H_eff( (u_ind-1)*Nr+1:u_ind*Nr , (u_ind-1)*Nr+1:u_ind*Nr ) );
                    [u s v] = svd( H_eff( (u_ind-1)*Nr+1:u_ind*Nr , start:sum(spatial_freedom(1:u_ind)) ) );
                    start = start + spatial_freedom(u_ind);
                    r = diag(s);
                    ind = find( diag(s)/sum(diag(s)) < 0.1 );
                    num_rank(m,u_ind) = Nr - length(ind);
                    u(:,ind) = [];
                    v(:,ind) = [];
                    r(ind) = [];
                    U = blkdiag(U, u);
                    V = blkdiag(V, v);
                    S = [S; r];
                    shape = blkdiag(shape, ones(1,num_rank(m,u_ind)));
                end
                Q = U';
                
                signel_power = S.^2/length(S);
                
                delta_error_m = Q*delta_H*M_BD;
                for t = 1:size(delta_error_m,1)
                    interference_power(t) = (norm(delta_error_m(t,:))^2)/sum(num_rank(m,:));
                end
                
                R_temp = zeros(1,numUsers);
                R_temp(1,[uepair(m,:)]) = shape*log2(1+signel_power*ebno_range);
                R = [R; R_temp];
                
                C_real(m) = sum(log2(1+(signel_power*ebno_range)./(1+interference_power'*ebno_range)));

                clear signel_power;
                clear interference_power;
            end    
        end
        C_err = sum(R,2);
        C_real([find(C_real == 0)]) = [];
        C_real = C_real';
        
        
%                 %   add channel error
%                 real_H = H{k};
%                 error_phase = pi*5/180*( phase_error(1) + (phase_error(2)-phase_error(1)).*rand(size(real_H,1), size(real_H,2)) );
%                 error_amp = 10.^( 0.1.*( amp_error(1) + (amp_error(2)-amp_error(1)).*rand(size(real_H,1), size(real_H,2)) ) );
%                 error_H = real_H.*error_amp.*exp(1i*error_phase);
%                 delta_H = real_H - error_H;
%                 
%                 [L,Q,G,PPP]= GMD_encoder(error_H,num_user2paired,syspar.Nr*ones(num_user2paired,1));
%                 [r(m,:), num_rank(m,:)] = unique_r(diag(G), num_user2paired);
%                 SINR(m,:) = (ebno_range.*r(m,:).^2)./sum(num_rank(m,:));
%                 R_temp = zeros(1,numUsers);
%                 R_temp(1,[uepair(m,:)]) = log2(1+SINR(m,:)).*num_rank(m,:);
%                 R = [R; R_temp];
%                 
% %                 if size(Q,2) ~= 2 || size(Q,1) ~= 2 
% %                     pause;
% %                 end
%                 
%                 delta_error_m = Q'*delta_H;
%                 
%                 for t = 1:size(delta_error_m,1)
%                     delta_error(t) = (ebno_range*norm(delta_error_m(t,:))^2)/sum(num_rank(m,:));
%                 end
%                 
%                 signel_power = SINR(m,1)*ones(num_rank(m,1),1);
%                 if num_rank(m,:) ~= 2
%                     delta_error(2) = []; 
%                 end
%                 
%                 C_real(m) = sum(log2(1+(signel_power)./(1+delta_error')));
%             end
%         end
%         C_err = sum(R,2);
%         C_real([find(C_real == 0)]) = [];
%         C_real = C_real';
end

function [C_err, C_real, R, sorted_uepair, SINR, num_rank] = user_paired_2(numUsers, RBnumPerSB, syspar, HRB, ebno_range, phase_error, amp_error)
        num_user2paired = 2;
        spatial_freedom = [4, 4];
        Nr = syspar.Nr;
        R = [];
        C_real = [];
        SINR = [];
        num_rank = [];
        sorted_uepair = [];
        %%%%%%%%%%%%%
        uepair=combntns(1:numUsers,num_user2paired);% different pair of UEs   % uepair=combntns(1:numUsers,numPair);
        for m=1:size(uepair,1) % a specific UE pairing
            u1=uepair(m,1);
            u2=uepair(m,2);
%             u3=uepair(m,3);
%             u4=uepair(m,4);
            %%%%%%%%%%%%
            H_u1 = HRB{u1};
            H_u2 = HRB{u2};
%             H_u3 = HRB{u3};
%             H_u4 = HRB{u4};
            for k = 1:RBnumPerSB
%                 H{k} = [H_u1(:,:,k); H_u2(:,:,k); H_u3(:,:,k); H_u4(:,:,k)];
                H{k} = [H_u1(:,:,k); H_u2(:,:,k)];
                H2user(:,:,m) = H{k};
                
                %   add channel error
                real_H = H{k};
                error_phase = pi*5/180*( phase_error(1) + (phase_error(2)-phase_error(1)).*rand(size(real_H,1), size(real_H,2)) );
                error_amp = 10.^( 0.1.*( amp_error(1) + (amp_error(2)-amp_error(1)).*rand(size(real_H,1), size(real_H,2)) ) );
                error_H = real_H.*error_amp.*exp(1i*error_phase);
                delta_H = real_H - error_H;
                
                if cond(error_H) > 10
                    continue;
                end
                
                [H_eff, M_BD] = BD_channel(error_H, num_user2paired, spatial_freedom);
                
                U = [];
                V = [];
                S = [];
                shape = [];
                start = 1;
                for u_ind = 1:num_user2paired
%                     [u s v] = svd( H_eff( (u_ind-1)*Nr+1:u_ind*Nr , (u_ind-1)*Nr+1:u_ind*Nr ) );
                    [u s v] = svd( H_eff( (u_ind-1)*Nr+1:u_ind*Nr , start:sum(spatial_freedom(1:u_ind)) ) );
                    start = start + spatial_freedom(u_ind);
                    r = diag(s);
                    ind = find( diag(s)/sum(diag(s)) < 0.1 );
                    num_rank(m,u_ind) = Nr - length(ind);
                    u(:,ind) = [];
                    v(:,ind) = [];
                    r(ind) = [];
                    U = blkdiag(U, u);
                    V = blkdiag(V, v);
                    S = [S; r];
                    shape = blkdiag(shape, ones(1,num_rank(m,u_ind)));
                end
                Q = U';
                
                signel_power = S.^2/length(S);
                
                delta_error_m = Q*delta_H*M_BD;
                for t = 1:size(delta_error_m,1)
                    interference_power(t) = (norm(delta_error_m(t,:))^2)/sum(num_rank(m,:));
                end
                
                R_temp = zeros(1,numUsers);
                R_temp(1,[uepair(m,:)]) = shape*log2(1+signel_power*ebno_range);
                R = [R; R_temp];
                sorted_uepair = [sorted_uepair; uepair(m,:)];
                
                C_real(m) = sum(log2(1+(signel_power*ebno_range)./(1+interference_power'*ebno_range)));

                clear signel_power;
                clear interference_power;
            end    
        end
        C_err = sum(R,2);
        C_real([find(C_real == 0)]) = [];
        C_real = C_real';
end

function [C_err, C_real, R, sorted_uepair, SINR, num_rank] = user_paired_3(numUsers, RBnumPerSB, syspar, HRB, ebno_range, phase_error, amp_error)
        num_user2paired = 3;
        spatial_freedom = [2 3 3];
        Nr = syspar.Nr;
        R = [];
        C_real = [];
        SINR = [];
        num_rank = [];
        sorted_uepair = [];
        %%%%%%%%%%%%%
        uepair=combntns(1:numUsers,num_user2paired);% different pair of UEs   % uepair=combntns(1:numUsers,numPair);
        for m=1:size(uepair,1) % a specific UE pairing
            u1=uepair(m,1);
            u2=uepair(m,2);
            u3=uepair(m,3);
%             u4=uepair(m,4);
            %%%%%%%%%%%%
            H_u1 = HRB{u1};
            H_u2 = HRB{u2};
            H_u3 = HRB{u3};
%             H_u4 = HRB{u4};
            for k = 1:RBnumPerSB
                H{k} = [H_u1(:,:,k); H_u2(:,:,k); H_u3(:,:,k)];
                H3user(:,:,m) = H{k};
                
                %   add channel error
                real_H = H{k};
                error_phase = pi*5/180*( phase_error(1) + (phase_error(2)-phase_error(1)).*rand(size(real_H,1), size(real_H,2)) );
                error_amp = 10.^( 0.1.*( amp_error(1) + (amp_error(2)-amp_error(1)).*rand(size(real_H,1), size(real_H,2)) ) );
                error_H = real_H.*error_amp.*exp(1i*error_phase);
                delta_H = real_H - error_H;
                
                if cond(error_H) > 10
                    continue;
                end
                
                [H_eff, M_BD] = BD_channel(error_H, num_user2paired, spatial_freedom);
                
                U = [];
                V = [];
                S = [];
                shape = [];
                start = 1;
                for u_ind = 1:num_user2paired
%                     [u s v] = svd( H_eff( (u_ind-1)*Nr+1:u_ind*Nr , (u_ind-1)*Nr+1:u_ind*Nr ) );
                    [u s v] = svd( H_eff( (u_ind-1)*Nr+1:u_ind*Nr , start:sum(spatial_freedom(1:u_ind)) ) );
                    start = start + spatial_freedom(u_ind);
                    r = diag(s);
                    ind = find( diag(s)/sum(diag(s)) < 0.1 );
                    num_rank(m,u_ind) = Nr - length(ind);
                    u(:,ind) = [];
                    v(:,ind) = [];
                    r(ind) = [];
                    U = blkdiag(U, u);
                    V = blkdiag(V, v);
                    S = [S; r];
                    shape = blkdiag(shape, ones(1,num_rank(m,u_ind)));
                end
                Q = U';
                
                signel_power = S.^2/length(S);
                
                delta_error_m = Q*delta_H*M_BD;
                for t = 1:size(delta_error_m,1)
                    interference_power(t) = (norm(delta_error_m(t,:))^2)/sum(num_rank(m,:));
                end
                
                R_temp = zeros(1,numUsers);
                R_temp(1,[uepair(m,:)]) = shape*log2(1+signel_power*ebno_range);
                R = [R; R_temp];
                sorted_uepair = [sorted_uepair; uepair(m,:)];
                
                C_real(m) = sum(log2(1+(signel_power*ebno_range)./(1+interference_power'*ebno_range)));

                clear signel_power;
                clear interference_power;
            end    
        end
        C_err = sum(R,2);
        C_real([find(C_real == 0)]) = [];
        C_real = C_real';
end

function [C_err, C_real, R, sorted_uepair, SINR, num_rank] = user_paired_4(numUsers, RBnumPerSB, syspar, HRB, ebno_range, phase_error, amp_error)
        num_user2paired = 4;
        spatial_freedom = [2 2 2 2];
        Nr = syspar.Nr;
        R = [];
        C_real = [];
        SINR = [];
        num_rank = [];
        sorted_uepair = [];
        %%%%%%%%%%%%%
        uepair=combntns(1:numUsers,num_user2paired);% different pair of UEs   % uepair=combntns(1:numUsers,numPair);
        for m=1:size(uepair,1) % a specific UE pairing
            u1=uepair(m,1);
            u2=uepair(m,2);
            u3=uepair(m,3);
            u4=uepair(m,4);
            %%%%%%%%%%%%
            H_u1 = HRB{u1};
            H_u2 = HRB{u2};
            H_u3 = HRB{u3};
            H_u4 = HRB{u4};
            for k = 1:RBnumPerSB
                H{k} = [H_u1(:,:,k); H_u2(:,:,k); H_u3(:,:,k);  H_u4(:,:,k)];
                H4user(:,:,m) = H{k};
                
                %   add channel error
                real_H = H{k};
                error_phase = pi*5/180*( phase_error(1) + (phase_error(2)-phase_error(1)).*rand(size(real_H,1), size(real_H,2)) );
                error_amp = 10.^( 0.1.*( amp_error(1) + (amp_error(2)-amp_error(1)).*rand(size(real_H,1), size(real_H,2)) ) );
                error_H = real_H.*error_amp.*exp(1i*error_phase);
                delta_H = real_H - error_H;
                
                if cond(error_H) > 10
                    continue;
                end
                
                [H_eff, M_BD] = BD_channel(error_H, num_user2paired, spatial_freedom);
                
                U = [];
                V = [];
                S = [];
                shape = [];
                start = 1;
                for u_ind = 1:num_user2paired
%                     [u s v] = svd( H_eff( (u_ind-1)*Nr+1:u_ind*Nr , (u_ind-1)*Nr+1:u_ind*Nr ) );
                    [u s v] = svd( H_eff( (u_ind-1)*Nr+1:u_ind*Nr , start:sum(spatial_freedom(1:u_ind)) ) );
                    start = start + spatial_freedom(u_ind);
                    r = diag(s);
                    ind = find( diag(s)/sum(diag(s)) < 0.1 );
                    num_rank(m,u_ind) = Nr - length(ind);
                    u(:,ind) = [];
                    v(:,ind) = [];
                    r(ind) = [];
                    U = blkdiag(U, u);
                    V = blkdiag(V, v);
                    S = [S; r];
                    shape = blkdiag(shape, ones(1,num_rank(m,u_ind)));
                end
                Q = U';
                
                signel_power = S.^2/length(S);
                
                delta_error_m = Q*delta_H*M_BD;
                for t = 1:size(delta_error_m,1)
                    interference_power(t) = (norm(delta_error_m(t,:))^2)/sum(num_rank(m,:));
                end
                
                R_temp = zeros(1,numUsers);
                R_temp(1,[uepair(m,:)]) = shape*log2(1+signel_power*ebno_range);
                R = [R; R_temp];
                sorted_uepair = [sorted_uepair; uepair(m,:)];
                
                C_real(m) = sum(log2(1+(signel_power*ebno_range)./(1+interference_power'*ebno_range)));

                clear signel_power;
                clear interference_power;
            end    
        end
        C_err = sum(R,2);
        C_real([find(C_real == 0)]) = [];
        C_real = C_real';
end