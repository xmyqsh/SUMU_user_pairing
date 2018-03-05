function capacity=sim_MUCB(ctrlpar)
% codebook based multi-user MIMO platform
% multiuser precoding and scheduling
%   Example:
%       ctrlpar=ctrlparset;
%       sim_MUCB(ctrlpar);
if nargin>1, error('Too many input arguments.');end
tic
addpath('SCME');
% addpath('CodebookProduce and CodewordChoice');

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
%             Frame_time{ueIdx} = Frame_time{ueIdx} + 1000;  %1e-3;
%             [Hf(:,:,ueIdx) chan_par_backup{ueIdx}] =calc_channel(syspar,scmpar,linkpar,antpar,Frame_time{ueIdx},chan_par_backup{ueIdx});
%         end
%         
%     end
%     channel(:,:,:,drop) = Hf;
% 
% end
%     save('channel_50u_8t_2r.mat', 'channel')
%     
    load channel_50u_8t_2r.mat;
for drop = [1 2 10 100 1000]
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
        
        CMax=zeros(length(ebno_range),1);
        % MU-MIMO scheduling procedure
        % multiplexing of two users
        [C1, uepair1, SINR1, num_rank1] = user_paired_1(numUsers, RBnumPerSB, syspar, HRB, ebno_range);
        [C2, uepair2, SINR2, num_rank2] = user_paired_2(numUsers, RBnumPerSB, syspar, HRB, ebno_range);
        [C3, uepair3, SINR3, num_rank3] = user_paired_3(numUsers, RBnumPerSB, syspar, HRB, ebno_range);
        [C4, uepair4, SINR4, num_rank4] = user_paired_4(numUsers, RBnumPerSB, syspar, HRB, ebno_range);
%         [C5, uepair5, SINR5, num_rank5] = user_paired_5(numUsers, RBnumPerSB, syspar, HRB, ebno_range);
%         [C6, uepair6, SINR6, num_rank6] = user_paired_6(numUsers, RBnumPerSB, syspar, HRB, ebno_range);
%         [C7, uepair7, SINR7, num_rank7] = user_paired_7(numUsers, RBnumPerSB, syspar, HRB, ebno_range);
%         nlayer=zeros(2,1);
%         %%%%%%%%%%%%%
%         num_user2paired = 1;
%         %%%%%%%%%%%%%
%         uepair=combntns(1:numUsers,num_user2paired);% different pair of UEs   % uepair=combntns(1:numUsers,numPair);
%         for m=1:size(uepair,1) % a specific UE pairing
%             u1=uepair(m,1);
% %             u2=uepair(m,2);
% %             u3=uepair(m,3);
% %             u4=uepair(m,4);
%             %%%%%%%%%%%%
%             H_u1 = HRB{u1};
% %             H_u2 = HRB{u2};
% %             H_u3 = HRB{u3};
% %             H_u4 = HRB{u4};
%             for k = 1:RBnumPerSB
% %                 H{k} = [H_u1(:,:,k); H_u2(:,:,k); H_u3(:,:,k); H_u4(:,:,k)];
% %                 H{k} = [H_u1(:,:,k); H_u2(:,:,k)];
%                 H{k} = [H_u1(:,:,k)];
%                 % to be update
%                 % a function to cope with this
% %                 SINR(m,k) = sum(sum(H{k}));
%                 [L,Q,G,PPP]= GMD_encoder(H{k},num_user2paired,syspar.Nr*ones(num_user2paired,1));
%                 [r(m,:), num_rank(m,:)] = unique_r(diag(G));
%                 SINR(m,:) = (ebno_range.*r(m,:).^2)./syspar.Nr;
%             end    
% %                 % calculate EESMs
% %                 betas =ctrlpar.EESMbetas;
% %                 SINR_vector = SINR{m}';
% %                 na = length(SINR_vector);
% %                 nb = length(betas);
% %                 [ja,jb] = meshgrid(1:na,1:nb);
% %                 effective_SINR{m} = -(betas.').*log(exp(-SINR_vector(ja)./betas(jb)));
% %                 max_effec_SINR(m) = max(abs(effective_SINR{m}));
%             C(m) = sum(log2(1+SINR(m,:)).*num_rank(m,:)); 
%         end

%         [max_SINR, pair_select] = max(max_effec_SINR);
%         disp('pairing users: ')
%         uepair(pair_select,:)
%         [~,CQI] = min(abs(max_SINR-ctrlpar.TargetSINR));

%         [max_C, pair_select] = max([C1 C2./2 C3./3 C4./4]);
        
        [max_C1, pair_select1] = max(C1);max_C1,pair_select1
        [max_C2, pair_select2] = max(C2);max_C2,pair_select2
        [max_C3, pair_select3] = max(C3);max_C3,pair_select3
        [max_C4, pair_select4] = max(C4);max_C4,pair_select4
%         [max_C5, pair_select5] = max(C5);max_C5
%         [max_C6, pair_select6] = max(C6);max_C6
%         [max_C7, pair_select7] = max(C7);max_C7
%         disp('pairing users: ')
%         uepair(pair_select,:)
%         save('C1.mat','C')
%         save('result_20_users_v1.mat', 'max_C1', 'max_C2', 'max_C3', 'max_C4', 'max_C5', 'max_C6', 'max_C7', 'pair_select1', 'pair_select2', 'pair_select3', 'pair_select4', 'pair_select5', 'pair_select6', 'pair_select7')
    end
end

capacity=capacity_vector/(num_drops*size(Hf,1)/syspar.Nr);
% plot(10*log10(ebno_range), capacity,'b');
% grid on; xlabel('E_b/N_o (dB)'); ylabel('Throughput (Mbps)');
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

function [C, uepair, SINR, num_rank] = user_paired_1(numUsers, RBnumPerSB, syspar, HRB, ebno_range)
        num_user2paired = 1;
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
                % to be update
                % a function to cope with this
%                 SINR(m,k) = sum(sum(H{k}));
                [L,Q,G,PPP]= GMD_encoder(H{k},num_user2paired,syspar.Nr*ones(num_user2paired,1));
                [r(m,:), num_rank(m,:)] = unique_r(diag(G), num_user2paired);
                SINR(m,:) = (ebno_range.*r(m,:).^2)./(syspar.Nt);
            end
            C(m) = sum(log2(1+SINR(m,:)).*num_rank(m,:));
        end
end

function [C, uepair, SINR, num_rank] = user_paired_2(numUsers, RBnumPerSB, syspar, HRB, ebno_range)
        num_user2paired = 2;
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

                % to be update
                % a function to cope with this
%                 SINR(m,k) = sum(sum(H{k}));
                [L,Q,G,PPP]= GMD_encoder(H{k},num_user2paired,syspar.Nr*ones(num_user2paired,1));
                [r(m,:), num_rank(m,:)] = unique_r(diag(G), num_user2paired);
                SINR(m,:) = (ebno_range.*r(m,:).^2)./(syspar.Nt);
            end    
            C(m) = sum(log2(1+SINR(m,:)).*num_rank(m,:));
        end
end

function [C, uepair, SINR, num_rank] = user_paired_3(numUsers, RBnumPerSB, syspar, HRB, ebno_range)
        num_user2paired = 3;
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

                % to be update
                % a function to cope with this
%                 SINR(m,k) = sum(sum(H{k}));
                [L,Q,G,PPP]= GMD_encoder(H{k},num_user2paired,syspar.Nr*ones(num_user2paired,1));
                [r(m,:), num_rank(m,:)] = unique_r(diag(G), num_user2paired);
                SINR(m,:) = (ebno_range.*r(m,:).^2)./(syspar.Nt);
            end    
            C(m) = sum(log2(1+SINR(m,:)).*num_rank(m,:));
        end
end

function [C, uepair, SINR, num_rank] = user_paired_4(numUsers, RBnumPerSB, syspar, HRB, ebno_range)
        num_user2paired = 4;
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

                % to be update
                % a function to cope with this
%                 SINR(m,k) = sum(sum(H{k}));
                [L,Q,G,PPP]= GMD_encoder(H{k},num_user2paired,syspar.Nr*ones(num_user2paired,1));
                [r(m,:), num_rank(m,:)] = unique_r(diag(G), num_user2paired);
                SINR(m,:) = (ebno_range.*r(m,:).^2)./(syspar.Nt);
            end    
            C(m) = sum(log2(1+SINR(m,:)).*num_rank(m,:));
        end
end

function [C, uepair, SINR, num_rank] = user_paired_5(numUsers, RBnumPerSB, syspar, HRB, ebno_range)
        num_user2paired = 5;
        %%%%%%%%%%%%%
        uepair=combntns(1:numUsers,num_user2paired);% different pair of UEs   % uepair=combntns(1:numUsers,numPair);
%         r = zeros(size(uepair,1),num_user2paired);
%         num_rank = zeros(size(uepair,1),num_user2paired);
%         SINR = zeros(size(uepair,1),num_user2paired);
        for m=1:size(uepair,1) % a specific UE pairing
            u1=uepair(m,1);
            u2=uepair(m,2);
            u3=uepair(m,3);
            u4=uepair(m,4);
            u5=uepair(m,5);
            %%%%%%%%%%%%
            H_u1 = HRB{u1};
            H_u2 = HRB{u2};
            H_u3 = HRB{u3};
            H_u4 = HRB{u4};
            H_u5 = HRB{u5};
            for k = 1:RBnumPerSB
                H{k} = [H_u1(:,:,k); H_u2(:,:,k); H_u3(:,:,k);  H_u4(:,:,k);  H_u5(:,:,k)];

                % to be update
                % a function to cope with this
%                 SINR(m,k) = sum(sum(H{k}));
                [L,Q,G,PPP]= GMD_encoder(H{k},num_user2paired,syspar.Nr*ones(num_user2paired,1));
                [r(m,:), num_rank(m,:)] = unique_r(diag(G), num_user2paired);
                SINR(m,:) = (ebno_range.*r(m,:).^2)./(syspar.Nt);
            end    
            C(m) = sum(log2(1+SINR(m,:)).*num_rank(m,:));
        end
end

function [C, uepair, SINR, num_rank] = user_paired_6(numUsers, RBnumPerSB, syspar, HRB, ebno_range)
        num_user2paired = 6;
        %%%%%%%%%%%%%
        uepair=combntns(1:numUsers,num_user2paired);% different pair of UEs   % uepair=combntns(1:numUsers,numPair);
%         r = zeros(size(uepair,1),num_user2paired);
%         num_rank = zeros(size(uepair,1),num_user2paired);
%         SINR = zeros(size(uepair,1),num_user2paired);
        for m=1:size(uepair,1) % a specific UE pairing
            u1=uepair(m,1);
            u2=uepair(m,2);
            u3=uepair(m,3);
            u4=uepair(m,4);
            u5=uepair(m,5);
            u6=uepair(m,6);
            %%%%%%%%%%%%
            H_u1 = HRB{u1};
            H_u2 = HRB{u2};
            H_u3 = HRB{u3};
            H_u4 = HRB{u4};
            H_u5 = HRB{u5};
            H_u6 = HRB{u6};
            for k = 1:RBnumPerSB
                H{k} = [H_u1(:,:,k); H_u2(:,:,k); H_u3(:,:,k);  H_u4(:,:,k);  H_u5(:,:,k);  H_u6(:,:,k)];

                % to be update
                % a function to cope with this
%                 SINR(m,k) = sum(sum(H{k}));
                [L,Q,G,PPP]= GMD_encoder(H{k},num_user2paired,syspar.Nr*ones(num_user2paired,1));
                [r(m,:), num_rank(m,:)] = unique_r(diag(G), num_user2paired);
                SINR(m,:) = (ebno_range.*r(m,:).^2)./(syspar.Nt);
            end    
            C(m) = sum(log2(1+SINR(m,:)).*num_rank(m,:));
        end
end

function [C, uepair, SINR, num_rank] = user_paired_7(numUsers, RBnumPerSB, syspar, HRB, ebno_range)
        num_user2paired = 7;
        %%%%%%%%%%%%%
        uepair=combntns(1:numUsers,num_user2paired);% different pair of UEs   % uepair=combntns(1:numUsers,numPair);
%         r = zeros(size(uepair,1),num_user2paired);
%         num_rank = zeros(size(uepair,1),num_user2paired);
%         SINR = zeros(size(uepair,1),num_user2paired);
        for m=1:size(uepair,1) % a specific UE pairing
            u1=uepair(m,1);
            u2=uepair(m,2);
            u3=uepair(m,3);
            u4=uepair(m,4);
            u5=uepair(m,5);
            u6=uepair(m,6);
            u7=uepair(m,7);
            %%%%%%%%%%%%
            H_u1 = HRB{u1};
            H_u2 = HRB{u2};
            H_u3 = HRB{u3};
            H_u4 = HRB{u4};
            H_u5 = HRB{u5};
            H_u6 = HRB{u6};
            H_u7 = HRB{u7};
            for k = 1:RBnumPerSB
                H{k} = [H_u1(:,:,k); H_u2(:,:,k); H_u3(:,:,k);  H_u4(:,:,k);  H_u5(:,:,k);  H_u6(:,:,k);  H_u7(:,:,k)];

                % to be update
                % a function to cope with this
%                 SINR(m,k) = sum(sum(H{k}));
                [L,Q,G,PPP]= GMD_encoder(H{k},num_user2paired,syspar.Nr*ones(num_user2paired,1));
                [r(m,:), num_rank(m,:)] = unique_r(diag(G), num_user2paired);
                SINR(m,:) = (ebno_range.*r(m,:).^2)./(syspar.Nt);
            end    
            C(m) = sum(log2(1+SINR(m,:)).*num_rank(m,:));
        end
end