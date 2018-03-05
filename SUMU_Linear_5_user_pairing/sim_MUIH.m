function capacity=sim_MUIH(ctrlpar)
%   Example:
%       ctrlpar=ctrlparset;
%       sim_MUIH(ctrlpar);
if nargin>1, error('Too many input arguments.');end
tic
addpath('SCME');

[syspar,scmpar,antpar]=argdeliv(ctrlpar);

ebno_range=syspar.ebno_range;
num_drops=syspar.ndrops;
numUsers=syspar.nusers;
capacity_vector=zeros(length(ebno_range),1);
FBgran = syspar.FBgran;
maxlayer = syspar.maxlayer;
for drop = 1:num_drops
    fprintf('drop: %d/%d\n', drop, num_drops);
    Hf=zeros(syspar.ntimesamples*syspar.Nr*syspar.ntones,syspar.Nt,numUsers);
    for ueIdx=1:numUsers
        linkpar=linkparset;
        if ueIdx==1
            d=linkpar.MsBsDistance;
        else
            % users are randomly dropped with the same geometry
            linkpar.MsBsDistance=d;
        end
        Hf(:,:,ueIdx) =calc_channel(syspar,scmpar,linkpar,antpar);
    end

    if FBgran == 'RB Level'
    for k=1:syspar.Nr:size(Hf,1)
        % fprintf('Drop:%.1d, Processing: %.2f%% \n', drop,100*k/size(Hf,1));
        m=k+(0:syspar.Nr-1);
        CMax=zeros(length(ebno_range),1);
        nlayer=zeros(2,1);% 仅支持2个UE进行MU-MIMO
        uepair=combntns(1:numUsers,2);
        for r=1:size(uepair,1)
            u1=uepair(r,1);
            u2=uepair(r,2);
            H1=Hf(m,:,u1);
            nlayer(1)=min(rank(H1/norm(H1),0.14),syspar.maxlayer);% 秩自适应
            H2=Hf(m,:,u2);
            nlayer(2)=min(rank(H2/norm(H2),0.14),syspar.maxlayer);
            for ebnodex=1:length(ebno_range)
                ebno=ebno_range(ebnodex);
                [w1 w2]=txweight(H1,H2,ebno,syspar.precoding,nlayer);
                l=size([w1 w2],2);% 信号发射功率归一化
                w1=w1./sqrt(l);
                w2=w2./sqrt(l);
                % trace([w1 w2]'*[w1 w2])
                hs1=H1*w1;
                hs2=H2*w2;
                hi1=H1*w2;
                hi2=H2*w1;
                sinr1=msinr(hs1,hi1,ebno);
                % 10*log10(sinr1)
                % pause
                sinr2=msinr(hs2,hi2,ebno);
                CSum =SinrC(sinr1)+ SinrC(sinr2);
                CMax(ebnodex)=max(CSum,CMax(ebnodex));
                
                if ctrlpar.AutoSwitch == 1
                [ aa bb cc] =svd(H1);
                hs = H1*cc(:,1:nlayer(1));
                snr = msnr(hs,ebno);
                CSU = SinrC(snr) ;  
                CMax(ebnodex)= max(CSU,CMax(ebnodex));  
                end
                
                
            end
        end
        capacity_vector=capacity_vector+CMax;
    end


    elseif FBgran == 'SB Level'
    RBFreqSB = 5 ;
    TTIgran  = 1 ;
    FBdelay  = 6 ; % 以6ms为单位,1表示1个6ms
    RBnumPerSB = RBFreqSB*TTIgran;
    
    for k=1:(RBnumPerSB*syspar.Nr):size(Hf,1) % 重复执行的次数与RB的数目相同
        % disp(sprintf('drop:%.1d, Processing: %.2f%% \n',...
        %    drop,100*k/size(Hf,1)));
        m=k+(0:RBnumPerSB*syspar.Nr-1);
        Hq=cell(numUsers,1);
        for u=1:numUsers
            Htmp(1:RBnumPerSB*syspar.Nr,:) = Hf(m,:,u);
            HRBtmp1 = reshape(Htmp.',[],RBnumPerSB);
            HRBtmp2 = reshape(HRBtmp1,syspar.Nt,syspar.Nr,[]);
            HRB = permute(HRBtmp2,[2 1 3]);
           
            RRB(:,:,u) = zeros(syspar.Nt);
            for RBIndex = 1 : RBnumPerSB
            RRB(:,:,u) = RRB(:,:,u) + HRB(:,:,RBIndex)'*HRB(:,:,RBIndex);
            end
            RRB(:,:,u) =RRB(:,:,u)/RBnumPerSB;
            nstream=min(rank(RRB(:,:,u)/norm(RRB(:,:,u)),0.14),maxlayer);% 秩自适应
            [uu dd vv]  =  svd(RRB(:,:,u));
            d = dd(1:nstream,1:nstream);
            v = vv(:,1:nstream);

            f=struct('v',v,'d',d);
            hq=struct('f',f,'nlayer',nstream);
            Hq{u} = hq;
        end
        
        CMax=zeros(length(ebno_range),1);
        nlayer=zeros(2,1);
        uepair=combntns(1:numUsers,2);% different pair of UEs
        for r=1:size(uepair,1) % a specific UE pairing
            u1=uepair(r,1);
            u2=uepair(r,2);
            % [u1 u2]
            if max(m)+FBdelay*syspar.Nr*syspar.ntones<=size(Hf,1)
            H1=Hf(m+FBdelay*syspar.Nr*syspar.ntones,:,u1);
            f1=Hq{u1}.f;
            nlayer(1)=Hq{u1}.nlayer;
            H2=Hf(m+FBdelay*syspar.Nr*syspar.ntones,:,u2);
            f2=Hq{u2}.f;
            nlayer(2)=Hq{u2}.nlayer;
            else
            H1=Hf(m,:,u1);
            f1=Hq{u1}.f;
            nlayer(1)=Hq{u1}.nlayer;
            H2=Hf(m,:,u2);
            f2=Hq{u2}.f;
            nlayer(2)=Hq{u2}.nlayer;
                
            end    
            for ebnodex=1:length(ebno_range)
                ebno=ebno_range(ebnodex);
                [w1 w2]=txweight2(f1,f2,ebno,syspar.precoding,nlayer);
                l=size([w1 w2],2);% 信号发射功率归一化
                w1=w1./sqrt(l);
                w2=w2./sqrt(l);
                for RBIndex = 1 : RBnumPerSB
                hs1=H1(syspar.Nr*(RBIndex-1)+1:syspar.Nr*RBIndex,:)*w1;
                hs2=H2(syspar.Nr*(RBIndex-1)+1:syspar.Nr*RBIndex,:)*w2;
                hi1=H1(syspar.Nr*(RBIndex-1)+1:syspar.Nr*RBIndex,:)*w2;% inter-user interference
                hi2=H2(syspar.Nr*(RBIndex-1)+1:syspar.Nr*RBIndex,:)*w1;
                sinr1=msinr(hs1,hi1,ebno);
                % 10*log10(sinr1)
                sinr2=msinr(hs2,hi2,ebno);
                CSum(RBIndex) =SinrC(sinr1)+ SinrC(sinr2);
                end
                C_MU = sum(CSum);
                CMax(ebnodex)=max(C_MU,CMax(ebnodex));
                
                if ctrlpar.AutoSwitch == 1

                for RBIndex = 1 : RBnumPerSB
                hs = H1(syspar.Nr*(RBIndex-1)+1:syspar.Nr*RBIndex,:)*f1.v;
                snr = msnr(hs,ebno);
                CSU(RBIndex) = SinrC(snr) ;
                end
                C_SU =sum(CSU);
                CMax(ebnodex)= max(C_SU,CMax(ebnodex));  
                
                f_1.v = f1.v(:,1);
                f_1.d = f1.d(1,1);
                f_2.v = f2.v(:,1);
                f_2.d = f2.d(1,1);   
                [w1 w2]=txweight2(f_1,f_2,ebno,syspar.precoding,[1;1]);
                l=size([w1 w2],2);% 信号发射功率归一化
                w1=w1./sqrt(l);
                w2=w2./sqrt(l);
                for RBIndex = 1 : RBnumPerSB
                hs1=H1(syspar.Nr*(RBIndex-1)+1:syspar.Nr*RBIndex,:)*w1;
                hs2=H2(syspar.Nr*(RBIndex-1)+1:syspar.Nr*RBIndex,:)*w2;
                hi1=H1(syspar.Nr*(RBIndex-1)+1:syspar.Nr*RBIndex,:)*w2;% inter-user interference
                hi2=H2(syspar.Nr*(RBIndex-1)+1:syspar.Nr*RBIndex,:)*w1;
                sinr1=msinr(hs1,hi1,ebno);
                % 10*log10(sinr1)
                sinr2=msinr(hs2,hi2,ebno);
                CRank1MU(RBIndex) =SinrC(sinr1)+ SinrC(sinr2);
                end
                C_Rank1MU = sum(CRank1MU);
                CMax(ebnodex)= max(C_Rank1MU,CMax(ebnodex));  
                
                end
            end
          end
        capacity_vector=capacity_vector+CMax;
    end   
        
    
    end
    
end

capacity=capacity_vector/(num_drops*size(Hf,1)/syspar.Nr);
% plot(10*log10(ebno_range), capacity,'b');
% grid on; xlabel('E_b/N_o (dB)'); ylabel('Throughput (Mbps)');
toc
rmpath('SCME');

function [w1 w2]=txweight(H1,H2,snr,precoding,nlayer)
if size(H1,2)~=size(H2,2);error('Nt必须相同');end
if length(nlayer)~=2;error('必须是2个UE进行配对');end
if snr<1e-4;error('SNR太小');end
switch precoding
    case 'PEV'
        [u s v1]=svd(H1);
        w1=v1(:,1:nlayer(1));
        [u s v2]=svd(H2);
        w2=v2(:,1:nlayer(2));
    case 'SLNR'
        R1=H1'*H1;
        R2=H2'*H2;
        [V D]=eig(inv(R2+1/snr*eye(size(R2,1)))*R1);
        a=diag(real(D));
        d=sort(a,'descend');
        m=logical(a>=d(nlayer(1)));
        w1=V(:,m);
        
        [V D]=eig(inv(R1+1/snr*eye(size(R1,1)))*R2);
        a=diag(real(D));
        d=sort(a,'descend');
        m=logical(a>=d(nlayer(2)));
        w2=V(:,m);
    case 'BD'
        % Block Diagonalization, R1-093996
        G1=eye(size(H2,2))-pinv(H2)*H2; % 参考R1-090596
        [u s v]=svd(H1*G1);
        T1=v(:,1:nlayer(1));
        w1=G1*T1;
        % w1'*w1
        % u'*H1*w1

        G2=eye(size(H1,2))-pinv(H1)*H1;
        [u s v]=svd(H2*G2);
        T2=v(:,1:nlayer(2));
        w2=G2*T2;
        % w2'*w2
        % u'*H2*w2
    case 'RZF'
        H=[H1;H2];
        [Nr Nt]=size(H);
        if Nr>Nt
           error('Nr必须不大于Nt');
        end
        W=H'*inv(H*H'+1/snr*eye(size(H,1)));
        W=W./repmat(sqrt(diag(W'*W))',size(W,1),1);% 对列矢量进行归一化
        W1=W(:,1:size(H1,1));
        W2=W(:,(1:size(H2,1))+size(H1,1));
        w1=W1(:,1:nlayer(1));
        w2=W2(:,1:nlayer(2));
    otherwise
        w1=0;
        w2=0;
        error('');
end


function [w1 w2]=txweight2(f1,f2,snr,precoding,nlayer)
% precoding matrices w1 and w2 are computed
if length(nlayer)~=2;error('必须是2个UE进行配对');end
if size(f1.v,2)~=nlayer(1)
    error('');
end
if size(f2.v,2)~=nlayer(2)
    error('');
end
switch precoding
    case 'PEV'
        w1=f1.v;
        w2=f2.v;
    case 'SLNR'
        R1=f1.v*f1.d*f1.v';
        R2=f2.v*f2.d*f2.v';
        [V D]=eig(inv(R2+1/snr*eye(size(R2,1)))*R1);
        % real([trace(R2) trace(1/snr*eye(size(R2,1)))])
        a=diag(real(D));
        d=sort(a,'descend');
        m=logical(a>=d(nlayer(1)));
        w1=V(:,m);
        
        [V D]=eig(inv(R1+1/snr*eye(size(R1,1)))*R2);
        a=diag(real(D));
        d=sort(a,'descend');
        m=logical(a>=d(nlayer(2)));
        w2=V(:,m);
        % [chord(f1.v,f2.v) chord(w1,w2)]
    case 'BD'
        H1=f1.v';% Block Diagonalization, R1-093996
        H2=f2.v';
        G1=eye(size(H2,2))-pinv(H2)*H2; % 参考R1-090596
        [u s v]=svd(H1*G1);
        T1=v(:,1:nlayer(1));
        w1=G1*T1;
        % w1'*w1
        % u'*H1*w1

        G2=eye(size(H1,2))-pinv(H1)*H1;
        [u s v]=svd(H2*G2);
        T2=v(:,1:nlayer(2));
        w2=G2*T2;
        % w2'*w2
        % u'*H2*w2
        % [chord(f1.v,f2.v) chord(w1,w2)]
    case 'RZF'
        H1=f1.v';
        H2=f2.v';
        H=[H1;H2];
        W=H'*inv(H*H'+1/snr*eye(size(H,1)));
        W=W./repmat(sqrt(diag(W'*W))',size(W,1),1);% 对列矢量进行归一化
        w1=W(:,1:nlayer(1));
        w2=W(:,nlayer(1)+(1:nlayer(2)));
        % [chord(f1.v,f2.v) chord(w1,w2)]
    otherwise
        w1=0;
        w2=0;
        error('');
end
