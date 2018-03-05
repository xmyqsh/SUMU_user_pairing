function  [Q,R,P] =BD_GMD_5_stream(H,UserNum,UserAntennaNum,flag)

[m,n]=size(H);
I=eye(n);

if UserNum~=2
    
    [Q1,R1,P1] = BD_GMD_5_stream(H,2,[UserAntennaNum(1) sum(UserAntennaNum(2:end))],0);%¼Ù¶þÓÃ»§
    H2=H(UserAntennaNum(1)+1:end,:);
    
    [Q2,R2,P2] = BD_GMD_5_stream(H2*(I-P1*P1'),UserNum-1,UserAntennaNum(2:end),1);
    
    Q=blkdiag(Q1,Q2);
    R=blkdiag(R1,R2);
    R(length(diag(R1))+1:end,1:length(diag(R1)))=Q2'*H2*P1;
    P=[P1 P2];
    
else
    if flag==0
        H1=H(1:UserAntennaNum(1),:);
        
        [Q,R,P]=gmd_zcy_streamreduce(H1);
        
    else
        
        H1=H(1:UserAntennaNum(1),:);
        H2=H(UserAntennaNum(1)+1:UserAntennaNum(1)+UserAntennaNum(2),:);
        
        [Q1,R1,P1]=gmd_zcy_streamreduce(H1);
        
        HH=H2*(I-P1*P1');
        [Q2,R2,P2]=gmd_zcy_streamreduce(HH);
        r=(Q2)'*H2*P1;
        R=blkdiag(R1,R2);
        [row col]=size(R1);
        R(length(diag(R1))+1:length(diag(R1))+length(diag(R2)),1:col)=r;
        
        P=[P1  P2];
        
        Q=blkdiag(Q1,Q2);
    end
    
end