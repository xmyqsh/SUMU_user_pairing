function  [L,Q,GG,P] = GMD_encoder(H,UserNum,UserAntennaNum) 
% H=QRP'

flag=1;
if UserNum > 1
    [Q,R,P] = BD_GMD_5_stream(H,UserNum,UserAntennaNum,flag);
else       %    UserNum == 1
    H1=H(1:UserAntennaNum(1),:);
    [Q,R,P]=gmd_zcy_streamreduce(H1);
end
GG=diag(diag(R));
G=diag(1./diag(R));
L=R;
B=G*R;

end


