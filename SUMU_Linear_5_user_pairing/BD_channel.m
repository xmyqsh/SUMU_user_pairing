%% BD algorithm precoder
%in this algorithm, Rx stands for the numbers of
%receiver antennas,similary to Tx stands for transmitter number
%%vision:
%%author:
%%goal :
%in this, Modmatrix stands for precode matrix of BD algorithm
clear 
clc
Tx=8;%the number of transimter
Rx=8;%the number of receiver
H_1=sqrt(0.5)*(randn(Rx,Tx)+sqrt(-1)*randn(Rx,Tx));%random channel

H=H_1;

%% 只需将H_estimate改为H_1即可
Users=4;%the number of users
sigma=0.6*randn(Rx,1)+j*0.8*randn(Rx,1);
sigma_2=sigma.*conj(sigma);%noise power
N_users=Rx/Users;% we assume all users have  same antennas 
%although in some situation this state would not practice
P=1:1:25;% power constranit in db
C_bd=zeros(length(P),1);

%************************************************************************%
%this section would get the Modmatrix_1 where the Modmatrix_1 stands for
%the left part of Modmatrix 
%*************************************************************************%

for i=1:Users
   H_k(:,:,i)= H(N_users*(i-1)+1:N_users*i,:);%H_k stands for block matrix

end
Hj(:,:,1)=H((Rx/Users)+1:Rx,:);
for q=2:Users-1
    Hj(:,:,q)=[H(1:N_users*(q-1),:);H((N_users*q+1):Rx,:)];%Hj means remain matrix that H eliminate H_k
end
Hj(:,:,Users)=H(1:N_users*(Users-1),:);
for k=1:Users
    [u(:,:,k),s(:,:,k),v(:,:,k)]=svd(Hj(:,:,k));%奇异值分解
    b=v(:,:,k);%将矩阵v赋予一个变量以便可以使用矩阵里面的一部分
    c=rank(Hj(:,:,k))+1;%求矩阵Hj（:,:,k）的秩，也是求v(:,:,k)非零列向量的个数。
     v_0(:,:,k)=b(:,c:Tx);%求零空间基
    Hs(:,:,k)=H_k(:,:,k)*v_0(:,:,k);
    [u_1(:,:,k),s_1(:,:,k),v_1(:,:,k)]=svd(Hs(:,:,k));%对Hs(:,:,k)做奇异值分解，以便得到矩阵
    %s_1(:,:,k)和v_1(:,:,k)，s_1(:,:,k)代表着奇异值矩阵，在以后用于注水功率分配，v_1(:,:,k)
    %代表右特征矩阵，其中的一部分用于做零空间的线性组合，具体参照文献的式(11)
    v_q=v_1(:,:,k);
    rhs=rank(Hs(:,:,k));
    v_1_1(:,:,k)=v_q(:,1:rhs);
    Ms(:,:,k)=v_0(:,:,k)*v_1_1(:,:,k);%对零空间的基做线性组合，以此可以求出预编码阵的一部分
end
temp=[];
for va=1:Users
  var_M_trix= Ms(:,:,va);
  temp=[temp var_M_trix];
end
  Modmatrix_1=temp;%这里就是将分块矩阵按行排列起来
%-*********************************************************************end
%% in this we get the left matrix 
%% using water-filling to power allocation
%**************************************************************************
%*this section would get the Modmatrix_2 where the Modmatrix_2 stands for
%the right of the part of Modmatrix  each element of the  Modmatrix_2 stands 
%for coefficients by water-filling methods
%**************************************************************************
eigen=[];
for var=1:Users
    matrix_cutted=block_diag_cut(s_1(:,:,var));
    var_vector=diag(matrix_cutted);
   eigen=[eigen;var_vector];
end %得到矩阵的非零奇异值排列成对角阵组成的列向量
eigenvalu_2=eigen.*conj(eigen);%求矩阵的非零奇异值排列成对角阵组成的列向量的平方
eigenvalue=diag(eigen);%得到矩阵的非零奇异值排列成对角阵的平方
diagmatrix_rank=rank(eigenvalue);%求非零奇异值排列成对角阵的秩
% water filling 
for temp=1:length(P)
[Total_rate,Power_per] = waterfill(P(temp),eigenvalu_2);%调用注水算法，函数参数为功率和所转化信道的奇异值
Modmatrix_2=sqrt(diag(Power_per));
%********************************************************************
%****************************************************************end
% get the Modmatrix
Modmatrix=Modmatrix_1*Modmatrix_2;
y=H_1*Modmatrix;% to check
C_bd(temp)=log2(det(eye(rank(H_1))+y*y'));
end
figure 
plot(P,C_bd,'*')








    
    
    

























