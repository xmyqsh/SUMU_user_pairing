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

%% ֻ�轫H_estimate��ΪH_1����
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
    [u(:,:,k),s(:,:,k),v(:,:,k)]=svd(Hj(:,:,k));%����ֵ�ֽ�
    b=v(:,:,k);%������v����һ�������Ա����ʹ�þ��������һ����
    c=rank(Hj(:,:,k))+1;%�����Hj��:,:,k�����ȣ�Ҳ����v(:,:,k)�����������ĸ�����
     v_0(:,:,k)=b(:,c:Tx);%����ռ��
    Hs(:,:,k)=H_k(:,:,k)*v_0(:,:,k);
    [u_1(:,:,k),s_1(:,:,k),v_1(:,:,k)]=svd(Hs(:,:,k));%��Hs(:,:,k)������ֵ�ֽ⣬�Ա�õ�����
    %s_1(:,:,k)��v_1(:,:,k)��s_1(:,:,k)����������ֵ�������Ժ�����עˮ���ʷ��䣬v_1(:,:,k)
    %�����������������е�һ������������ռ��������ϣ�����������׵�ʽ(11)
    v_q=v_1(:,:,k);
    rhs=rank(Hs(:,:,k));
    v_1_1(:,:,k)=v_q(:,1:rhs);
    Ms(:,:,k)=v_0(:,:,k)*v_1_1(:,:,k);%����ռ�Ļ���������ϣ��Դ˿������Ԥ�������һ����
end
temp=[];
for va=1:Users
  var_M_trix= Ms(:,:,va);
  temp=[temp var_M_trix];
end
  Modmatrix_1=temp;%������ǽ��ֿ��������������
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
end %�õ�����ķ�������ֵ���гɶԽ�����ɵ�������
eigenvalu_2=eigen.*conj(eigen);%�����ķ�������ֵ���гɶԽ�����ɵ���������ƽ��
eigenvalue=diag(eigen);%�õ�����ķ�������ֵ���гɶԽ����ƽ��
diagmatrix_rank=rank(eigenvalue);%���������ֵ���гɶԽ������
% water filling 
for temp=1:length(P)
[Total_rate,Power_per] = waterfill(P(temp),eigenvalu_2);%����עˮ�㷨����������Ϊ���ʺ���ת���ŵ�������ֵ
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








    
    
    

























