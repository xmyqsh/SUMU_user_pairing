function [cw EigenValue] = CodewordChoiceSB(ScenarioParam,FeedbackBit,Rank,H,RBnum)
% *************************************************************************
%  根据场景类型产
% 输入参数
%       ScenarioParam ---- Structural Element
%             AntennaNum        : Tx Antenna Num
%             AntennaConfig     : 0->ULA 1->DP
%             CorrFlag          : 0->High correlation 1->Non correlation
%             UEPolarizeSlant   : 45 or 90
%             
%       FeedbackBit             : Codewords Num =2^FeedbackBit
%       Rank                    : The rank of channel matrix

% 输出参数
%             Codebook          : 码本
%             INDEX             : 索引
%             cw               : Precoding Vector/Matrix

% --------------------------------------------------------------------------
% 修改日期        版本号     修改人	     更改原因
% --------------------------------------------------------------------------
% 2009-09-09      v1.0      陈艺戬 杨勋     创建
% 2009-09-11      v2.0      陈艺戬      增加码字选择算法
% **************************************************************************
AntennaNum = ScenarioParam.AntennaNum;
AntennaConfig = ScenarioParam.AntennaConfig;
CorrFlag = ScenarioParam.CorrFlag ;
UEPolarizeSlant = ScenarioParam.UEPolarizeSlant;
N = 2^FeedbackBit;
if 4 ==  AntennaNum                                           % 4Tx codebook
    switch AntennaConfig
       case 0                                                 % ULA Antenna Config
           switch Rank                                        % Different Rank,different codebook design
               case 1                                      
                  switch CorrFlag  
                     case 0                                   % High correlation 
                     Codebook = GetCodeBook_Antenna4_Rank1_ULA(FeedbackBit);  
%                      Codebook = uncorrCodebookRank1Tx4(FeedbackBit);
                                                              % Produce High correlation Rank1 ULA Codebook   
                                                              % saved as 2-dimension matrix:Codebook(:,index)
                     case 1                                   % Low correlation 
                     Codebook = uncorrCodebookRank1Tx4(FeedbackBit);   
                                                              % Produce High correlation Rank1 ULA Codebook     
                                                              % saved as 2-dimension matrix:Codebook(:,index)
                  end    
                     %     add  codeword choice Algorithm     %
%                      %%%%%%%%%% codeword choice start%%%%%%%%%%
%                      lamda = zeros(1,N);
%                      R = H'*H;
%                      for n = 1:N       
%                      lamda(n)=Codebook(:,n)'*R*Codebook(:,n) ;% Codeword choice rule
%                      end
%                      [TempValue Index] = max(lamda);
%                      EigenValue(1) = TempValue1/2;
% 
%                      cw = Codebook(:,Index);                 % return rank 1 precoding vector
%                      INDEX = Index;                           % return rank 1 precoding codeword index
%                      %%%%%%%%%% codeword choice end %%%%%%%%%%%

                     %%%%%%%%%% codeword choice start%%%%%%%%%%
                     

                     for RBIndex =1 : RBnum
                     [Utmp Dtmp Vtmp]= svd(H(:,:,RBIndex));
                     EVector(:,:,RBIndex) = Vtmp;
                     Evalue(:,RBIndex) =Dtmp(1,1)^2;
                     end

                     lamda = zeros(1,N);
                     for n = 1:N       
                       for RBIndex =1 : RBnum
                       lamda(n) = lamda(n) + sqrt(sum(sum(abs(Codebook(:,n)*Codebook(:,n)'-EVector(:,1,RBIndex)*EVector(:,1,RBIndex)').^2)))/sqrt(2);
                       end
                     end
                     
                     [TempValue Index] = min(lamda);
                     EigenValue = sum(Evalue)/RBnum;

                     cw = Codebook(:,Index);                 % return rank 1 precoding vector
                     INDEX = Index;                           % return rank 1 precoding codeword index
                     %%%%%%%%%% codeword choice end %%%%%%%%%%%

            
               case 2         
                   
                  
                   if mod(FeedbackBit,2)==1
                        error( 'Rank2码本开销暂限定为偶数');
                   elseif FeedbackBit<8
                        error( 'Rank2码本开销太小');
                   elseif FeedbackBit>18
                        error( 'Rank2码本开销太大');
                   end 
                   switch CorrFlag  
                      case 0        % High correlation 
                        error('4Tx ULA Rank2没有设计相关信道码字');
                      case 1        % Low correlation 
                     [Codebook CBTemp] = uncorrCodebookRank2Tx4(FeedbackBit);   
                                                              % saved as 3-dimension matrix:Codebook(:,:,index)
                   end
                   % Codebook

                            
                     %%%%%%%%%% codeword choice start%%%%%%%%%%
%                      [aa bb cc]=svd(H);
%                      lamda1 = zeros(1,sqrt(N));
%                      R = H'*H;
%                      for n = 1:sqrt(N)       
%                      lamda1(n) = CBTemp(:,n)'*R*CBTemp(:,n);  % codeword choice rule
%                      end
%                      [TempValue1 Index1] = max(lamda1);
%                      PMI1 = CBTemp(:,Index1); 
%                      EigenValue(1) = TempValue1/2;
%                      
%                      n = Index1;
%                      Vector_a = PMI1;
%                      for m = 1:sqrt(N)
%                          Vector_b = CBTemp(:,m);
%                          OrthonormalV_Temp(:,m) = (Vector_b-(Vector_a'* Vector_b)* Vector_a)/sqrt(sum(abs(Vector_b-(Vector_a'* Vector_b)* Vector_a).^2));
%                      end                     
%                      
%                      lamda2 = zeros(1,sqrt(N));
%                      for m=1:sqrt(N)       
%                      lamda2(m)=OrthonormalV_Temp(:,m)'* R * OrthonormalV_Temp(:,m) ;
%                                                              % codeword choice rule
%                      end
%                      lamda2(Index1)=0;
%                      [TempValue2 Index2] = max(lamda2);
%                      PMI2 = OrthonormalV_Temp(:,Index2);                     
%                      EigenValue(2) = TempValue2/2;
%                      cw = [PMI1 PMI2];
%                      INDEX = (Index1-1)*sqrt(N)+Index2;
              % ------------------------------------------------      
                     for RBIndex =1 : RBnum
                     [Utmp Dtmp Vtmp]= svd(H(:,:,RBIndex));
                     EVector(:,:,RBIndex) = Vtmp;
                     Evalue(1,RBIndex) =Dtmp(1,1)^2;
                     Evalue(2,RBIndex) =Dtmp(2,2)^2;
                     end

                     lamda1 =zeros(1,sqrt(N));
                     for n = 1:sqrt(N)
                       for RBIndex =1 : RBnum
                        lamda1(n) = lamda1(n)+ sqrt(sum(sum(abs(CBTemp(:,n)*CBTemp(:,n)'-EVector(:,1,RBIndex)*EVector(:,1,RBIndex)').^2)))/sqrt(2);
                       end
                     end
                     [TempValue1 Index1] = min(lamda1);
                     PMI1 = CBTemp(:,Index1); 
                     
                     n = Index1;
                     Vector_a = PMI1;
                     for m = 1:sqrt(N)
                         Vector_b = CBTemp(:,m);
                         OrthonormalV_Temp(:,m) = (Vector_b-(Vector_a'* Vector_b)* Vector_a)/sqrt(sum(abs(Vector_b-(Vector_a'* Vector_b)* Vector_a).^2));
                     end                      
                     lamda2 = zeros(1,sqrt(N));
                     for m=1:sqrt(N)       
                        for RBIndex =1 : RBnum
                        lamda2(m)=lamda2(m)+sqrt(sum(sum(abs(OrthonormalV_Temp(:,m)*OrthonormalV_Temp(:,m)'-EVector(:,2,RBIndex)*EVector(:,2,RBIndex)').^2)))/sqrt(2);
                        end
                     end
                     lamda2(Index1)=1000;
                     [TempValue2 Index2] = min(lamda2);
                     PMI2 = OrthonormalV_Temp(:,Index2);                     
                     cw = [PMI1 PMI2];
                     INDEX = (Index1-1)*sqrt(N)+Index2;
                     EigenValue(1,1) = sum(Evalue(1,:))/RBnum;
                     EigenValue(2,1) = sum(Evalue(2,:))/RBnum;
                     
                     %%%%%%%%%% codeword choice end %%%%%%%%%%%
               otherwise
                   error('暂不支持Rank>2');
            end      %  end rank switch      

       case 1                                                 % Dual Polarized(DP) Antenna Config
           switch Rank                                        % Different Rank,different codebook design
               case 1                                        
                  switch CorrFlag   
                     case 0                                   % High correlation 
                     [Codebook CBTemp] = highcorrCodebookRank1Tx4DP(FeedbackBit,UEPolarizeSlant);  
                                                              % Produce High correlation Rank1 DP Codebook
                                                              % Codewords design based on DP channel character                             
                     case 1                                   % Low correlation 
                     Codebook = uncorrCodebookRank1Tx4(FeedbackBit);   
                                                              % Produce Low correlation  Rank1 DP Codebook
                                                              % Grassmannian Codebook, the same as ULA antenna config      
                  end
                     %     add  codeword choice Algorithm     %
                     %%%%%%%%%% codeword choice start%%%%%%%%%%
                     lamda = zeros(1,N);
                     R = H'*H;
                     for n = 1:N       
                     lamda(n)=Codebook(:,n)'*R*Codebook(:,n) ;% Codeword choice rule
                     end
                     [TempValue Index] = max(lamda);
                     cw = Codebook(:,Index);                 % return rank 1 precoding vector
                     INDEX = Index;                           % return rank 1 precoding codeword index
                     %%%%%%%%%% codeword choice end %%%%%%%%%%%
             case 2   
                  switch CorrFlag  
                     case 0                                   % High correlation 
                     Codebook = GetCodeBook_Antenna4_Rank2_DP(FeedbackBit);   
                                                              % [w1 w1;-w1 w1];       
                     %%%%%%%%%% codeword choice start%%%%%%%%%%     
                     lamda = zeros(1,N);
                     R = H'*H;
                     for n = 1:N       
                     lamda(n)=Codebook(:,1,n)'*R*Codebook(:,1,n) ;% Codeword choice rule
                     end
                     [TempValue Index] = max(lamda);
                     cw = Codebook(:,:,Index);               % return rank 2 precoding matirx
                     INDEX = Index;                           % return rank 2 precoding codeword index
                     %%%%%%%%%% codeword choice end %%%%%%%%%%%
                     
                     case 1                                   % Low correlation 
                     [Codebook CBTemp]= uncorrCodebookRank2Tx4DP(FeedbackBit,UEPolarizeSlant);    
                                                              % [w1 w2;-w1 w2]  
                     %%%%%%%%%% codeword choice start%%%%%%%%%%
                     lamda1 = zeros(1,sqrt(N));
                     R = H'*H;
                     for n = 1:sqrt(N)  
                     Temp1(:,n) = [CBTemp(:,n);-CBTemp(:,n)];
                     lamda1(n)=Temp1(:,n)'*R*Temp1(:,n) ;        % Codeword choice rule
                     end
                     [TempValue1 Index1] = max(lamda1);
                     PMI1 = Temp1(:,Index1);
                    
                     n = Index1;
                     Vector_a = CBTemp(:,Index1);
                     for m = 1:sqrt(N)
                         Vector_b = CBTemp(:,m);
                         OrthonormalV_Temp(:,m) = (Vector_b-(Vector_a'* Vector_b)* Vector_a)/sqrt(sum(abs(Vector_b-(Vector_a'* Vector_b)* Vector_a).^2));
                         Temp2(:,m) = [OrthonormalV_Temp(:,m);OrthonormalV_Temp(:,m)];
                     end      
                     
                     lamda2 = zeros(1,sqrt(N));
                     for m=1:sqrt(N)       
                     lamda2(m)=Temp2(:,m)'* R * Temp2(:,m) ; % codeword choice rule
                     end   
                     lamda2(Index1)=0;
                     [TempValue2 Index2] = max(lamda2);
                     PMI2 = Temp2(:,Index2);
                     
                     cw = [PMI1 PMI2];
                     INDEX = (Index1-1)*sqrt(N)+Index2;
                     %%%%%%%%%% codeword choice end %%%%%%%%%%%                
                  end
              otherwise
                   error('暂不支持Rank>2');
           end
   end
end  
   
   
if(8== AntennaNum)                                        % 8Tx codebook
    switch AntennaConfig
       case 0                                                 % ULA Antenna Config
           switch Rank                                        % Different Rank,different codebook design
               case 1                                         
                  switch CorrFlag  
                     case 0                                   % High correlation 
                     Codebook = GetCodeBook_Antenna8_Rank1_ULA(FeedbackBit);   
                                                              % Produce High correlation Rank1 ULA Codebook   
                                                              % saved as 2-dimension matrix:Codebook(:,index)
                     case 1                                   % Low correlation 
                     Codebook = uncorrCodebookRank1Tx8(FeedbackBit);   
                                                              % Produce High correlation Rank1 ULA Codebook     
                                                              % saved as 2-dimension matrix:Codebook(:,index)
                  end    
                     %     add  codeword choice Algorithm     %
                     %%%%%%%%%% codeword choice start%%%%%%%%%%
                     lamda = zeros(1,N);
                     R = H'*H;
                     for n = 1:N       
                     lamda(n)=Codebook(:,n)'*R*Codebook(:,n) ;% Codeword choice rule
                     end
                     [TempValue Index] = max(lamda);
                     cw = Codebook(:,Index);                 % return rank 1 precoding vector
                     INDEX = Index;                           % return rank 1 precoding codeword index
                     %%%%%%%%%% codeword choice end %%%%%%%%%%%                 
                case 2                                       
                  
                   if mod(FeedbackBit,2)==1
                        error( 'Rank2码本开销暂限定为偶数');
                   elseif FeedbackBit<8
                        error( 'Rank2码本开销太小');
                   elseif FeedbackBit>18
                        error( 'Rank2码本开销太大');
                   end 
                 switch CorrFlag  
                    case 0        % High correlation 
                        error('8Tx ULA Rank2没有设计相关信道码字');
                    case 1        % Low correlation 
                     [Codebook CBTemp] = uncorrCodebookRank2Tx8(FeedbackBit);   
                                                              % saved as 3-dimension matrix:Codebook(:,:,index)
                 end
                            
                    %%%%%%%%%% codeword choice start%%%%%%%%%%
                     lamda1 = zeros(1,sqrt(N));
                     R = H'*H;
                     for n = 1:sqrt(N)       
                     lamda1(n) = CBTemp(:,n)'*R*CBTemp(:,n);  % codeword choice rule
                     end
                     [TempValue1 Index1] = max(lamda1);
                     PMI1 = CBTemp(:,Index1); 
                     
                     n = Index1;
                     Vector_a = PMI1;
                     for m = 1:sqrt(N)
                         Vector_b = CBTemp(:,m);
                         OrthonormalV_Temp(:,m) = (Vector_b-(Vector_a'* Vector_b)* Vector_a)/sqrt(sum(abs(Vector_b-(Vector_a'* Vector_b)* Vector_a).^2));
                     end                     
                     
                     lamda2 = zeros(1,sqrt(N));
                     for m=1:sqrt(N)       
                     lamda2(m)=OrthonormalV_Temp(:,m)'* R * OrthonormalV_Temp(:,m) ;
                                                             % codeword choice rule
                     end
                     lamda2(Index1)=0;
                     [TempValue2 Index2] = max(lamda2);
                     PMI2 = OrthonormalV_Temp(:,Index2);                     
                     
                     cw = [PMI1 PMI2];
                     INDEX = (Index1-1)*sqrt(N)+Index2;
                     
                     %%%%%%%%%% codeword choice end %%%%%%%%%%%
               otherwise
                   error('暂不支持Rank>2');                 
           end
           
           
       case 1                                                 % Dual Polarized(DP) Antenna Config
           switch Rank                                        % Different Rank,different codebook design
               case 1                                        
                  switch CorrFlag  
                     case 0                                   % High correlation 
                     Codebook = highcorrCodebookRank1Tx8DP(FeedbackBit,UEPolarizeSlant);  
                                                              % Produce High correlation Rank1 DP Codebook
                                                              % Codewords design based on DP channel character                             
                     case 1                                   % Low correlation 
                     Codebook = uncorrCodebookRank1Tx8(FeedbackBit);  
                                                              % Produce Low correlation  Rank1 DP Codebook
                                                              % Grassmannian Codebook, the same as ULA antenna config      
                  end
                     %     add  codeword choice Algorithm     %
                     %%%%%%%%%% codeword choice start%%%%%%%%%%     
                     lamda = zeros(1,N);
                     R = H'*H;
                     for n = 1:N       
                     lamda(n)=Codebook(:,n)'*R*Codebook(:,n) ;% Codeword choice rule
                     end
                     [TempValue Index] = max(lamda);
                     cw = Codebook(:,Index);                 % return rank 2 precoding matirx
                     INDEX = Index;                           % return rank 2 precoding codeword index
                     %%%%%%%%%% codeword choice end %%%%%%%%%%%
                  
              case 2   
                  switch CorrFlag  
                     case 0                                   % High correlation 
                     Codebook = GetCodeBook_Antenna8_Rank2_DP(FeedbackBit);  
                                                              % [w1 w1;-w1 w1];   
                     %%%%%%%%%% codeword choice start%%%%%%%%%%     
                     lamda = zeros(1,N);
                     R = H'*H;
                     for n = 1:N       
                     lamda(n)=Codebook(:,1,n)'*R*Codebook(:,1,n) ;% Codeword choice rule
                     end
                     [TempValue Index] = max(lamda);
                     cw = Codebook(:,:,Index);               % return rank 2 precoding matirx
                     INDEX = Index;                           % return rank 2 precoding codeword index
                     %%%%%%%%%% codeword choice end %%%%%%%%%%%                                                              
                                  
                                                              
                     case 1                                   % Low correlation 
                     [Codebook CBTemp] = uncorrCodebookRank2Tx8DP(FeedbackBit);    
                                                              % [w1 w2;-w1 w2]  
                     % Codebook                                         
                     %%%%%%%%%% codeword choice start%%%%%%%%%%
                     lamda1 = zeros(1,sqrt(N));
                     R = H'*H;
                     for n = 1:sqrt(N)  
                     Temp1(:,n) = [CBTemp(:,n);-CBTemp(:,n)];
                     lamda1(n)=Temp1(:,n)'*R*Temp1(:,n) ;        % Codeword choice rule
                     end
                     [TempValue1 Index1] = max(lamda1);
                     PMI1 = Temp1(:,Index1);
                    
                     n = Index1;
                     Vector_a = CBTemp(:,Index1);
                     for m = 1:sqrt(N)
                         Vector_b = CBTemp(:,m);
                         OrthonormalV_Temp(:,m) = (Vector_b-(Vector_a'* Vector_b)* Vector_a)/sqrt(sum(abs(Vector_b-(Vector_a'* Vector_b)* Vector_a).^2));
                         Temp2(:,m) = [OrthonormalV_Temp(:,m);OrthonormalV_Temp(:,m)];
                     end   
                     
                     lamda2 = zeros(1,sqrt(N));
                     for m=1:sqrt(N)       
                     lamda2(m)=Temp2(:,m)'* R * Temp2(:,m) ; % codeword choice rule
                     end              
                     lamda2(Index1)
                     [TempValue2 Index2] = max(lamda2);
                     PMI2 = Temp2(:,Index2);
                     
                     cw = [PMI1 PMI2]/sqrt(2);
                     INDEX = (Index1-1)*sqrt(N)+Index2;
                     %%%%%%%%%% codeword choice end %%%%%%%%%%% 
                   end
                 
                otherwise
                   error('暂不支持Rank>2');
           end
   end
end