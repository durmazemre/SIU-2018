% Step size sabit
%Nakagami deneme
%1-prob + semilog

clear; clc; close all;dbstop if error;
%% Defining the System
p = 1; %Signal power 
MaxIt = 2000;
SNR = [0:1:20];  
epsilon = 1/100;  
f1 = figure;
f2 = figure;
min_node = 3;
max_node = 8;

%% This is for naming the data
number = [min_node:1:max_node]';
names = int2str(number);
                                
Allowed_Error = 1/10000; 
for nodes = min_node:1:max_node
    
N = nodes^2 ;
Edge_num = nchoosek(nodes,2);
monte_counter = 0;
counter = 0;
switch nodes
      case 3
          Measured = [-6 ; -3 ; 7];   
      case 4
          Measured = [-6 ; -3 ; 7 ; 14 ];
      case 5
          Measured = [-6 ; -3 ; 7 ; 14; 21];
      case 6
          Measured = [-6; -3; 7; 14; 21; -12.5];
      case 7
          Measured = [-6; -3; 7; 14; 17; -12.5; 8.5];
      case 8
          Measured = [-6; -3; 7; 14; 17; -12.5; 8.5; -9.5];
    end     
IterMeasured = Measured;
%% Calculating Average
sz = size(Measured);
Mysum = sum(Measured);
Average = Mysum/sz(1);
%% Starting to Monte Carlo
montemax = 1e5;                                                
L_Networks = zeros(montemax,nodes,nodes); 
Ranks = zeros(montemax,6);              
for SNR_counter = 1:length(SNR) 
    tic 
for monte = 1:montemax
    %Nakagami başlangıcı
    R = 1 ; 
    Naka_m=2;
    M=100000;
    h=[];
for d=1: 2*Naka_m
 h(d,:)= (((randn(1,Edge_num))/sqrt(2*Naka_m)));

end
      h=h.^2;
      
if Naka_m==0.5    %%Gauss
        h_Naka=sqrt(h);
else
        h_Naka=sqrt( sum(h));
end

gama = 10^(SNR(SNR_counter)/10); 
threshold = (2^R - 1) / gama ;
 Edge_con = abs(h_Naka).^2 > threshold ;

StateContainer = zeros(size(Measured,1), MaxIt + 1 ); 
     
     StateContainer(:,1) = Measured(:); 
     
     FeedbackContainer = zeros(MaxIt,size(Measured,2));  


      
     
     [ii,jj] = ndgrid(1:nodes);            
      A = zeros(nodes);                       
      A(jj>ii) =  Edge_con;                   
     
      
      A = A + A';                             
      D = diag(sum(A));                        
      L = D - A ;                             
           
      L_Networks(monte,:,:) = L ; 

     
      Ranks(monte,1) = rank(L) ;  
      
      %%  monte carlo
      
      if rank(L) == nodes - 1 
     
      Ranks(monte,4) = 1; 
      monte_counter = monte_counter + 1;
      
      for k = 1:MaxIt

      P_epsilon = eye(nodes) - epsilon * L;
      IterMeasured = P_epsilon * IterMeasured;  
      StateContainer(1:end,k) = IterMeasured(1:end);
      %% Defining an Error Rule
      if abs(max(IterMeasured) - min(IterMeasured)) > Allowed_Error 
      counter = counter + 1;
      else
          Ranks(monte,5) = counter;
          IterMeasured = Measured;
          counter = 0;
          break 
      end
            end
            else
      Ranks(monte,4) = 0;
      Ranks(monte,5) = 0; 
      end
      
end
%% Calculating the Probabilities of Ranks

   Average_It_per_SNR = sum(Ranks(:,5))/sum(Ranks(:,4));
   Average_It_Container(1,SNR_counter) = Average_It_per_SNR;
   Succ_Probability = sum(Ranks(:,4)) / montemax;
   Probability_Container(1,SNR_counter) =  Succ_Probability;
  toc 
end
figure(f1);
hold on;
semilogy(SNR,1-Probability_Container,'x -'); 
legend(names);
title('Propability of Reaching Consensus vs. Different SNR');

figure(f2);
hold on;
semilogy(SNR,Average_It_Container,'x -');
legend(names);
title('Average Number of Iterations vs. Different SNR');

end
