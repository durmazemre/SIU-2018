% Step size sabit
%Rician deneme
%1-prob + semilog
%Hata verirse row-column deðiþtir rician ile ilgili olan yerde


clear; clc; close all;dbstop if error;
%% Defining the System
p = 1; %Signal power 
MaxIt = 2000;
SNR = [1:1:20];  
epsilon = 1/100;  %step size determined random
f1 = figure;
f2 = figure;
min_node = 3;
max_node = 8;

%% This is for naming the data
number = [min_node:1:max_node]';
names = int2str(number);
                                
Allowed_Error = 1/1000; 
for nodes = min_node:1:max_node
    
N = nodes^2 ;
Edge_num = nchoosek(nodes,2);%for us to be able to define any number of edges
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
montemax = 1e5;                           %how many times monte carlo                         
L_Networks = zeros(montemax,nodes,nodes); %thil is L matrix container
Ranks = zeros(montemax,6);                %this will keep everything about matrices
for SNR_counter = 1:length(SNR) 
    tic %alttaki yazý 20 defa tekrarlanacak.Word'e kopyala ve onlarý say. Kalan zamaný hesapla.
for monte = 1:montemax
 % Rician baþlangýcý    
    R = 1 ;                            % Data rate = 1 alindi
    K=3;					
    mu = sqrt( K/(2*(K+1)) );
    s = sqrt( 1/(2*(K+1)) );

    %h_Rician=abs( s*randn(N,1) + mu ) + 1i*( s*randn(N,1) + mu );         %Rician fading    % N deðiþkeni Edge_num olabilir mi?   
    
    StateContainer = zeros(size(Measured,1), MaxIt + 1 ); %for every monte, this will keep states
    StateContainer(:,1) = Measured(:); %içine alamýyor tranpoz sorunu  
    FeedbackContainer = zeros(MaxIt,size(Measured,2));   %for every monte this will keep feedbacks   

    gama = 10^(SNR(SNR_counter)/10);  %gama deðiþmiyor (shannon, kapasite)
    threshold = (2^R - 1) / gama ;           % R=1 ?
    h_Rician=abs( s*randn(1,Edge_num) + mu ) + 1i*( s*randn(1,Edge_num) + mu ); 
    Edge_con = abs(h_Rician).^2 > threshold ;      
     
     [ii,jj] = ndgrid(1:nodes);              % Used to choose the upper triangle of The Matrix
      A = zeros(nodes);                       % Adjacency Matrix's for all
      A(jj>ii) =  Edge_con;                   % Filling the upper Triangle
      
      
      A = A + A';                             % Adjacency Matrix created
      D = diag(sum(A));                       % In-Degree Matrix created  
      L = D - A ;                             % Laplacian Matrix
           
      L_Networks(monte,:,:) = L ; %her bir monte deðiþkeni birer L A ve D oluþturacak

     
      Ranks(monte,1) = rank(L) ;  %her bir monte deðiþkeni için matrix rank toplamaca
      
      %% Now let's start iterating for every matrix in monte carlo
      
      if rank(L) == nodes - 1 %bu saðlanmýyorsa alphaya ulaþýlamadý ranks5 = 0
     
      Ranks(monte,4) = 1; %we can reach alpha,consensus
      monte_counter = monte_counter + 1;
      
      for k = 1:MaxIt

      P_epsilon = eye(nodes) - epsilon * L;
      IterMeasured = P_epsilon * IterMeasured;  
      StateContainer(1:end,k) = IterMeasured(1:end);
      %% Defining an Error Rule
      if max(IterMeasured) - min(IterMeasured) > Allowed_Error 
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
      Ranks(monte,5) = 0; %bunun bu þekilde olmamasý gerekiyor. iterasyonu önle
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