% Step size sabit
%Nakagami deneme
%1-prob + semilog

clear; clc; close all;dbstop if error;
%% Defining the System
p = 1; %Signal power 
MaxIt = 2000;
SNR = [0:1:20];  
epsilon = 1/100;  %step size determined random
f1 = figure;
f2 = figure;
min_node = 3;
max_node = 4;

%% This is for naming the data
number = [min_node:1:max_node]';
names = int2str(number);
                                
Allowed_Error = 1/100; 
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
montemax = 1e3;                           %how many times monte carlo                         
L_Networks = zeros(montemax,nodes,nodes); %thil is L matrix container
Ranks = zeros(montemax,6);                %this will keep everything about matrices
for SNR_counter = 1:length(SNR) 
    tic %alttaki yazı 20 defa tekrarlanacak.Word'e kopyala ve onları say. Kalan zamanı hesapla.
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
%%%%mean(h_Naka.^2);  kullanmak gerekir mi?
gama = 10^(SNR(SNR_counter)/10); 
threshold = (2^R - 1) / gama ;
 Edge_con = abs(h_Naka).^2 > threshold ;

StateContainer = zeros(size(Measured,1), MaxIt + 1 ); %for every monte, this will keep states
     
     StateContainer(:,1) = Measured(:); %içine alamıyor tranpoz sorunu
     
     FeedbackContainer = zeros(MaxIt,size(Measured,2));   %for every monte this will keep feedbacks


      
     
     [ii,jj] = ndgrid(1:nodes);              % Used to choose the upper triangle of The Matrix
      A = zeros(nodes);                       % Adjacency Matrix's for all
      A(jj>ii) =  Edge_con;                   % Filling the upper Triangle
     
      
      A = A + A';                             % Adjacency Matrix created
      D = diag(sum(A));                       % In-Degree Matrix created  
      L = D - A ;                             % Laplacian Matrix
           
      L_Networks(monte,:,:) = L ; %her bir monte değişkeni birer L A ve D oluşturacak

     
      Ranks(monte,1) = rank(L) ;  %her bir monte değişkeni için matrix rank toplamaca
      
      %% Now let's start iterating for every matrix in monte carlo
      
      if rank(L) == nodes - 1 %bu sağlanmıyorsa alphaya ulaşılamadı ranks5 = 0
     
      Ranks(monte,4) = 1; %we can reach alpha,consensus
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
      Ranks(monte,5) = 0; %bunun bu şekilde olmaması gerekiyor. iterasyonu önle
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
