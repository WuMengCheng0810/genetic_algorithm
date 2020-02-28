%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%程式名稱-基因演算法
%功用-搜尋最佳解
%作者-吳孟城/D1054162010
%變數說明-
%         M:幾條基因(必須是偶數，才能交配) 
%         N:幾位元(當num_var>2時，N必須是偶數) 
%         num_var:幾個變數 
%         pm:突變率
%         move_value:每次移動的單位 
%         range:目標函數範圍
%         generation:代數
%         point_x基因的X座標
%         point_y基因的Y座標
%版本及說明-V5
%               複製是複製適應值較好的前50%，這50%乘2是下一代
%               交配是讓適應值好的先選，選的對象是與它最近的
%與V4版差異-初始座標改為隨機產生，每代座標會移動
%撰寫時間-2018/12/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;  clc;
ans=1; 
M=100;
N=8;
num_var=2; %改變數時:range、fitness要換
precision=[4 4]; %精度，第一個元素代表第一個變數用幾位元表示,總合必須=Ｎ.
pm=0.08;
move_value=1;
range=[0 10;0 10]; 
generation=1;
gene=round(rand(M,N));
gene_new=gene;
p=1;%for 存每代最好的fitness值，以plot
if mod(N,num_var)>0
    disp('N無法整除num_var，請改變N或num_var的大小')
end
point_x=range(1,2)*rand(1,M)-range(1,1);
point_y=range(1,2)*rand(1,M)-range(1,1);
for generation=2:100 %複製交配突變做100次
fitness=[];
binary_gene=cell(M,num_var);
if num_var==1 %一個變數時轉換十進制、算fitness
    for i=1:M
    integer=polyval(gene_new(i,:),2);
    decimal(i,:)=integer*((range(1,2)-range(1,1))/(2^N-1))+range(1,1);
    fitness(i,:)=decimal(i,:)*decimal(i,:);
    end
else %2個變數以上時轉換十進制、算fitness
%fit是一行一行算，'每一行的第一個變數的十進位'跟'每一行的第二個以上變數的十進位'是分開轉換
    for i=1:M %括號內是數字1，此for loop是第一個變數轉十進制
        integer=polyval(gene_new(i,1:precision),2); 
        decimal=integer*((range(1,2)-range(1,1))/(2^(N/2)-1))+range(1,1);
        variable(i,1)=decimal;
        for k=2:num_var % 此for loop是第二個到第num_var個變數轉十進制
            precision_add=0;
            for u=1:k
                precision_add=precision_add+precision(u);
            end
            disp(precision_add);
            binary_gene(i,k)=gene_new(i,(precision_add+1:precision_add+precision(k)));
            integer=polyval(binary_gene(i,k),2);
            decimal=integer*((range(k,2)-range(k,1))/(2^(N/2)-1))+range(k,1);
            variable(i,k)=decimal;
        end
        %sum=0;
        %for j=1:num_var
        %    sum=sum+variable(i,j)^2;
        %   fitness(i)=sum;
        %end
    end
    %fitness在這一行算的話 x是一個矩陣，代表每個變數的十進位
end

fitness=evaluate_fiteness(variable);
%---------reproduction---------
[sort_fitness,sort_fitness_index]=sort(fitness,'descend'); %降冪排列
gene=gene_new;
elite(1,:)=gene(sort_fitness_index(1),:);%存最好的兩個基因(菁英)
elite(2,:)=gene(sort_fitness_index(2),:);
gene_new=[];
k=1;
for i=1:M/2  %把好的50%複製到交配池前50%
    gene_new(k,:)=gene(sort_fitness_index(i),:);
    k=k+1;
end
for i=1:M/2  %把好的50%複製到交配池後50%
    gene_new(k,:)=gene(sort_fitness_index(i),:);
    k=k+1;
end
%-----存每代最好的fitness值-------
best_fitness_index=sort_fitness_index(1);
best_fitness(p,:)=fitness(best_fitness_index);
p=p+1;

%-------------crossover------------
min_distance=10^9;
ready_crossover_first=0;
ready_crossover_second=0;
crossover_ed_second=0;
gene=gene_new; 
crossovered=0;
was_crossover=zeros(1,100);
count=0;
for i=1:M %移動座標並交配
    for k=1:M
        xi_crossovered=ismember(was_crossover,point_x(sort_fitness_index(i)));
        xk_crossovered=ismember(was_crossover,point_x(k));
        yi_crossovered=ismember(was_crossover,point_y(sort_fitness_index(i)));
        yk_crossovered=ismember(was_crossover,point_y(k));
        if i~=k & xi_crossovered~=1 & xk_crossovered~=1 & yi_crossovered~=1 & yk_crossovered~=1 %如果沒交配過才算距離，看是否交配
            subtract_x=point_x(sort_fitness_index(i))-point_x(k);
            subtract_y=point_y(sort_fitness_index(i))-point_y(k);
            distance=sqrt(subtract_x^2+subtract_y^2);
            if distance<min_distance
                min_distance=distance;
                ready_crossover_first=i;
                ready_crossover_second=k;
            end %if
        end %if
    end %for k
    count=count+1;
    was_crossover(count*2-1)=i;
    was_crossover(count*2)=k;
    min_distance=10^9;
    dad_new=[]; %建立一個空矩陣for要當父親的基因
    mother_new=[]; %建立一個空矩陣for要當母親的基因
    k=round(rand()*N); %交配點
    if k<1
       k=k+1;
    end
    gene_new=[];
    child_1=[];
    child_2=[];
    dad=gene(ready_crossover_first,:); %把第一條的基因來當爸
    mother=gene(ready_crossover_second,:); %把第二條的基因來當媽
    for j=1:2:M
        child_1(1:k)=dad(1:k);
        child_1(k+1:N)=mother(k+1:N);
        child_2(1:k)=mother(1:k);
        child_2(k+1:N)=dad(k+1:N);
        gene_new(j,:)=child_1;
        gene_new(j+1,:)=child_2;
    end %for j
end %for i
%移動座標
decide_move=rand(M,2); %第一行>0.4代表X移動，第二行>0.4代表Y移動
move_x=rand(M,2); %第一行>0.4代表X加move_value，第二行>0.4代表X減move_value
move_y=rand(M,2); %第一行>0.4代表Y加move_value，第二行>0.4代表Y減move_value
for i=1:M
    if decide_move(i,1)>0.4 %X移動
        if point_x(i)>range(1,1)+move_value & point_x(i)<range(1,2)-move_value %如果移動不超出範圍才移動
            if move_x(i,1)>0.4 % x+move_value
                point_x(i)=point_x(i)+move_value;
            elseif move_x(i,2)>0.4 %X減move_value
                    point_x(i)=point_x(i)-move_value;
            else
                break
            end %if move_x
        end%if 沒有超出範圍
    elseif decide_move(i,2)>0.4 %y移動
        if point_y(i)>range(2,1)+move_value & point_y(i)<range(2,2)-move_value %如果移動不超出範圍才移動
            if move_y(i,1)>0.4 % y+move_value
                point_y(i)=point_y(i)+move_value;
            elseif move_y(i,2)>0.4 %y減move_value
                    point_y(i)=point_y(i)-move_value;
            else
                break
            end %if move_y
        end%if 沒有超出範圍
    else
        break
    end
end%for i
%---------mutation---------
gene=gene_new;
mutation_gene=0;
   for i=1:M %看每一條gene是否要突變
        ready_mutation=gene(i,:); %可能會突變的那一條是gene_old
        k=rand(); %產生一個0~0.8的數字，看是否突變
         if  k<pm %如果達到突變率(小於突變率):突變
             mutation_point=ceil(rand()*N); %取得突變點(0~n的數字)
             %開始突變
             if mutation_point==0
                 mutation_point=mutation_point+round(rand()*N);
             end%if 
             a=~ready_mutation(mutation_point); %把突變點的1變0，0變1
             gene_new(i,mutation_point)=a; %將突變後的值(0或1)放進突變後的第j個，也就是說新的突變點覆蓋掉舊的突變點(0變1or1變0)
             mutation_gene=mutation_gene+1;
         end;
    end;
    gene_new(1,:)=elite(1,:);
    gene_new(2,:)=elite(2,:);
    end;  %for generation 
    
disp('generation is')
disp(generation);
disp('突變的基因有')
disp(mutation_gene)
disp('最好的適應值')
disp(best_fitness)
plot(best_fitness)
disp('最好的基因')
disp(elite)
disp('基因數量(條)=')
disp(M)

