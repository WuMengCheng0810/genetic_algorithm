%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�{���W��-��]�t��k
%�\��-�j�M�̨θ�
%�@��-�d�s��/D1054162010
%�ܼƻ���-
%         M:�X����](�����O���ơA�~���t) 
%         N:�X�줸(��num_var>2�ɡAN�����O����) 
%         num_var:�X���ܼ� 
%         pm:���ܲv
%         move_value:�C�����ʪ���� 
%         range:�ؼШ�ƽd��
%         generation:�N��
%         point_x��]��X�y��
%         point_y��]��Y�y��
%�����λ���-V5
%               �ƻs�O�ƻs�A���ȸ��n���e50%�A�o50%��2�O�U�@�N
%               ��t�O���A���Ȧn������A�諸��H�O�P���̪�
%�PV4���t��-��l�y�Чאּ�H�����͡A�C�N�y�з|����
%���g�ɶ�-2018/12/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;  clc;
ans=1; 
M=100;
N=8;
num_var=2; %���ܼƮ�:range�Bfitness�n��
precision=[4 4]; %��סA�Ĥ@�Ӥ����N��Ĥ@���ܼƥδX�줸���,�`�X����=��.
pm=0.08;
move_value=1;
range=[0 10;0 10]; 
generation=1;
gene=round(rand(M,N));
gene_new=gene;
p=1;%for �s�C�N�̦n��fitness�ȡA�Hplot
if mod(N,num_var)>0
    disp('N�L�k�㰣num_var�A�Ч���N��num_var���j�p')
end
point_x=range(1,2)*rand(1,M)-range(1,1);
point_y=range(1,2)*rand(1,M)-range(1,1);
for generation=2:100 %�ƻs��t���ܰ�100��
fitness=[];
binary_gene=cell(M,num_var);
if num_var==1 %�@���ܼƮ��ഫ�Q�i��B��fitness
    for i=1:M
    integer=polyval(gene_new(i,:),2);
    decimal(i,:)=integer*((range(1,2)-range(1,1))/(2^N-1))+range(1,1);
    fitness(i,:)=decimal(i,:)*decimal(i,:);
    end
else %2���ܼƥH�W���ഫ�Q�i��B��fitness
%fit�O�@��@���A'�C�@�檺�Ĥ@���ܼƪ��Q�i��'��'�C�@�檺�ĤG�ӥH�W�ܼƪ��Q�i��'�O���}�ഫ
    for i=1:M %�A�����O�Ʀr1�A��for loop�O�Ĥ@���ܼ���Q�i��
        integer=polyval(gene_new(i,1:precision),2); 
        decimal=integer*((range(1,2)-range(1,1))/(2^(N/2)-1))+range(1,1);
        variable(i,1)=decimal;
        for k=2:num_var % ��for loop�O�ĤG�Ө��num_var���ܼ���Q�i��
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
    %fitness�b�o�@��⪺�� x�O�@�ӯx�}�A�N��C���ܼƪ��Q�i��
end

fitness=evaluate_fiteness(variable);
%---------reproduction---------
[sort_fitness,sort_fitness_index]=sort(fitness,'descend'); %�����ƦC
gene=gene_new;
elite(1,:)=gene(sort_fitness_index(1),:);%�s�̦n����Ӱ�](�׭^)
elite(2,:)=gene(sort_fitness_index(2),:);
gene_new=[];
k=1;
for i=1:M/2  %��n��50%�ƻs���t���e50%
    gene_new(k,:)=gene(sort_fitness_index(i),:);
    k=k+1;
end
for i=1:M/2  %��n��50%�ƻs���t����50%
    gene_new(k,:)=gene(sort_fitness_index(i),:);
    k=k+1;
end
%-----�s�C�N�̦n��fitness��-------
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
for i=1:M %���ʮy�Шå�t
    for k=1:M
        xi_crossovered=ismember(was_crossover,point_x(sort_fitness_index(i)));
        xk_crossovered=ismember(was_crossover,point_x(k));
        yi_crossovered=ismember(was_crossover,point_y(sort_fitness_index(i)));
        yk_crossovered=ismember(was_crossover,point_y(k));
        if i~=k & xi_crossovered~=1 & xk_crossovered~=1 & yi_crossovered~=1 & yk_crossovered~=1 %�p�G�S��t�L�~��Z���A�ݬO�_��t
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
    dad_new=[]; %�إߤ@�Ӫůx�}for�n����˪���]
    mother_new=[]; %�إߤ@�Ӫůx�}for�n����˪���]
    k=round(rand()*N); %��t�I
    if k<1
       k=k+1;
    end
    gene_new=[];
    child_1=[];
    child_2=[];
    dad=gene(ready_crossover_first,:); %��Ĥ@������]�ӷ�
    mother=gene(ready_crossover_second,:); %��ĤG������]�ӷ��
    for j=1:2:M
        child_1(1:k)=dad(1:k);
        child_1(k+1:N)=mother(k+1:N);
        child_2(1:k)=mother(1:k);
        child_2(k+1:N)=dad(k+1:N);
        gene_new(j,:)=child_1;
        gene_new(j+1,:)=child_2;
    end %for j
end %for i
%���ʮy��
decide_move=rand(M,2); %�Ĥ@��>0.4�N��X���ʡA�ĤG��>0.4�N��Y����
move_x=rand(M,2); %�Ĥ@��>0.4�N��X�[move_value�A�ĤG��>0.4�N��X��move_value
move_y=rand(M,2); %�Ĥ@��>0.4�N��Y�[move_value�A�ĤG��>0.4�N��Y��move_value
for i=1:M
    if decide_move(i,1)>0.4 %X����
        if point_x(i)>range(1,1)+move_value & point_x(i)<range(1,2)-move_value %�p�G���ʤ��W�X�d��~����
            if move_x(i,1)>0.4 % x+move_value
                point_x(i)=point_x(i)+move_value;
            elseif move_x(i,2)>0.4 %X��move_value
                    point_x(i)=point_x(i)-move_value;
            else
                break
            end %if move_x
        end%if �S���W�X�d��
    elseif decide_move(i,2)>0.4 %y����
        if point_y(i)>range(2,1)+move_value & point_y(i)<range(2,2)-move_value %�p�G���ʤ��W�X�d��~����
            if move_y(i,1)>0.4 % y+move_value
                point_y(i)=point_y(i)+move_value;
            elseif move_y(i,2)>0.4 %y��move_value
                    point_y(i)=point_y(i)-move_value;
            else
                break
            end %if move_y
        end%if �S���W�X�d��
    else
        break
    end
end%for i
%---------mutation---------
gene=gene_new;
mutation_gene=0;
   for i=1:M %�ݨC�@��gene�O�_�n����
        ready_mutation=gene(i,:); %�i��|���ܪ����@���Ogene_old
        k=rand(); %���ͤ@��0~0.8���Ʀr�A�ݬO�_����
         if  k<pm %�p�G�F����ܲv(�p����ܲv):����
             mutation_point=ceil(rand()*N); %���o�����I(0~n���Ʀr)
             %�}�l����
             if mutation_point==0
                 mutation_point=mutation_point+round(rand()*N);
             end%if 
             a=~ready_mutation(mutation_point); %������I��1��0�A0��1
             gene_new(i,mutation_point)=a; %�N���ܫ᪺��(0��1)��i���ܫ᪺��j�ӡA�]�N�O���s�������I�л\���ª������I(0��1or1��0)
             mutation_gene=mutation_gene+1;
         end;
    end;
    gene_new(1,:)=elite(1,:);
    gene_new(2,:)=elite(2,:);
    end;  %for generation 
    
disp('generation is')
disp(generation);
disp('���ܪ���]��')
disp(mutation_gene)
disp('�̦n���A����')
disp(best_fitness)
plot(best_fitness)
disp('�̦n����]')
disp(elite)
disp('��]�ƶq(��)=')
disp(M)

