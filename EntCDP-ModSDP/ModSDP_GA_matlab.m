function maxpop = ModSDP_GA_matlab(A,C,k,exclusion)
%
%If IBM ILOG CPLEX optimizer is not added to your matlab, genetic algorithms can be used instead
% This function uses genetic algorithm to identify specific driver gene sets in cancer (A relative to C).
%
% A and C are all mutation matrices with the same number of columns: m_i (sample) x n (gene).
% For example, to identify specific driver gene sets of certain r cancer types relative to other s cancer types, A ans C can be represented as follows: 
% A : A{1}, A{2}, ..., A{r}
% C : C{1}, C{2}, ..., C{s}
%
% k : number of desired genes 
%
% exclusion : the genes which are to be excluded from the results
%
% maxpop: output vector with k+r+s+3 columns (r and s denote the numbers of corresponding cancer types considered);
%         it records the information of the optimal solution; 
%         maxpop(1:k) records the selected genes;
%         maxpop(k+1) records the weight of this gene set;
%         maxpop(k+2:k+r+s+1) records the significance level of this gene set in each cancer type;
%         maxpop(k+r+s+2) records the overall significance level of this gene set;
%         maxpop(k+r+s+3) records the number of selected genes (i.e., k);
% 

num=length(A);
num2=length(C);


if length(exclusion)>=1
    for r=1:num
        A{r}(:,exclusion)=[];
    end
    for r=1:num2
        C{r}(:,exclusion)=[];
    end
    [~,n]=size(A{1});
end
[~,n]=size(A{1});
popsize=100*(floor(n/100)); % the number of individual in one population
nger=1000; % maximal number of iterations, can be adjusted by users

object_value=zeros(1,nger);
pop=zeros(popsize*2,k+1);%% % population matrix, every row record an individual

% Generate initial population
for i=1:popsize*2
    temp=randperm(n);
    pop(i,1:k)=temp(1:k);
    pop(i,k+1)=fitmulti(A,C,pop(i,1:k),k);
end
[~,I]=sort(pop(:,k+1),'descend');
pop=pop(I,:);  % sort the population in descending order

R=0; % record whether current optimal solution is trapped in one local optimal solution
i=1;
while(R<10 && i<=nger) % termination
    fit_vector=pop(:,k+1);
    temp_maxweight=max(fit_vector);
    j=1;
    while(j<=popsize)  % generate one new individual
        [index1, index2]=select_order_fitness(fit_vector);
        pop(popsize+j,1:k)=crossover(pop(index1,1:k),pop(index2,1:k),n); % crosscove
        pop(popsize+j,1:k)=mutation_SA(A,C,pop(popsize+j,1:k),n,1,k); % mutation
        pop(popsize+j,k+1)=fitmulti(A,C,pop(popsize+j,1:k),k);
        j=j+1;
     end
    [~,I]=sort(pop(:,k+1),'descend');
    pop=pop(I,:);   % sort the population in descending order, the maximal n individual are transfered to next generation
    object_value(i)=pop(1,k+1);  
    
    
    %use local search to improve the current optimal solution when it is trapped in one local optimal solution
    if R==2
        temp=sort(unique(pop(:,k+1)),'Descend');
        index=find(pop(:,k+1)>=temp(1));
        for j=1:min(length(index),1)
            pop(index(j),1:k)=mutation_SA(A,C,pop(index(j),1:k),n,sqrt(n),k);
            pop(index(j),k+1)=fitmulti(A,C,pop(index(j),1:k),k);
        end
    end
    
    if R==5
        temp=sort(unique(pop(:,k+1)),'Descend');
        index=find(pop(:,k+1)>=temp(1));
        for j=1:min(length(index),1)
            pop(index(j),1:k)=mutation_SA(A,C,pop(index(j),1:k),n,n,k);
            pop(index(j),k+1)=fitmulti(A,C,pop(index(j),1:k),k);
        end
    end
    
    maxweight=max(pop(:,k+1));

    if maxweight==temp_maxweight
        R=R+1;
    else
        R=0;
    end
    
    i=i+1;
end
[~,I]=sort(pop(:,k+1),'descend');
pop=pop(I,:);
maxpop=pop(1:popsize,:);

%delete the repetitive solution
[m,~]=size(maxpop);
i=1;
while(i<m)
    index=zeros(1,m);
    for j=i+1:m
        if all(maxpop(j,:)==maxpop(i,:))
            index(j)=1;
        end
    end
    index=logical(index);
    maxpop(index,:)=[];
    [m,~]=size(maxpop);
    i=i+1;
end

%significance test
[m,~]=size(maxpop);
for i=1:m
    for j=1:num
        maxpop(i,k+j+1)=significance_A(A{j},maxpop(i,1:k));
    end
    for j=1:num2
        maxpop(i,k+num+j+1)=significance_A(C{j},maxpop(i,1:k));
    end
    maxpop(i,k+num+num2+2)=significance_multi(A,C,maxpop(i,1:k),k);
    maxpop(i,k+num+num2+3)=k;
end



 %The funcitons used in the above main function
function f=fitmulti(A,C,x,k)        % fitness calculation function

num1=length(A);
n=size(A{1},2);
temp=zeros(1,n);
temp(x)=1;
x=temp;
index = x==1;
A_index=cell(1,num1);
A_index_sum=cell(1,num1);

num2=length(C);
C_index=cell(1,num1);
C_index_sum=cell(1,num1);

f1=0;f2=0;
for i=1:num1
    A_index{i}=A{i}(:,index);
    A_index_sum{i}=sum(A_index{i},2);
    f1=f1+(2*sum(A_index_sum{i}>0)-sum(A_index_sum{i}));
end
for i=1:num2
    C_index{i}=C{i}(:,index);
    C_index_sum{i}=sum(C_index{i},2);
    f2=f2+(k*sum(C_index_sum{i}>0)-sum(C_index_sum{i}))/k;%/m(i);(min(m)/m(i));
end

f=f1/num1-f2/num2;




function f=fit_A(A,x)        % fitness calculation function

[m,n]=size(A);
temp=zeros(1,n);
temp(x)=1;
x=temp;
index = x==1;
A_index=A(:,index);
A_index_sum=sum(A_index,2);
%w=(sum(A_index_sum>0))/m;%%%%%%%%%%%%%%
f=(2*sum(A_index_sum>0)-sum(A_index_sum))*1;

function [index1, index2]=select_order_fitness(fit_vector)   % select two individuals based on rank fitness
n=length(fit_vector);
[~,I]=sort(fit_vector);
p=zeros(1,n);
for i=1:n
    p(I(i))=2*i/(n*(n+1));
end
pp=cumsum(p); 
random_data=rand(1,1);
temp=find(pp>=random_data);
index1=temp(1);
random_data=rand(1,1);
temp=find(pp>=random_data);
index2=temp(1);


function newpop=crossover(parent1,parent2,n)   % crossover function
temp=zeros(1,n);
temp(parent1)=1;
parent1=temp;
temp=zeros(1,n);
temp(parent2)=1;
parent2=temp;
k=sum(parent1);
newpop=zeros(1,n);
index=(parent1+parent2==2);
newpop(index)=1;
parent1(index)=0;
parent2(index)=0;
temp=find(parent1+parent2==1);
index=randperm(sum(parent1+parent2==1));
newpop(temp(index(1:(k-sum(newpop)))))=1;
newpop=find(newpop==1);


function m_x=mutation(x,n)             % mutation function
temp=zeros(1,n);
temp(x)=1;
x=temp;
k=sum(x);
index1=randi(n,1);
while(x(index1)==1)
    index1=randi(n,1);
end
x_nonzero=find(x==1);
index2=x_nonzero(randi(k,1));
m_x=x;
m_x(index1)=1;
m_x(index2)=0;
m_x=find(m_x==1);


function m_x=mutation_SA(A,C,pop_i,n,N,k)  % local search
for i=1:N
    pop_j=mutation(pop_i,n);
    if fitmulti(A,C,pop_j,k)>=fitmulti(A,C,pop_i,k)
        pop_i=pop_j;
    end
end
m_x=pop_i;


function p=significance_A(A,subset)   % significance test function
[m,~]=size(A);
w=zeros(1,1000);
n=length(subset);
for j=1:1000
    A_temp=A;
    A_temp(:,subset)=0;
    for i=1:n
        temp=sum(A(:,subset(i)));
        index = randperm(m,temp);
        A_temp(index,subset(i))=1;
    end
    w(j)=fit_A(A_temp,subset);
end
p=sum(w>=fit_A(A,subset))/1000;


function p=significance_multi(A,C,subset,K)
mul=length(A);
mul2=length(C);
for k=1:mul
    [r(k),~]=size(A{k});
end
for k=1:mul2
    [s(k),~]=size(C{k});
end
t=length(subset);
w=zeros(1,1000);
for j=1:1000
    for k=1:mul
        A_temp{k}=A{k};
        A_temp{k}(:,subset)=0;
        for i=1:t
            temp=sum(A{k}(:,subset(i)));
            index = randperm(r(k),temp);
            A_temp{k}(index,subset(i))=1;
        end
    end
    for k=1:mul2
        C_temp{k}=C{k};
        C_temp{k}(:,subset)=0;
        for i=1:t
            temp=sum(C{k}(:,subset(i)));
            index = randperm(s(k),temp);
            C_temp{k}(index,subset(i))=1;
        end
    end
    w(j)=fitmulti(A_temp,C_temp,subset,K);
end
p=sum(w>=fitmulti(A,C,subset,K))/1000;
