function max_geneset=ModSDP_matlab(A,C,k,exclusion)
%
% This function is designed to identify specific driver gene sets in cancer (A relative to C).
%
% This function uses the MATLAB interface of CPLEX to solve the Binary Linear Program (BLP) problem. 
% Before running this funciton, please make sure you have installed CPLEX software in your computer.
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
% max_geneset : output vector with k+r+s+3 columns (r and s denote the numbers of corresponding cancer types considered);
%               it records the information of the optimal solution; 
%               max_geneset(1:k) records the selected genes;
%               max_geneset(k+1) records the weight of this gene set;
%               max_geneset(k+2:k+r+s+1) records the significance level of this gene set in each cancer type;
%               max_geneset(k+r+s+2) records the overall significance level of this gene set;
%               max_geneset(k+r+s+3) records the number of selected genes (i.e., k);
% 

num1=length(A);
for i=1:num1
    if ~isempty(exclusion)
        A{i}(:,exclusion)=[];
    end
    [m(i),n]=size(A{i});
end

num2=length(C);
for i=1:num2
    if ~isempty(exclusion)
        C{i}(:,exclusion)=[];
    end
    [t(i),n]=size(C{i});
end

f=zeros(1,sum(m)+sum(t)+(num1+num2)*n);   % Vector containing the coefficients of the linear objective function
tem1=0;
for i=1:num1
    %f(tem1+1:tem1+m(i))=-4;
    %f(tem1+1:tem1+m(i))=-k;   %%%
    f(tem1+1:tem1+m(i))=-2/num1;   %%%
    %f(tem1+m(i)+1:tem1+m(i)+n)=2*sum(A{i});
    %f(tem1+m(i)+1:tem1+m(i)+n)=sum(A{i});
    f(tem1+m(i)+1:tem1+m(i)+n)=sum(A{i})/num1;
    tem1=tem1+m(i)+n;
end
common_s=sum(m)+num1*n;
tem2=0;
for i=1:num2
    %f(common_s+tem2+1:common_s+tem2+t(i))=k;   %%%
    f(common_s+tem2+1:common_s+tem2+t(i))=k/num2/(k-1);   %%%
    %f(common_s+tem2+t(i)+1:common_s+tem2+t(i)+n)=-sum(C{i});
    f(common_s+tem2+t(i)+1:common_s+tem2+t(i)+n)=-sum(C{i})/num2/(k-1);
    tem2=tem2+t(i)+n;
end

clear tem1; clear tem2;

B=zeros(sum(m)+2*sum(t),sum(m)+sum(t)+(num1+num2)*n);  
%B=zeros(sum(m)+2*sum(t)+num1+num2,sum(m)+sum(t)+(num1+num2)*n);  % Matrix containing the coefficients of the linear inequality constraints
   
tem1=0; tem2=0;
for i=1:num1
    B(tem1+1:tem1+m(i),tem2+1:tem2+m(i))=diag(ones(1,m(i)));%x^(i)_1~mi
    B(tem1+1:tem1+m(i),tem2+m(i)+1:tem2+m(i)+n)=-A{i};
    tem1=tem1+m(i);
    tem2=tem2+m(i)+n;
end
clear tem1; clear tem2;
tem1=0; tem2=0;
for i=1:num2
    B(sum(m)+tem1+1:sum(m)+tem1+t(i),common_s+tem2+1:common_s+tem2+t(i))=diag(ones(1,t(i)));
    B(sum(m)+tem1+1:sum(m)+tem1+t(i),common_s+tem2+t(i)+1:common_s+tem2+t(i)+n)=-C{i};
    B(sum(m)+tem1+t(i)+1:sum(m)+tem1+2*t(i),common_s+tem2+1:common_s+tem2+t(i))=-n*diag(ones(1,t(i)));   %%%
    B(sum(m)+tem1+t(i)+1:sum(m)+tem1+2*t(i),common_s+tem2+t(i)+1:common_s+tem2+t(i)+n)=C{i};
    tem1=tem1+2*t(i);
    tem2=tem2+t(i)+n;
end
clear tem1; clear tem2;

% tem1=0; tem2=0;
% for i=1:num1
%     B(sum(m)+2*sum(t)+tem1+1,tem2+1:tem2+m(i))=-2*(ones(1,m(i)));
%     B(sum(m)+2*sum(t)+tem1+1,tem2+t(i)+1:tem2+t(i)+n)=sum(A{i});
%     tem1=tem1+1;
%     tem2=tem2+m(i)+n;
% end

% clear tem1; clear tem2;

% tem1=0; tem2=0;
% for i=1:num2
%     B(sum(m)+2*sum(t)+num1+tem1+1,common_s+tem2+1:common_s+tem2+t(i))=-2*(ones(1,t(i)));
%     B(sum(m)+2*sum(t)+num1+tem1+1,common_s+tem2+t(i)+1:common_s+tem2+t(i)+n)=sum(C{i});
%     tem1=tem1+1;
%     tem2=tem2+t(i)+n;
% end
% clear tem1; clear tem2;

b=zeros(sum(m)+2*sum(t),1);
%b=zeros(sum(m)+2*sum(t)+num1+num2,1);     % Vector corresponding to the right-hand side of the linear inequality constraints

Aeq=zeros((num1+num2-1)*n+1,sum(m)+sum(t)+(num1+num2)*n);     % Matrix containing the coefficients of the linear equality constraints
Aeq(1,m(1)+1:m(1)+n)=1;
tem1=0; tem2=0;
for i=1:num1-1
    Aeq(tem1+2:tem1+n+1,m(1)+1:m(1)+n)=diag(ones(1,n));
    Aeq(tem1+2:tem1+n+1,tem2+m(i)+m(i+1)+n+1:tem2+m(i)+m(i+1)+2*n)=-diag(ones(1,n));
    tem1=tem1+n;
    tem2=tem2+m(i)+n;
end
clear tem1; clear tem2;
equal_f=(num1-1)*n+1;
tem1=0;
tem2=0;
for i=1:num2
    Aeq(equal_f+tem1+1:equal_f+tem1+n,m(1)+1:m(1)+n)=diag(ones(1,n));
    Aeq(equal_f+tem1+1:equal_f+tem1+n,common_s+tem2+t(i)+1:common_s+tem2+t(i)+n)=-diag(ones(1,n));
    tem1=tem1+n;
    tem2=tem2+t(i)+n;
end
clear tem1; clear tem2;
beq=zeros((num1+num2-1)*n+1,1);       % Vector containing the constants of the linear equality constraints
beq(1)=k;
% size(f)
% size(B)
% size(b)
[X,maxweight]= cplexbilp(f,B,b,Aeq,beq);

k=sum(X(m(1)+1:m(1)+n)>0.5);
max_geneset(1:k)=find(X(m(1)+1:m(1)+n)>0.5);
max_geneset(k+1)=-maxweight;
for i=1:num1
    max_geneset(k+i+1)=significance_A(A{i},max_geneset(1:k));   %%%
end
for i=1:num2
    max_geneset(k+num1+i+1)=significance_NA(C{i},max_geneset(1:k));   %%%
end
max_geneset(k+num1+num2+2)=significance_specific(A,C,max_geneset(1:k),k);   %%%
max_geneset(k+num1+num2+3)=k;


function p=significance_A(A,subset)   %%%
%[m,~]=size(A);
[r,~]=size(A);
w=zeros(1,1000);
%n=length(subset);
t=length(subset);
for j=1:1000
    A_temp=A;
    A_temp(:,subset)=0;
    for i=1:t
        temp=sum(A(:,subset(i)));
        index = randperm(r,temp);
        A_temp(index,subset(i))=1;
    end
    w(j)=fit_A(A_temp,subset);   %%%
end
p=sum(w>=fit_A(A,subset))/1000;   %%%

function p=significance_NA(A,subset)   %%%
%[m,~]=size(A);
[r,~]=size(A);
w=zeros(1,1000);
%n=length(subset);
t=length(subset);
for j=1:1000
    A_temp=A;
    A_temp(:,subset)=0;
    for i=1:t
        temp=sum(A(:,subset(i)));
        index = randperm(r,temp);
        A_temp(index,subset(i))=1;
    end
    w(j)=fit_A(A_temp,subset);   %%%
end
%p=sum((w<fit_A(A,subset))&(w>=0))/1000;   %%%
p=sum(w>=fit_A(A,subset))/1000;   %%%%%%

function p=significance_specific(A,C,subset,coef)   %%%
mul_1=length(A);
mul_2=length(C);
for k=1:mul_1
    [r(k),~]=size(A{k});
end
for k=1:mul_2
    [s(k),~]=size(C{k});
end
t=length(subset);
w=zeros(1,1000);
for j=1:1000
    for k=1:mul_1
        A_temp{k}=A{k};
        A_temp{k}(:,subset)=0;
        for i=1:t
            temp=sum(A{k}(:,subset(i)));
            index = randperm(r(k),temp);
            A_temp{k}(index,subset(i))=1;
        end
    end
    for k=1:mul_2
        C_temp{k}=C{k};
        C_temp{k}(:,subset)=0;
        for i=1:t
            temp=sum(C{k}(:,subset(i)));
            index = randperm(s(k),temp);
            C_temp{k}(index,subset(i))=1;
        end
    end    
    w(j)=fit_specific(A_temp,C_temp,subset,coef);   %%%
end
p=sum(w>=fit_specific(A,C,subset,coef))/1000;   %%%


function q=significance_pair(A,subset1,subset2)
[m,~]=size(A);
v=zeros(1,1000);
subset=[subset1 subset2];
n=length(subset);
for j=1:1000
    A_temp=A;
    A_temp(:,subset)=0;
    for i=1:n
        temp=sum(A(:,subset(i)));
        index = randperm(m,temp);
        A_temp(index,subset(i))=1;
    end
    e=fit_pair(A_temp,subset1,subset2);
    v(j)=e(5);
    clear e;
end
h=fit_pair(A,subset1,subset2);
u=h(5);
clear h;
q=sum(v>=u)/1000;


function f=fit_A(A,x)        % fitness calculation function
[~,n]=size(A);
temp=zeros(1,n);
temp(x)=1;
x=temp;
index = x==1;
A_index=A(:,index);
A_index_sum=sum(A_index,2);
f=(2*sum(A_index_sum>0)-sum(A_index_sum));   %%%

function f=fit_specific(A,C,x,coef)   %%%
mul_1=length(A);
mul_2=length(C);
[~,n]=size(A{1});
temp=zeros(1,n);
temp(x)=1;
x=temp;
index = x==1;
f=0;
for i=1:mul_1
    A_index{i}=A{i}(:,index);
    A_index_sum{i}=sum(A_index{i},2);
    g(i)=(2*sum(A_index_sum{i}>0)-sum(A_index_sum{i}));   %%%
end
for i=1:mul_2
    C_index{i}=C{i}(:,index);
    C_index_sum{i}=sum(C_index{i},2);
    h(i)=(coef*sum(C_index_sum{i}>0)-sum(C_index_sum{i}))/coef;   %%%
end
%f=sum(g)-sum(h);
f=sum(g)/mul_1-sum(h)/mul_2;   %%%


function g=fit_pair(A,x,y)
[~,n]=size(A);
temp1=zeros(1,n);
temp2=zeros(1,n);
temp1(x)=1;
temp2(y)=1;
x=temp1;
y=temp2;
index1 = x==1;
index2 = y==1;
A_index1=A(:,index1);
A_index2=A(:,index2);
A_index1_sum=sum(A_index1,2);
A_index2_sum=sum(A_index2,2);
intersec = (A_index1_sum>0)&(A_index2_sum>0);
union = (A_index1_sum>0)|(A_index2_sum>0);
g(1)=sum(A_index1_sum>0);
g(2)=sum(A_index2_sum>0);
g(3)=sum(intersec);
g(4)=sum(union);
g(5)=sum(intersec)/sum(union);
