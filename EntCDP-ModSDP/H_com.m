function w=H_com(A,x,exclusion)
num=length(A);
if ~isempty(exclusion)
    for r=1:num
        A{r}(:,exclusion)=[];
    end
    [~,n]=size(A{1});
end
n=size(A{1},2);
%x=find(popuCode==1);
temp=zeros(1,n);
temp(x)=1;
x=temp;
index = x==1;
%f=0;
A_index=cell(1,num);
A_index_sum=cell(1,num);
w=zeros(1,num);
f=0;
for i=1:num
    m=size(A{i},1);
    A_index{i}=A{i}(:,index);
    A_index_sum{i}=sum(A_index{i},2);
    w(i)=(sum(A_index_sum{i}>0))/m;
end
