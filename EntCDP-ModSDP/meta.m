function [B,metagene,gene_type_new]=meta(A,gene_type)
%B=A;
num=length(A);
for r=1:num
    [m(r),~]=size(A{r});
end
%s=zeros(num,n);
%u=zeros(1,n);
%for r=1:num
s=sum(A{1},1);

u=unique(s);

t=zeros(100,100);
A1=A;
zhou=1;
for i=1:length(u)
    w=find(s==u(i));
    if length(w)==1
        continue;
    else
        for j=1:length(w)-1
            wen=1;
            for l=j+1:length(w)
                if isequal(A{1}(:,w(j)),A{1}(:,w(l)))==1  && w(j)~=w(l)
                    %[w(j),w(l)]
                    tep=1;
                    for r=2:num
                        if isequal(A{r}(:,w(j)),A{r}(:,w(l)))==0
                            break;
                        else
                            %[w(j),w(l)]
                            tep=tep+1;
                        end
                    end
                    if tep==num
                        t(zhou,wen*2-1:wen*2)=[w(j),w(l)];
                        A1{1}(:,w(l))=ones(m(1),1);
                        if i<length(u)
                            ww=find(s==u(i+1));
                            w(l)=ww(1);
                        else 
                        	ww=find(s==u(1));
                            w(l)=ww(1);
                        end
                        wen=wen+1;
                    end
                end
            end
            if wen>1
                zhou=zhou+1;
            end
        end
    end
end
for r=2:num
    k=1:num;
    k(find(k==r))=[];
    s=sum(A{r},1);
    w=find(s==0);
    if length(w)==1
        continue;
    else
        for j=1:length(w)-1
            wen=1;
            for l=j+1:length(w)
                if isequal(A{r}(:,w(j)),A{1}(:,w(l)))==1  && w(j)~=w(l)
                    %[w(j),w(l)]
                    tep=1;
                    for h=1:num-1

                        if isequal(A{k(h)}(:,w(j)),A{k(h)}(:,w(l)))==0
                            break;
                        else
                            %[w(j),w(l)]
                            tep=tep+1;
                        end
                    end
                    if tep==num
                        t(zhou,wen*2-1:wen*2)=[w(j),w(l)];
                        A1{r}(:,w(l))=ones(m(r),1);
                        ww=find(s==1);
                        w(l)=ww(1);
                        wen=wen+1;
                    end
                end
            end
            if wen>1
                zhou=zhou+1;
            end
        end
    end
end
ss=[];
for r=1:num
    ss=[ss,find(sum(A1{r},1)==m(r))];
end
su=unique(ss);

for r=1:num
    A1{r}(:,su)=[];
end
gene_type_new=gene_type;
gene_type_new(su)=[];
B=A1;
metagene=[];
for i=1:zhou-1
    v=unique(t(i,:));
    v(find(v==0))=[];
    for j=1:length(v)
        metagene{i,j}=gene_type{v(j)};
    end
end

                    
            