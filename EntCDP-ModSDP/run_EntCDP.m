clc;clear
load('.\example\ALL_ST.mat');
load('.\example\AML_TARGET.mat');
A{1}=ALL_ST;
A{2}=AML_TARGET;

[~,gene_type]=xlsread('.\drivergene.xlsx');
[B,mg,gene_type_new]=meta(A,gene_type);

gene_type_new0=gene_type_new;
for k=2:10
    k
    gene_type_new=gene_type_new0;
    
    num=length(B);
    exclusion=[];
    if exclusion>=1
        gene_type_new(exclusion)=[];
    end
    
    max_value=1;
    for t=1:3
        max_geneset0=EntCDP_matlab(B,k,exclusion);
    
        if (max(max_geneset0(:,k+1))>max_value) || (max(max_geneset0(:,k+1))==max_value && length(find(max_geneset0(:,k+1)==max_value))>length(find(max_geneset(:,k+1)==max_value)))
            max_value=max(max_geneset0(:,k+1));
            max_geneset=max_geneset0;
        end
    end
    m0=find(max_geneset(:,k+1)==max_value);

    com_gene=cell(1,k);
    for j=1:k
        com_gene{j}=gene_type_new{max_geneset(m0(1),j)};
        for u=1:size(mg,1)
            s=0;
            for v=1:size(mg,2)
                if isempty(mg{u,v})==1
                    break;
                end
                s=s+1;
            end
            %find(categorical(mg(u,1:s))==com_gene{j})
            if length(find(categorical(mg(u,1:s))==com_gene{j}))==1
                disp('meta')
                disp(com_gene{j})
                break
            end
        end
    end

    if sum(max_geneset(m0(1),k+2:k+2+num)<0.05)==num+1
        disp(com_gene);
        %disp(max_geneset(i,1:k));disp(max_geneset(m0(i),k+1));
        disp(max_geneset(m0(1),k+2:k+2+num));
        
    else
        disp('Not significant');
    end
    clear com_gene

    
    clear max_geneset com_gene;
end

