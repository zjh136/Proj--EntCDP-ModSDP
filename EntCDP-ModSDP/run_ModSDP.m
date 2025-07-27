clc;clear

load('.\example\ALL_ST.mat');
load('.\example\AML_TARGET.mat');
D{1}=ALL_ST;
D{2}=AML_TARGET;

num=length(D);
[~,txt]=xlsread('.\drivergene.xlsx');
gene_type=txt;   
[B,metagene,gene_type_new]=meta(D,gene_type);
mg=metagene;

gene_type_new0=gene_type_new;
f1=1;f2=1;f0=1;
for k=2:10
    k
    exclusion=[];    
    gene_type_new=gene_type_new0;
    max_value=1;
    %-----------Com-----------%
    for t=1:3
        max_geneset00=EntCDP_matlab(B,k,exclusion);
    
        if (max(max_geneset00(:,k+1))>max_value) || (max(max_geneset00(:,k+1))==max_value && length(find(max_geneset00(:,k+1)==max_value))>length(find(max_geneset(:,k+1)==max_value)))
            max_value=max(max_geneset00(:,k+1));
            max_geneset=max_geneset00;
        end
    end
    
    m0=find(max_geneset(:,k+1)==max_value);
    
    com_gene=cell(1,k);
    for i=1:k
        com_gene{i}=gene_type_new{max_geneset(m0(1),i)};
    end
    w=H_com(B,max_geneset(m0(1),1:k),exclusion);
    
    disp('Com')
    if sum(max_geneset(m0(1),k+2:k+2+num)<0.05)==num+1
        exclusion=max_geneset(m0(1),1:k);disp(com_gene);%disp(max_geneset(m0(1),k+1));disp(max_geneset(m0(1),k+2:k+num+2));disp(w)
        allgene0(1,f0:f0+k-1)=com_gene;f0=f0+k;
    else
        exclusion=[];disp('Not significant');
    end
    
    
    gene_type_new=gene_type_new0;
    %-----------Spe-----------%
    if exclusion>=1
        gene_type_new(exclusion)=[];
    end
    
    A{1}=B{1};
    for r=2:length(B)
        C{r-1}=B{r};
    end
    % A{1}=B{1};A{2}=B{2};%A{3}=B{3};
    % for r=3:length(B)
    %     C{r-2}=B{r};
    % end
    
    num1=length(A);num2=length(C);
    max_value=0;
    
    for t=1:3
        max_geneset0=ModSDP_GA_matlab(A,C,k,exclusion); 
        %max_geneset0=ModSDP_matlab(A,C,k,exclusion);
        %If IBM ILOG CPLEX optimizer is not added to your matlab, genetic algorithms can be used instead
        if (max(max_geneset0(:,k+1))>max_value) || (max(max_geneset0(:,k+1))==max_value && length(find(max_geneset0(:,k+1)==max_value))>length(find(max_geneset1(:,k+1)==max_value)))
            max_value=max(max_geneset0(:,k+1));
            max_geneset1=max_geneset0;
        end
    end
    
    
    m0=find(max_geneset1(:,k+1)==max_value);
    disp('spe1/2')

    spe_gene1=cell(1,k);
    for j=1:k
        spe_gene1{j}=gene_type_new{max_geneset1(m0(1),j)};
        for u=1:size(mg,1)
            s=0;
            for v=1:size(mg,2)
                if isempty(mg{u,v})==1
                    break;
                end
                s=s+1;
            end
            if length(find(categorical(mg(u,1:s))==spe_gene1{j}))==1
                disp('meta')
                disp(spe_gene1{j})
                break
            end
        end
    end
        
    if  sum(max_geneset1(m0(1),k+2:k+1+num1)<0.05)==num1 && sum(max_geneset1(m0(1),k+2+num1:k+1+num1+num2)>0.05)==num2 && sum(max_geneset1(m0(1),k+2+num1+num2)<0.05)==1
        disp(spe_gene1)
        %disp(max_geneset1(m0(1),1:k));
        %disp(max_geneset1(m0(1),k+1));
        disp(max_geneset1(m0(1),k+2:k+2+num1+num2));
        allgene1(1,f1:f1+k-1)=spe_gene1;f1=f1+k;
        
    else
        disp('Not significant');%disp(com_gene);%disp(w)
    end
    
    clear A C max_geneset0 max_geneset00 max_geneset1 spe_gene1;
    
    
    C{1}=B{1};
    for r=2:length(B)
        A{r-1}=B{r};
    end
    % C{1}=B{1};C{2}=B{2};%C{3}=B{3};
    % for r=3:length(B)
    %     A{r-2}=B{r};
    % end
    
    num1=length(A);num2=length(C);
    max_value=0;
    
    for t=1:3
        max_geneset0=ModSDP_GA_matlab(A,C,k,exclusion);
        if (max(max_geneset0(:,k+1))>max_value) || (max(max_geneset0(:,k+1))==max_value && length(find(max_geneset0(:,k+1)==max_value))>length(find(max_geneset2(:,k+1)==max_value)))
            max_value=max(max_geneset0(:,k+1));
            max_geneset2=max_geneset0;
        end
    end
    
    m0=find(max_geneset2(:,k+1)==max_value);
    disp('spe2/1')

    spe_gene2=cell(1,k);
    for j=1:k
        spe_gene2{j}=gene_type_new{max_geneset2(m0(1),j)};
        for u=1:size(mg,1)
            s=0;
            for v=1:size(mg,2)
                if isempty(mg{u,v})==1
                    break;
                end
                s=s+1;
            end
            if length(find(categorical(mg(u,1:s))==spe_gene2{j}))==1
                disp('meta')
                disp(spe_gene2{j})
                break
            end
        end
    end
    if  sum(max_geneset2(m0(1),k+2:k+1+num1)<0.05)==num1 && sum(max_geneset2(m0(1),k+2+num1:k+1+num1+num2)>0.05)==num2 && sum(max_geneset2(m0(1),k+2+num1+num2)<0.05)==1
        disp(spe_gene2)
        %disp(max_geneset2(m0(1),1:k));
        %disp(max_geneset2(m0(1),k+1));
        disp(max_geneset2(m0(1),k+2:k+2+num1+num2));
        allgene2(1,f2:f2+k-1)=spe_gene2;f2=f2+k;
    
    else
        disp('Not significant');%disp(com_gene);%disp(w)
    end
    clear max_geneset max_geneset1 max_geneset2 spe_gene1 spe_gene2 A C max_geneset00 max_geneset0

end



