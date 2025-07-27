function H=shang(x)

p=x/sum(x);
if sum(find(p==0))==0
    H=-sum(p.*log2(p));
else
    p1=p;
    p1(find(p1==0))=[];
    H=-sum(p1.*log2(p1));
end