addpath(genpath('mRMR_0.9_compiled2'));%添加计算熵函数所在文件夹路径
load GSE2034_ma2;
load Background_vector;
%选择有效的病人数据
data=ma2(:,co>=0);
co=co(co>=0); 
%分别求基因和筛选出来的病人的总数
gene_num=length(data(:,1));
patient_num=length(data(1,:));

medians=median(data,2);%求出ma2每一行的中值
for i=1:gene_num
    ma(i,:)=imbinarize(data(i,:),medians(i,1));%使数据二值化
end
ma=double(ma);

%n用来保存3个子网络有多少基因对
n1=0;
n2=0;
n3=0;
%对所有可能的基因对进行计算
for i=1:gene_num-1
    for j=i+1:gene_num 
        A=ma(i,:);
        B=ma(j,:);
        for k=1:patient_num
            if A(k)==1&B(k)==1
                C1(k)=1;  
            else
                C1(k)=0;
            end
            if A(k)==1|B(k)==1
                C2(k)=1;
            else
                C2(k)=0;
            end
            if A(k)==1&B(k)~=1
                C3(k)=1;
            else
                C3(k)=0;
            end
            if A(k)~=1&B(k)==1
                C4(k)=1;
            else
                C4(k)=0;
            end
        end
        U_a=(entropy(co')+entropy(A)-jointentropy(co',A))/entropy(co');
        p_a=sum(U_a<=u)/10000;
        U_b=(entropy(co')+entropy(B)-jointentropy(co',B))/entropy(co');
        p_b=sum(U_b<=u)/10000;
        U_c1=(entropy(co')+entropy(C1)-jointentropy(co',C1))/entropy(co');
        U_c2=(entropy(co')+entropy(C2)-jointentropy(co',C2))/entropy(co');
        U_c3=(entropy(co')+entropy(C3)-jointentropy(co',C3))/entropy(co');
        U_c4=(entropy(co')+entropy(C4)-jointentropy(co',C4))/entropy(co');
        p_c1=sum(U_c1<=u1)/10000;
        p_c2=sum(U_c2<=u2)/10000;
        p_c3=sum(U_c3<=u3)/10000;
        p_c4=sum(U_c4<=u4)/10000;
        %relation保存三种网络的边
        if p_a>0.001&p_b>0.001
            if p_c1<=0.001
                n1=n1+1;
                relation1(n1,1)=ma2_geneId(i);
                relation1(n1,2)=ma2_geneId(j);
                relation1(n1,3)=p_c1;
            end
             if p_c2<=0.001
                n2=n2+1;
                relation2(n2,1)=ma2_geneId(i);
                relation2(n2,2)=ma2_geneId(j);
                relation2(n2,3)=p_c2;
             end
             if p_c3<=0.001|p_c4<=0.001
                n3=n3+1;
                relation3(n3,1)=ma2_geneId(i);
                relation3(n3,2)=ma2_geneId(j);
                relation3(n3,3)=p_c3;
             end
        end
    end
end


