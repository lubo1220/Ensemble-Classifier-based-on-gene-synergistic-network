addpath(genpath('mRMR_0.9_compiled2'));%添加计算熵函数所在文件夹路径
load GSE2034_ma2;
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

%获得单个基因背景向量
for i=1:10000
    temp1=co';
    temp2=ma(1,:);
    RandIndex=randperm(length(temp2));
    temp2=temp2(RandIndex);
    u(i)=(entropy(temp1)+entropy(temp2)-jointentropy(temp1,temp2))/entropy(temp1);
end

%获得A∩B背景向量
for i=1:10000
    temp=co';
    temp1=ma(1,:);
    temp2=ma(2,:);
    RandIndex1=randperm(length(temp1));
    RandIndex2=randperm(length(temp2));
    temp1=temp1(RandIndex1);
    temp2=temp2(RandIndex2);
    for k=1:patient_num
            if temp1(k)==1&temp2(k)==1
                temp3(k)=1;
            else
                temp3(k)=0;
            end
    end
    u1(i)=(entropy(temp)+entropy(temp3)-jointentropy(temp,temp3))/entropy(temp);
end

%获得A∪B背景向量
for i=1:10000
    temp=co';
    temp1=ma(1,:);
    temp2=ma(2,:);
    RandIndex1=randperm(length(temp1));
    RandIndex2=randperm(length(temp2));
    temp1=temp1(RandIndex1);
    temp2=temp2(RandIndex2);
    for k=1:patient_num
            if temp1(k)==1|temp2(k)==1
                temp3(k)=1;
            else
                temp3(k)=0;
            end
    end
    u2(i)=(entropy(temp)+entropy(temp3)-jointentropy(temp,temp3))/entropy(temp);
end
%获得A∩?B背景向量
for i=1:10000
    temp=co';
    temp1=ma(1,:);
    temp2=ma(2,:);
    RandIndex1=randperm(length(temp1));
    RandIndex2=randperm(length(temp2));
    temp1=temp1(RandIndex1);
    temp2=temp2(RandIndex2);
    for k=1:patient_num
            if temp1(k)==1&temp2(k)~=1
                temp3(k)=1;
            else
                temp3(k)=0;
            end
    end
    u3(i)=(entropy(temp)+entropy(temp3)-jointentropy(temp,temp3))/entropy(temp);
end
%获得?A∩B背景向量
for i=1:10000
    temp=co';
    temp1=ma(1,:);
    temp2=ma(2,:);
    RandIndex1=randperm(length(temp1));
    RandIndex2=randperm(length(temp2));
    temp1=temp1(RandIndex1);
    temp2=temp2(RandIndex2);
    for k=1:patient_num
            if temp1(k)~=1&temp2(k)==1
                temp3(k)=1;
            else
                temp3(k)=0;
            end
    end
    u4(i)=(entropy(temp)+entropy(temp3)-jointentropy(temp,temp3))/entropy(temp);
end