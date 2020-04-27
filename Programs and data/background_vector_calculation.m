addpath(genpath('mRMR_0.9_compiled2'));%��Ӽ����غ��������ļ���·��
load GSE2034_ma2;
%ѡ����Ч�Ĳ�������
data=ma2(:,co>=0);
co=co(co>=0); 
%�ֱ�������ɸѡ�����Ĳ��˵�����
gene_num=length(data(:,1));
patient_num=length(data(1,:));

medians=median(data,2);%���ma2ÿһ�е���ֵ
for i=1:gene_num
    ma(i,:)=imbinarize(data(i,:),medians(i,1));%ʹ���ݶ�ֵ��
end
ma=double(ma);

%��õ������򱳾�����
for i=1:10000
    temp1=co';
    temp2=ma(1,:);
    RandIndex=randperm(length(temp2));
    temp2=temp2(RandIndex);
    u(i)=(entropy(temp1)+entropy(temp2)-jointentropy(temp1,temp2))/entropy(temp1);
end

%���A��B��������
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

%���A��B��������
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
%���A��?B��������
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
%���?A��B��������
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