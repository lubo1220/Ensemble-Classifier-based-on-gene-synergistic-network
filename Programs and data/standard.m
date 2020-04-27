load Test_example;%test_example.mat保存了5个数据集的基因表达和真实预后
for i=1:length(ma2(1,:))
    index1(:,i)=(1:length(ma2(:,1)))';
    temp1(:,1)=index1(:,i);
    temp1(:,2)=ma2(:,i);
    temp1=sortrows(temp1,2);
    index1(:,i)=temp1(:,1);
    data1(:,i)=temp1(:,2);
end

for i=1:length(GSE4922_ma2(1,:))
    index2(:,i)=(1:length(GSE4922_ma2(:,1)))';
    temp2(:,1)=index2(:,i);
    temp2(:,2)=GSE4922_ma2(:,i);
    temp2=sortrows(temp2,2);
    index2(:,i)=temp2(:,1);
    data2(:,i)=temp2(:,2); 
end

mid=mean(data1,2);
for i=1:length(data1(:,1))
    data1(i,:)=mid(i);
    data2(i,:)=mid(i);
end

for i=1:length(data1(1,:))
    temp1(:,1)=index1(:,i);
    temp1(:,2)=data1(:,i);
    temp1=sortrows(temp1);
    data1(:,i)=temp1(:,2);
end


for i=1:length(data2(1,:))
    temp2(:,1)=index2(:,i);
    temp2(:,2)=data2(:,i);
    temp2=sortrows(temp2);
    data2(:,i)=temp2(:,2);
end
