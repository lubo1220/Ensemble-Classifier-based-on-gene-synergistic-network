load Subnetwork;
load R; 
load GSE2034_ma2; 

temp_network=subnetwork1;
m=max(temp_network(:,2));
index=(1:m)';
index=index(r1(:,1)>0.6);
m=length(index);
Id1=zeros(1000,m);
for i=1:m     %i表示第i个密集子网络
    id=temp_network(temp_network(:,2)==index(i),1);  %找到第i个密集子网络的所有基因id
    Id1(1:length(id),i)=id;
end

temp_network=subnetwork2;
m=max(temp_network(:,2));
index=(1:m)';
index=index(r2(:,1)>0.6);
m=length(index);
Id2=zeros(1000,m);
for i=1:m     %i表示第i个密集子网络
    id=temp_network(temp_network(:,2)==index(i),1);  %找到第i个密集子网络的所有基因id
    Id2(1:length(id),i)=id;
end

temp_network=subnetwork3_new;
m=max(temp_network(:,2));
index=(1:m)';
index=index(r3(:,1)>0.6);
m=length(index);
Id3=zeros(1000,m);
for i=1:m     %i表示第i个密集子网络
    id=temp_network(temp_network(:,2)==index(i),1);  %找到第i个密集子网络的所有基因id
    Id3(1:length(id),i)=id;
end