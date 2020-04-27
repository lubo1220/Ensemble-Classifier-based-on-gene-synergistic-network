load Subnetwork;
load R;%R保存子分类器独立验证AUC,MCC,ACC等评价指标
co_1=co;
co_2=GSE4922_co;
%A∩B
temp_network=subnetwork1;
m=max(temp_network(:,2));
index=(1:m)';
index=index(r1(:,1)>0.60);
m=length(index);
temp1=zeros(length(co_2(co_2>=0)),m);
for i=1:m    %i表示第i个密集子网络
    id=temp_network(temp_network(:,2)==index(i),1);  %找到第i个密集子网络的所有基因id
    [~,n1]=ismember(id,ma2_geneId);  %n是一个竖向量，表示id在ma2_geneId中的位置，即在ma2中的位置
    [~,n2]=ismember(id,ma2_geneId);
    train=data1(n1,co_1>=0);
    train_label=co_1(co_1>=0);
    test=data2(n2,co_2>=0);
    test_label=co_2(co_2>=0);
    [res,score]=centroid_classifier(train',train_label,test');
    temp1(:,i)=res;
end
%%A∪B
temp_network=subnetwork2;
m=max(temp_network(:,2));
index=(1:m)';
index=index(r2(:,1)>0.60);
m=length(index);
temp2=zeros(length(co_2(co_2>=0)),m);
for i=1:m     %i表示第i个密集子网络
    id=temp_network(temp_network(:,2)==index(i),1);  %找到第i个密集子网络的所有基因id
    [~,n1]=ismember(id,ma2_geneId);  %n是一个竖向量，表示id在ma2_geneId中的位置，即在ma2中的位置
    [~,n2]=ismember(id,ma2_geneId);
    train=data1(n1,co_1>=0);
    train_label=co_1(co_1>=0);
    test=data2(n2,co_2>=0);
    test_label=co_2(co_2>=0);
    [res,score]=centroid_classifier(train',train_label,test');
    temp2(:,i)=res;
end
%%(A∩?B,?A∩B)
temp_network=subnetwork3;
m=max(temp_network(:,2));
index=(1:m)';
index=index(r3(:,1)>0.60);
m=length(index);
temp3=zeros(length(co_2(co_2>=0)),m);
for i=1:m     %i表示第i个密集子网络
    id=temp_network(temp_network(:,2)==index(i),1);  %找到第i个密集子网络的所有基因id
        [~,n1]=ismember(id,ma2_geneId);  %n是一个竖向量，表示id在ma2_geneId中的位置，即在ma2中的位置
        [~,n2]=ismember(id,ma2_geneId);
        train=data1(n1,co_1>=0);
        train_label=co_1(co_1>=0);
        test=data2(n2,co_2>=0);
        test_label=co_2(co_2>=0);
        [res,score]=centroid_classifier(train',train_label,test');
        temp3(:,i)=res;
end



for i=1:length(co_2(co_2>=0))
    predict_score1(i,2)=sum(temp1(i,:))/length(temp1(i,:));
    predict_score1(i,1)=1-predict_score1(i,2);
    predict_score2(i,2)=sum(temp2(i,:))/length(temp2(i,:));
    predict_score2(i,1)=1-predict_score2(i,2);
    predict_score3(i,2)=sum(temp3(i,:))/length(temp3(i,:));
    predict_score3(i,1)=1-predict_score3(i,2);

    predict_score(i,2)=(sum(temp1(i,:))+sum(temp2(i,:))+sum(temp3(i,:)))/(length(temp1(i,:))+length(temp2(i,:))+length(temp3(i,:)));
    predict_score(i,1)=1-predict_score(i,2);
    if predict_score(i,1)>0.5
         predict(i,1)=0;
    else
         predict(i,1)=1;
    end
end
[~,~,~,s_auc1]=perfcurve(co_2(co_2>=0),predict_score1(:,1),'0'); 
[~,~,~,s_auc2]=perfcurve(co_2(co_2>=0),predict_score2(:,1),'0');     
[~,~,~,s_auc3]=perfcurve(co_2(co_2>=0),predict_score3(:,1),'0');     
[~,~,~,s_auc]=perfcurve(co_2(co_2>=0),predict_score(:,1),'0');      
TP=sum(predict(co_2(co_2>=0)==0)==0);
FN=sum(predict(co_2(co_2>=0)==0)==1);
TN=sum(predict(co_2(co_2>=0)==1)==1);
FP=sum(predict(co_2(co_2>=0)==1)==0);
ACC=(TP+TN)/(TP+TN+FP+FN);
SN=TP/(TP+FP);
SP=TN/(TN+FN);
MCC=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));