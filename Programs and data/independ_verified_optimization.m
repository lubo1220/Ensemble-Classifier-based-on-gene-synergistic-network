load Subnetwork;
load R;%R保存子分类器独立验证AUC,MCC,ACC等评价指标
co_1=co;
co_2=GSE4922_co;
N1=length(co_1(co_1>=0));
N2=length(co_2(co_2>=0));
%A∩B
temp_network=subnetwork1;
m=max(temp_network(:,2));
index=(1:m)';
index=index(r1(:,1)>0.6);
m=length(index);
k=0;
for i=1:m     %i表示第i个密集子网络
    id=temp_network(temp_network(:,2)==index(i),1);  %找到第i个密集子网络的所有基因id
    [~,n1]=ismember(id,ma2_geneId);  %n是一个列向量，表示id在ma2_geneId中的位置，即在ma2中的位置
    [~,n2]=ismember(id,ma2_geneId);
    for j=1:100
        Indices=crossvalind('Kfold', N1, 100);
        temp_train=data1(n1,co_1>=0);
        train=temp_train(:,Indices<=50);
        test=temp_train(:,Indices>50);
        temp_trainlabel=co_1(co_1>=0);
        train_label=temp_trainlabel(Indices<=50);
        test_label=temp_trainlabel(Indices>50);
        [res,score]=centroid_classifier(train',train_label,test');
        [~,~,~,auc]=perfcurve(test_label,score(:,1),'0');
        if auc>0.6
            test=data2(n2,co_2>=0);
            test_label=co_2(co_2>=0);
            [res,score]=centroid_classifier(train',train_label,test');
            k=k+1;
            tmp1(:,k)=res;
        end
    end
end
%%A∪B
temp_network=subnetwork2;
m=max(temp_network(:,2));
index=(1:m)';
index=index(r2(:,1)>0.6);
m=length(index);
k=0;
for i=1:m     %i表示第i个密集子网络
    id=temp_network(temp_network(:,2)==index(i),1);  %找到第i个密集子网络的所有基因id
    [~,n1]=ismember(id,ma2_geneId);  %n是一个竖向量，表示id在ma2_geneId中的位置，即在ma2中的位置
    [~,n2]=ismember(id,ma2_geneId);
    for j=1:100
        Indices = crossvalind('Kfold', N1, 100);
        temp_train=data1(n1,co_1>=0);
        train=temp_train(:,Indices<=50);
        test=temp_train(:,Indices>50);
        temp_trainlabel=co_1(co_1>=0);
        train_label=temp_trainlabel(Indices<=50);
        test_label=temp_trainlabel(Indices>50);
        [res,score]=centroid_classifier(train',train_label,test');
        [~,~,~,auc]=perfcurve(test_label,score(:,1),'0');
        if auc>0.6
            test=data2(n2,co_2>=0);
            test_label=co_2(co_2>=0);
            [res,score]=centroid_classifier(train',train_label,test');
            k=k+1;
            tmp2(:,k)=res;
        end
    end
end
%%(A∩?B,?A∩B)
temp_network=subnetwork3;
m=max(temp_network(:,2));
index=(1:m)';
index=index(r3(:,1)>0.6);
m=length(index);
k=0;
for i=1:m     %i表示第i个密集子网络
    id=temp_network(temp_network(:,2)==index(i),1);  %找到第i个密集子网络的所有基因id
    [~,n1]=ismember(id,ma2_geneId);  %n是一个竖向量，表示id在ma2_geneId中的位置，即在ma2中的位置
    [~,n2]=ismember(id,ma2_geneId);
    for j=1:100
        Indices = crossvalind('Kfold', N1, 100);
        temp_train=data1(n1,co_1>=0);
        train=temp_train(:,Indices<=50);
        test=temp_train(:,Indices>50);
        temp_trainlabel=co_1(co_1>=0);
        train_label=temp_trainlabel(Indices<=50);
        test_label=temp_trainlabel(Indices>50);
        [res,score]=centroid_classifier(train',train_label,test');
        [~,~,~,auc]=perfcurve(test_label,score(:,1),'0');
        if auc>0.6
            test=data2(n2,co_2>=0);
            test_label=co_2(co_2>=0);
            [res,score]=centroid_classifier(train',train_label,test');
            k=k+1;
            tmp3(:,k)=res;
        end
    end
end



for i=1:N2
    predict_score1(i,2)=sum(tmp1(i,:))/length(tmp1(i,:));
    predict_score1(i,1)=1-predict_score1(i,2);
    predict_score2(i,2)=sum(tmp2(i,:))/length(tmp2(i,:));
    predict_score2(i,1)=1-predict_score2(i,2);
    predict_score3(i,2)=sum(tmp3(i,:))/length(tmp3(i,:));
    predict_score3(i,1)=1-predict_score3(i,2);
    predict_score(i,2)=(sum(tmp1(i,:))+sum(tmp2(i,:))+sum(tmp3(i,:)))/(length(tmp1(i,:))+length(tmp2(i,:))+length(tmp3(i,:)));
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