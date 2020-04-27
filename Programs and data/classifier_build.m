load Subnetwork;
load GSE2034_ma2;
co=co(co>=0,1);

%%用网络A∩B的密集子网络建立子分类器
temp_network=subnetwork1;
m=max(temp_network(:,2));
for i=1:m    %i表示第i个密集子网络
    id=temp_network(temp_network(:,2)==i,1);  %找到第i个密集子网络的所有基因id
    if length(id)>=5
        [~,n]=ismember(id,ma2_geneId);  %n是一个列向量，表示id在ma2_geneId中的位置，即在ma2中的位置
        data=ma2(n,co>=0);
        N=length(data(1,:));
        attr_num=length(data(:,1));     
        index=(1:attr_num)';          
        Accuracy_sum=0;  %初始化变量，存储分类器的准确率
        auc_sum=0;       %初始化变量，存储分类器的auc
        sensitivity_sum=0; %初始化变量，存储分类器的灵敏度
        specificity_sum=0; %初始化变量，存储分类器的特异度
        mcc_sum=0;
        %10次5倍交叉验证
        for k=1:10
            Indices = crossvalind('Kfold', N, 5);%五倍交叉验证，将样本分成五份
            for l=1:5
                train=data(:,Indices~=l);     %训练集数据
                train_label=co(Indices~=l,1); %训练集标签
                test=data(:,Indices==l);      %测试集数据
                test_label=co(Indices==l,1);  %测试集标签
                [res,score]=centroid_classifier(train',train_label,test');
                temp1(i,Indices==l,k)=res;%保存子分类器对所有病人每一次的预测结果
                if var(score)~=0
                    [~,~,~,auc]=perfcurve(test_label,score(:,1),'0');
                else
                    auc=0.5;
                end
                auc_sum=auc_sum+auc;
                cpre=classperf(test_label,res);  %计算分类器的准确率等指标
                Accuracy_sum=Accuracy_sum+cpre.CorrectRate;  %准确率
                sensitivity_sum=sensitivity_sum+cpre.Sensitivity;  %灵敏度
                specificity_sum=specificity_sum+cpre.Specificity;  %特异度
                XX_train=cpre.DiagnosticTable;%混淆矩阵
                fenmu=(XX_train(1,1)+XX_train(1,2))*(XX_train(1,1)+XX_train(2,1))*(XX_train(2,2)+XX_train(1,2))*(XX_train(2,1)+XX_train(2,2));
                if (fenmu==0)
                    MCC_train=0;
                else
                    MCC_train=(XX_train(1,1)*XX_train(2,2)-XX_train(1,2)*XX_train(2,1))/sqrt(fenmu);
                end
                mcc_sum=mcc_sum+MCC_train;
            end
        end
        r1(i,1)=auc_sum/50;%r1保存子分类器的AUC,MCC,ACC等指标
        r1(i,2)=Accuracy_sum/50;
        r1(i,3)=sensitivity_sum/50;
        r1(i,4)=specificity_sum/50;
        r1(i,5)=mcc_sum/50;
    else
        r1(i,1)=0;
        r1(i,2)=0;
        r1(i,3)=0;
        r1(i,4)=0;
        r1(i,5)=0;
    end
end
%predict1保存各子分类器对每个病人的预测结果
for i=1:length(r1(:,1))
    for j=1:length(co)
        predict1(i,j)=sum(temp1(i,j,:));
        if(predict1(i,j)>5)
            predict1(i,j)=1;
        else
            predict1(i,j)=0;
        end
    end
end

%%用网络A∪B的密集子网络建立子分类器
temp_network=subnetwork2;
m=max(temp_network(:,2));
for i=1:m     %i表示第i个密集子网络
    id=temp_network(temp_network(:,2)==i,1);  %找到第i个密集子网络的所有基因id
    if length(id)>=5
        [~,n]=ismember(id,ma2_geneId);  %n是一个列向量，表示id在ma2_geneId中的位置，即在ma2中的位置
        data=ma2(n,co>=0);
        N=length(data(1,:));
        attr_num=length(data(:,1));     
        index=(1:attr_num)';            
        Accuracy_sum=0;  %初始化变量，存储分类器的准确率
        auc_sum=0;       %初始化变量，存储分类器的auc
        sensitivity_sum=0; %初始化变量，存储分类器的灵敏度
        specificity_sum=0; %初始化变量，存储分类器的特异度
        mcc_sum=0;
        %10次5倍交叉验证
        for k=1:10
            Indices = crossvalind('Kfold', N, 5);%五倍交叉验证，将样本分成五份
            for l=1:5
                train=data(:,Indices~=l);     %训练集数据
                train_label=co(Indices~= l,1); %训练集标签
                test=data(:,Indices==l);      %测试集数据
                test_label=co(Indices==l,1);  %测试集标签
                [res,score]=centroid_classifier(train',train_label,test');
                temp2(i,Indices==l,k)=res;%保存子分类器对所有病人每一次的预测结果
                if var(score)~=0
                    [~,~,~,auc]=perfcurve(test_label,score(:,1),'0');
                else
                    auc=0.5;
                end
                auc_sum=auc_sum+auc;
                cpre=classperf(test_label,res);  %计算分类器的准确率等指标
                Accuracy_sum=Accuracy_sum+cpre.CorrectRate;  %准确率
                sensitivity_sum=sensitivity_sum+cpre.Sensitivity;  %灵敏度
                specificity_sum=specificity_sum+cpre.Specificity;  %特异度
                XX_train=cpre.DiagnosticTable;%混淆矩阵
                fenmu=(XX_train(1,1)+XX_train(1,2))*(XX_train(1,1)+XX_train(2,1))*(XX_train(2,2)+XX_train(1,2))*(XX_train(2,1)+XX_train(2,2));
                if (fenmu==0)
                    MCC_train=0;
                else
                    MCC_train=(XX_train(1,1)*XX_train(2,2)-XX_train(1,2)*XX_train(2,1))/sqrt(fenmu);
                end
                mcc_sum=mcc_sum+MCC_train;
            end
        end
        r2(i,1)=auc_sum/50;%r2保存子分类器的AUC,MCC,ACC等指标
        r2(i,2)=Accuracy_sum/50;
        r2(i,3)=sensitivity_sum/50;
        r2(i,4)=specificity_sum/50;
        r2(i,5)=mcc_sum/50;
    else
        r2(i,1)=0;
        r2(i,2)=0;
        r2(i,3)=0;
        r2(i,4)=0;
        r2(i,5)=0;
    end
end
%predict2保存各子分类器对每个病人的预测结果
for i=1:length(r2(:,1))
    for j=1:length(co)
        predict2(i,j)=sum(temp2(i,j,:));
        if(predict2(i,j)>5)
            predict2(i,j)=1;
        else
            predict2(i,j)=0;
        end
    end
end

%用网络(A∩?B,?A∩B)的密集子网络建立子分类器
temp_network=subnetwork3;
m=max(temp_network(:,2));
for i=1:m     %i表示第i个密集子网络
    id=temp_network(temp_network(:,2)==i,1);  %找到第i个密集子网络的所有基因id
    if length(id)>=5
        [~,n]=ismember(id,ma2_geneId);  %n是一个列向量，表示id在ma2_geneId中的位置，即在ma2中的位置
        data=ma2(n,co>=0);
        N=length(data(1,:));
        attr_num=length(data(:,1));     
        index=(1:attr_num)';            
        Accuracy_sum=0;  %初始化变量，存储分类器的准确率
        auc_sum=0;       %初始化变量，存储分类器的auc
        sensitivity_sum=0; %初始化变量，存储分类器的灵敏度
        specificity_sum=0; %初始化变量，存储分类器的特异度
        mcc_sum=0;
        %10次5倍交叉验证
        for k=1:10
            Indices = crossvalind('Kfold', N, 5);%五倍交叉验证，将样本分成五份
            for l=1:5
                train=data(:,Indices~=l);     %训练集数据
                train_label=co(Indices~=l,1); %训练集标签
                test=data(:,Indices==l);      %测试集数据
                test_label=co(Indices==l,1);  %测试集标签
                [res,score]=centroid_classifier(train',train_label,test');
                temp3(i,Indices==l,k)=res;%保存子分类器对所有病人每一次的预测结果
                if var(score)~=0
                    [~,~,~,auc]=perfcurve(test_label,score(:,1),'0');
                else
                    auc=0.5;
                end
                auc_sum=auc_sum+auc;
                cpre=classperf(test_label,res);  %计算分类器的准确率等指标
                Accuracy_sum=Accuracy_sum+cpre.CorrectRate;  %准确率
                sensitivity_sum=sensitivity_sum+cpre.Sensitivity;  %灵敏度
                specificity_sum=specificity_sum+cpre.Specificity;  %特异度
                XX_train=cpre.DiagnosticTable;%混淆矩阵
                fenmu=(XX_train(1,1)+XX_train(1,2))*(XX_train(1,1)+XX_train(2,1))*(XX_train(2,2)+XX_train(1,2))*(XX_train(2,1)+XX_train(2,2));
                if (fenmu==0)
                    MCC_train=0;
                else
                    MCC_train=(XX_train(1,1)*XX_train(2,2)-XX_train(1,2)*XX_train(2,1))/sqrt(fenmu);
                end
                mcc_sum=mcc_sum+MCC_train;
            end
        end
        r3(i,1)=auc_sum/50;%r3保存子分类器的AUC,MCC,ACC等指标
        r3(i,2)=Accuracy_sum/50;
        r3(i,3)=sensitivity_sum/50;
        r3(i,4)=specificity_sum/50;
        r3(i,5)=mcc_sum/100;
    else
        r3(i,1)=0;
        r3(i,2)=0;
        r3(i,3)=0;
        r3(i,4)=0;
        r3(i,5)=0;
    end
end
%predict保存各子分类器对每个病人的预测结果
for i=1:length(r3(:,1))
    for j=1:length(co)
        predict3(i,j)=sum(temp3(i,j,:));
        if(predict3(i,j)>5)
            predict3(i,j)=1;
        else
            predict3(i,j)=0;
        end
    end
end

%%用多数投票策略集成分类器
for i=1:length(co)
    s1=sum(predict1(r1(:,1)>0.6,i));
    s2=sum(predict2(r2(:,1)>0.6,i));
    s3=sum(predict3(r3(:,1)>0.6,i));
    l1=length(predict1(r1(:,1)>0.6,i));
    l2=length(predict2(r2(:,1)>0.6,i)); 
    l3=length(predict3(r3(:,1)>0.6,i));
    predict_score1(i,2)=s1/l1;
    predict_score1(i,1)=1-predict_score1(i,2);
    predict_score2(i,2)=s2/l2;
    predict_score2(i,1)=1-predict_score2(i,2);
    predict_score3(i,2)=s3/l3;
    predict_score3(i,1)=1-predict_score3(i,2);

    predict_score(i,2)=(s1+s2+s3)/(l1+l2+l3);
    predict_score(i,1)=1-predict_score(i,2);
    if predict_score(i,1)>0.5
        predict(i,1)=0;
    else
        predict(i,1)=1;
    end
end

[~,~,~,s_auc1]=perfcurve(co,predict_score1(:,1),'0');
[~,~,~,s_auc2]=perfcurve(co,predict_score2(:,1),'0');
[~,~,~,s_auc3]=perfcurve(co,predict_score3(:,1),'0');
[~,~,~,s_auc]=perfcurve(co,predict_score(:,1),'0');
TP=sum(predict(co==0)==0);
FN=sum(predict(co==0)==1);
TN=sum(predict(co==1)==1);
FP=sum(predict(co==1)==0);
ACC=(TP+TN)/(TP+TN+FP+FN);
SN=TP/(TP+FP);
SP=TN/(TN+FN);
MCC=(TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));