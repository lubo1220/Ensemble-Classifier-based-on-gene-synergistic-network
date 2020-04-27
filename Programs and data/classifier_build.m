load Subnetwork;
load GSE2034_ma2;
co=co(co>=0,1);

%%������A��B���ܼ������罨���ӷ�����
temp_network=subnetwork1;
m=max(temp_network(:,2));
for i=1:m    %i��ʾ��i���ܼ�������
    id=temp_network(temp_network(:,2)==i,1);  %�ҵ���i���ܼ�����������л���id
    if length(id)>=5
        [~,n]=ismember(id,ma2_geneId);  %n��һ������������ʾid��ma2_geneId�е�λ�ã�����ma2�е�λ��
        data=ma2(n,co>=0);
        N=length(data(1,:));
        attr_num=length(data(:,1));     
        index=(1:attr_num)';          
        Accuracy_sum=0;  %��ʼ���������洢��������׼ȷ��
        auc_sum=0;       %��ʼ���������洢��������auc
        sensitivity_sum=0; %��ʼ���������洢��������������
        specificity_sum=0; %��ʼ���������洢�������������
        mcc_sum=0;
        %10��5��������֤
        for k=1:10
            Indices = crossvalind('Kfold', N, 5);%�屶������֤���������ֳ����
            for l=1:5
                train=data(:,Indices~=l);     %ѵ��������
                train_label=co(Indices~=l,1); %ѵ������ǩ
                test=data(:,Indices==l);      %���Լ�����
                test_label=co(Indices==l,1);  %���Լ���ǩ
                [res,score]=centroid_classifier(train',train_label,test');
                temp1(i,Indices==l,k)=res;%�����ӷ����������в���ÿһ�ε�Ԥ����
                if var(score)~=0
                    [~,~,~,auc]=perfcurve(test_label,score(:,1),'0');
                else
                    auc=0.5;
                end
                auc_sum=auc_sum+auc;
                cpre=classperf(test_label,res);  %�����������׼ȷ�ʵ�ָ��
                Accuracy_sum=Accuracy_sum+cpre.CorrectRate;  %׼ȷ��
                sensitivity_sum=sensitivity_sum+cpre.Sensitivity;  %������
                specificity_sum=specificity_sum+cpre.Specificity;  %�����
                XX_train=cpre.DiagnosticTable;%��������
                fenmu=(XX_train(1,1)+XX_train(1,2))*(XX_train(1,1)+XX_train(2,1))*(XX_train(2,2)+XX_train(1,2))*(XX_train(2,1)+XX_train(2,2));
                if (fenmu==0)
                    MCC_train=0;
                else
                    MCC_train=(XX_train(1,1)*XX_train(2,2)-XX_train(1,2)*XX_train(2,1))/sqrt(fenmu);
                end
                mcc_sum=mcc_sum+MCC_train;
            end
        end
        r1(i,1)=auc_sum/50;%r1�����ӷ�������AUC,MCC,ACC��ָ��
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
%predict1������ӷ�������ÿ�����˵�Ԥ����
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

%%������A��B���ܼ������罨���ӷ�����
temp_network=subnetwork2;
m=max(temp_network(:,2));
for i=1:m     %i��ʾ��i���ܼ�������
    id=temp_network(temp_network(:,2)==i,1);  %�ҵ���i���ܼ�����������л���id
    if length(id)>=5
        [~,n]=ismember(id,ma2_geneId);  %n��һ������������ʾid��ma2_geneId�е�λ�ã�����ma2�е�λ��
        data=ma2(n,co>=0);
        N=length(data(1,:));
        attr_num=length(data(:,1));     
        index=(1:attr_num)';            
        Accuracy_sum=0;  %��ʼ���������洢��������׼ȷ��
        auc_sum=0;       %��ʼ���������洢��������auc
        sensitivity_sum=0; %��ʼ���������洢��������������
        specificity_sum=0; %��ʼ���������洢�������������
        mcc_sum=0;
        %10��5��������֤
        for k=1:10
            Indices = crossvalind('Kfold', N, 5);%�屶������֤���������ֳ����
            for l=1:5
                train=data(:,Indices~=l);     %ѵ��������
                train_label=co(Indices~= l,1); %ѵ������ǩ
                test=data(:,Indices==l);      %���Լ�����
                test_label=co(Indices==l,1);  %���Լ���ǩ
                [res,score]=centroid_classifier(train',train_label,test');
                temp2(i,Indices==l,k)=res;%�����ӷ����������в���ÿһ�ε�Ԥ����
                if var(score)~=0
                    [~,~,~,auc]=perfcurve(test_label,score(:,1),'0');
                else
                    auc=0.5;
                end
                auc_sum=auc_sum+auc;
                cpre=classperf(test_label,res);  %�����������׼ȷ�ʵ�ָ��
                Accuracy_sum=Accuracy_sum+cpre.CorrectRate;  %׼ȷ��
                sensitivity_sum=sensitivity_sum+cpre.Sensitivity;  %������
                specificity_sum=specificity_sum+cpre.Specificity;  %�����
                XX_train=cpre.DiagnosticTable;%��������
                fenmu=(XX_train(1,1)+XX_train(1,2))*(XX_train(1,1)+XX_train(2,1))*(XX_train(2,2)+XX_train(1,2))*(XX_train(2,1)+XX_train(2,2));
                if (fenmu==0)
                    MCC_train=0;
                else
                    MCC_train=(XX_train(1,1)*XX_train(2,2)-XX_train(1,2)*XX_train(2,1))/sqrt(fenmu);
                end
                mcc_sum=mcc_sum+MCC_train;
            end
        end
        r2(i,1)=auc_sum/50;%r2�����ӷ�������AUC,MCC,ACC��ָ��
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
%predict2������ӷ�������ÿ�����˵�Ԥ����
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

%������(A��?B,?A��B)���ܼ������罨���ӷ�����
temp_network=subnetwork3;
m=max(temp_network(:,2));
for i=1:m     %i��ʾ��i���ܼ�������
    id=temp_network(temp_network(:,2)==i,1);  %�ҵ���i���ܼ�����������л���id
    if length(id)>=5
        [~,n]=ismember(id,ma2_geneId);  %n��һ������������ʾid��ma2_geneId�е�λ�ã�����ma2�е�λ��
        data=ma2(n,co>=0);
        N=length(data(1,:));
        attr_num=length(data(:,1));     
        index=(1:attr_num)';            
        Accuracy_sum=0;  %��ʼ���������洢��������׼ȷ��
        auc_sum=0;       %��ʼ���������洢��������auc
        sensitivity_sum=0; %��ʼ���������洢��������������
        specificity_sum=0; %��ʼ���������洢�������������
        mcc_sum=0;
        %10��5��������֤
        for k=1:10
            Indices = crossvalind('Kfold', N, 5);%�屶������֤���������ֳ����
            for l=1:5
                train=data(:,Indices~=l);     %ѵ��������
                train_label=co(Indices~=l,1); %ѵ������ǩ
                test=data(:,Indices==l);      %���Լ�����
                test_label=co(Indices==l,1);  %���Լ���ǩ
                [res,score]=centroid_classifier(train',train_label,test');
                temp3(i,Indices==l,k)=res;%�����ӷ����������в���ÿһ�ε�Ԥ����
                if var(score)~=0
                    [~,~,~,auc]=perfcurve(test_label,score(:,1),'0');
                else
                    auc=0.5;
                end
                auc_sum=auc_sum+auc;
                cpre=classperf(test_label,res);  %�����������׼ȷ�ʵ�ָ��
                Accuracy_sum=Accuracy_sum+cpre.CorrectRate;  %׼ȷ��
                sensitivity_sum=sensitivity_sum+cpre.Sensitivity;  %������
                specificity_sum=specificity_sum+cpre.Specificity;  %�����
                XX_train=cpre.DiagnosticTable;%��������
                fenmu=(XX_train(1,1)+XX_train(1,2))*(XX_train(1,1)+XX_train(2,1))*(XX_train(2,2)+XX_train(1,2))*(XX_train(2,1)+XX_train(2,2));
                if (fenmu==0)
                    MCC_train=0;
                else
                    MCC_train=(XX_train(1,1)*XX_train(2,2)-XX_train(1,2)*XX_train(2,1))/sqrt(fenmu);
                end
                mcc_sum=mcc_sum+MCC_train;
            end
        end
        r3(i,1)=auc_sum/50;%r3�����ӷ�������AUC,MCC,ACC��ָ��
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
%predict������ӷ�������ÿ�����˵�Ԥ����
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

%%�ö���ͶƱ���Լ��ɷ�����
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