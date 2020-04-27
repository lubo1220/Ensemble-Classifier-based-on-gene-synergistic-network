function [cpre,score]=centroid_classifier(train_data,train_label,test_data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Use centroid clssifier to predict the labels of the  test_data, based on the train_data and train_label.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
train_data=train_data';

c_p=mean(train_data(:,train_label==0),2);

c_n=mean(train_data(:,train_label==1),2);
c_m=(c_p+c_n)/2;
w=c_p-c_n;
test_data=test_data';
n=length(test_data(1,:));
for i=1:n
    score(i,1)=sum((test_data(:,i)-c_m).*w);
    if score(i,1)>0
        cpre(i,1)=0;
    else
        cpre(i,1)=1;
    end;
end;
end
