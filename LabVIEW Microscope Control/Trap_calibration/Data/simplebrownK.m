Data= lvm_import();
Data=Data.Segment1.data;
LR=Data(:,1);
BT=Data(:,2);
SUM=Data(:,3);
nLR=LR./SUM;
nBT=BT./SUM;
znLR=nLR-mean(nLR);
znBT=nBT-mean(nBT);
meanznLR=movmean(znLR,100);
meanznBT=movmean(znBT,100);
plot(meanznLR,meanznBT)

fitLR = fitdist(meanznLR,'normal');
fitBT = fitdist(meanznBT,'normal');
