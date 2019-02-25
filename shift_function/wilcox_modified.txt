## ===================================================
## The functions in this file were modified from the
## original code from Rand Wilcox:
## http://dornsife.usc.edu/labs/rwilcox/software/
## ===================================================

shifthd<-function(x,y,nboot=200){
#
#   Compute confidence intervals for the difference between deciles
#   of two independent groups. The simultaneous probability coverage is .95.
#   The Harrell-Davis estimate of the qth quantile is used.
#   The default number of bootstrap samples is nboot=200
#
#   The results are stored and returned in a 9 by 5 matrix,
#   the ith row corresponding to the i/10 quantile.
#   The first column is the lower end of the confidence interval.
#   The second column is the upper end.
#   The third column is the estimated difference between the deciles
#   (second group minus first).
#
# Modified from Rallfun-v31: GAR - University of Glasgow - 2016-06-22
# simplified inputs/outputs to make independent figures
# changed difference to x-y instead of y-x for consistency across functions
# now outputs confidence intervals, similarly to shiftdhd
# changed labels to reflect nature of data: 'groups' here, 'conditions' for shiftdhd
# now returns a data frame for use with ggplot2
x<-x[!is.na(x)]
y<-y[!is.na(y)]
crit<-80.1/(min(length(x),length(y)))^2+2.73
m<-matrix(0,9,5)
for (d in 1:9){
q<-d/10
data<-matrix(sample(x,size=length(x)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,hd,q)
sex<-var(bvec)
data<-matrix(sample(y,size=length(y)*nboot,replace=TRUE),nrow=nboot)
bvec<-apply(data,1,hd,q)
sey<-var(bvec)
m[d,1]<-hd(x,q)
m[d,2]<-hd(y,q)
m[d,3]<-m[d,1]-m[d,2]
m[d,4]<-m[d,3]-crit*sqrt(sex+sey)
m[d,5]<-m[d,3]+crit*sqrt(sex+sey)
}
out<-data.frame(m)
names(out)<-c('group1','group2','difference','ci_lower','ci_upper')
out
}

shiftdhd<-function(x,y,nboot=200){
#
#   Compute confidence intervals for the difference between deciles
#   of two dependent groups. The simultaneous probability coverage is .95.
#   The Harrell-Davis estimate of the qth quantile is used.
#   The default number of bootstrap samples is nboot=200
#
#   The results are stored and returned in a 9 by 5 matrix,
#   the ith row corresponding to the i/10 quantile.
#   The first column is the lower end of the confidence interval.
#   The second column is the upper end.
#   The third column is the estimated difference between the deciles
#   (second group minus first).
#   The fourth column contains the estimated standard error.
#
#   No missing values are allowed.
#
#   If the goal is to use an alpha value different from .05,
#   use the function qcomdhd or  qdec2ci
#
# Modified from Rallfun-v31: GAR - University of Glasgow - 2016-06-22
# simplified inputs/outputs for use with new shift_function_figure()
# changed labels to reflect nature of data: 'conditions' here, 'groups' for shifthd
# now returns a data frame for use with ggplot2
xy=elimna(cbind(x,y))
x=xy[,1]
y=xy[,2]
crit<-37/length(x)^(1.4)+2.75
m<-matrix(0,9,5)
data<-matrix(sample(length(x),size=length(x)*nboot,replace=TRUE),nrow=nboot)
xmat<-matrix(x[data],nrow=nboot,ncol=length(x))
ymat<-matrix(y[data],nrow=nboot,ncol=length(x))
for (d in 1:9){
q<-d/10
bvec<-apply(xmat,1,hd,q)-apply(ymat,1,hd,q)
se<-sqrt(var(bvec))
m[d,1]=hd(x,q)
m[d,2]=hd(y,q)
m[d,3]<-m[d,1]-m[d,2]
m[d,4]<-m[d,3]-crit*se
m[d,5]<-m[d,3]+crit*se
}
out<-data.frame(m)
names(out)<-c('condition1','condition2','difference','ci_lower','ci_upper')
out
}
