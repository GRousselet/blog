#
#License: USC-RL v1.0
#The Software is made available for academic or non-commercial purposes only. The license is for
#a copy of the program for an unlimited term. Individuals requesting a license for commercial use must pay for a commercial license.
# USC Stevens Institute for Innovation University of Southern California
#1150 S. Olive Street, Suite 2300
#Los Angeles, CA 90115, USA
#ATTN: Accounting
#DISCLAIMER.  USC MAKES NO EXPRESS OR IMPLIED WARRANTIES, EITHER IN FACT OR BY
#OPERATION OF LAW, BY STATUTE OR OTHERWISE, AND USC SPECIFICALLY AND EXPRESSLY
#DISCLAIMS ANY EXPRESS OR IMPLIED WARRANTY OF MERCHANTABILITY OR FITNESS FOR A
#PARTICULAR PURPOSE, VALIDITY OF THE SOFTWARE OR ANY OTHER INTELLECTUAL PROPERTY
#RIGHTS OR NON-INFRINGEMENT OF THE INTELLECTUAL PROPERTY OR OTHER RIGHTS OF ANY
#THIRD PARTY. SOFTWARE IS MADE AVAILABLE AS-IS.
#LIMITATION OF LIABILITY.  TO THE MAXIMUM EXTENT PERMITTED BY LAW, IN NO EVENT WILL
#USC BE LIABLE TO ANY USER OF THIS CODE FOR ANY INCIDENTAL, CONSEQUENTIAL, EXEMPLARY
#OR PUNITIVE DAMAGES OF ANY KIND, LOST GOODWILL, LOST PROFITS, LOST BUSINESS AND/OR
#ANY INDIRECT ECONOMIC DAMAGES WHATSOEVER, REGARDLESS OF WHETHER SUCH DAMAGES
#ARISE FROM CLAIMS BASED UPON CONTRACT, NEGLIGENCE, TORT (INCLUDING STRICT LIABILITY
#OR OTHER LEGAL THEORY), A BREACH OF ANY WARRANTY OR TERM OF THIS AGREEMENT, AND
#REGARDLESS OF WHETHER USC WAS ADVISED OR HAD REASON TO KNOW OF THE POSSIBILITY OF
#INCURRING SUCH DAMAGES IN ADVANCE.
#For commercial license pricing and annual commercial update and support pricing, please
#contact:
#<Licensing Associate Name>
#USC Stevens Institute for Innovation
#University of Southern California
#1150 S. Olive Street, Suite 2300
#Los Angeles, CA 90015, USA
#Tel: <Licensing Associate phone number>
#Fax: +1 213-821-5001
#Email: a
#<Licensing Associate Email>
#and cc to:
#accounting@stevens.usc.edu


#  Last update:
#  July., 2022

comvar2 <- function(x,y,nboot=1000){
  #
  #  Compare the variances of two independent groups.
  #
  x<-x[!is.na(x)]  # Remove missing values in x
  y<-y[!is.na(y)]  # Remove missing values in y
  # set seed of random number generator so that
  # results can be duplicated.
  est1=var(x)
  est2=var(y)
  sig<-est1-est2
  # if(SEED)set.seed(2)
  nmin<-min(length(x),length(y))
  datax<-matrix(sample(x,size=nmin*nboot,replace=TRUE),nrow=nboot)
  datay<-matrix(sample(y,size=nmin*nboot,replace=TRUE),nrow=nboot)
  v1<-apply(datax,1,FUN=var)
  v2<-apply(datay,1,FUN=var)
  boot<-v1-v2
  boot<-sort(boot)
  ilow <- 15
  ihi <- 584
  if(nmin < 250) {
    ilow <- 13
    ihi <- 586
  }
  if(nmin < 180) {
    ilow <- 10
    ihi <- 589
  }
  if(nmin < 80) {
    ilow <- 7
    ihi <- 592
  }
  if(nmin < 40) {
    ilow <- 6
    ihi <- 593
  }
  ilow<-round((ilow/599)*nboot)
  ihi<-round((ihi/599)*nboot)
  ci<-c(boot[ilow+1],boot[ihi])
  list(n=c(length(x),length(y)),ci=ci,est.1=est1,est.2=est2,vardif=sig,ratio=est1/est2)
}

comvar2_wrap <- function(x,y,n,nboot=1000){
  #
  #  Compare the variances of two independent groups.
  #  Calls comvar2_cpp to speed up bootstrap.
  #  Only compute ci
  #  Assume equal sample sizes
  bootdiff <- comvar2_cpp(x,y,n,nboot)
  ilow <- 15
  ihi <- 584
  nmin <- n
  if(nmin < 250) {
    ilow <- 13
    ihi <- 586
  }
  if(nmin < 180) {
    ilow <- 10
    ihi <- 589
  }
  if(nmin < 80) {
    ilow <- 7
    ihi <- 592
  }
  if(nmin < 40) {
    ilow <- 6
    ihi <- 593
  }
  ilow <- round((ilow/599)*nboot)
  ihi <- round((ihi/599)*nboot)
  ci <- c(bootdiff[ilow+1],bootdiff[ihi])
  ci
  # list(n=c(length(x),length(y)),ci=ci,est.1=est1,est.2=est2,vardif=sig,ratio=est1/est2)
}

# library(rbenchmark)
# benchmark("R" = {comvar2(x,y)}, "C++" = {comvar2_wrap(x,y,n,nboot)}, 
#           replications = 500,
#           columns = c("test", "replications", "elapsed", "relative"))
