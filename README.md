# pharmcovigialanceprogramme-of-india
r <- read.csv("G:/Downloads/New folder/phvid/data oct 2612.csv")
library("dplyr", lib.loc="~/R/win-library/3.2")
library("PhViD", lib.loc="~/R/win-library/3.2")# install the dplyr and PhViD package
p <- as.PhViD(r)#create a dataset for PhViD package
c <- BCPNN(p,DECISION = 3,DECISION.THRES = 0,RANKSTAT = 2)# find out the IC value 
p1 <- PRR(DATABASE = p,DECISION = 3,RANKSTAT = 2,DECISION.THRES = 1)#find out the PRR value
n=p$N
#Create a data set for Chi square value
s  <- data.frame(c$ALLSIGNALS$`drug code`,c$ALLSIGNALS$`event effect`,c$ALLSIGNALS$count,c$ALLSIGNALS$`drug margin`,c$ALLSIGNALS$`event margin`)
#find out the Chi square value mutate value
m <- mutate(s,n2=n-s$c.ALLSIGNALS..drug.margin,n.2 = n-s$c.ALLSIGNALS..event.margin,n21=s$c.ALLSIGNALS..event.margin-s$c.ALLSIGNALS.count,n22=n2-n21,n12=s$c.ALLSIGNALS..drug.margin-s$c.ALLSIGNALS.count,d1=n*(s$c.ALLSIGNALS.count*n22-n12*n21)*(s$c.ALLSIGNALS.count*n22-n12*n21),d2=s$c.ALLSIGNALS..drug.margin*s$c.ALLSIGNALS..event.margin*n2*n.2,chi=d1/d2)
c1 <- data.frame(c$ALLSIGNALS,m)# create a data set on drug adr combination basis of the vale of IC.
p2 <- data.frame(p1$ALLSIGNALS,lbprr=exp(p1$ALLSIGNALS$`LB95(log(PRR))`))#create a data set and the find out Prr lower bond.
write.csv(c1,"G:/Downloads/New folder/phvid/c1.csv")
write.csv(p2,"G:/Downloads/New folder/phvid/p2.csv")
