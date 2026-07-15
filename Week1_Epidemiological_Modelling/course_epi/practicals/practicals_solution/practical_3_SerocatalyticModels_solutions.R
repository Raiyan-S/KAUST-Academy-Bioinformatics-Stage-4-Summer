## Bring in data
dat <-read.csv("data/SeroStudy.csv")

## Look at the top 6 rows to see what the is in this data set
head(dat)

### What is the range of ages and mean age
range(dat[,3])
mean(dat$Age)
summary(dat$Age)

### What is mean seropositivity by sex
mean(dat$IgG[dat$Sex=="Male"])
mean(dat$IgG[dat$Sex=="Female"])

### What is mean seropositivity by age
ages      <- 1:15
PropByAge <- rep(NA,15)
for(i in 1:15)
{
	PropByAge[i]<-mean(dat$IgG[dat$Age==i])
}

plot(ages, PropByAge)

### For a given theta calculate the expected proportion seropositive by age
### assuming a constant force of infection
theta    <- 0.05
ages     <- 1:15
pSeropos <- 1 - exp(-theta * ages)
plot(ages, pSeropos, pch = 20)

### Calculate the likelihood for a range of theta values
thetas    <- seq(0.001,0.1,0.001)
llByTheta <- rep(NaN,length(thetas))

for(i in 1:length(thetas))
{
	ll <- (1 - dat$IgG) * (-thetas[i] * dat$Age) +
	  dat$IgG * log(1 - exp(-thetas[i] * dat$Age))
	llByTheta[i]<-sum(ll)
}
plot(thetas, llByTheta)

## Identify the maximum likelihood estimate
MLE <- thetas[which.max(llByTheta)]

plot(ages, PropByAge, pch = 20)

for (k in 1:length(thetas))
{
	pSeropos <- 1 - exp(-thetas[k] * ages)
	lines(ages,pSeropos,col="red")
}

pSeropos.MLE <- 1 - exp(-MLE*ages)
lines(ages, pSeropos.MLE, col="blue", lwd = 3)
