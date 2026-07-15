library(tidyverse)

dat <- read.csv("~/Course_Materials/Week1_Epidemiological_Modelling/course_epi/data/SeroStudy.csv")
head(dat)

# 1
summary(dat)
# Mean = 8.0 | Min = 1.0 | Max 15.0

# 2 & 3
dat |> 
  group_by(Sex) |> 
  summarise(mean(IgG, na.rm = TRUE))
# Female = 0.389  Male = 0.375

# 4
dat |> 
  group_by(Age) |> 
  summarise(mean_IgG = mean(IgG, na.rm = TRUE)) |> 
  ggplot(aes(x = Age, y = mean_IgG)) +
  geom_point()

mod <- dat |> 
  group_by(Age) |> 
  summarise(mean_IgG = mean(IgG, na.rm = TRUE))
# Age 14 is the HigAge# Age 14 is the Highest, Age 1 is the Lowest

# 5
theta <- 0.05
ages  <- 1:15
equa <- 1-(exp(-theta*ages))
model_data <- data.frame(age = ages, predicted_val = equa)

ggplot(data = model_data, aes(age,predicted_val)) +
  geom_point() 

# 6
mod$probsero <- 1-(exp(-theta*ages))
ggplot(data = mod, aes(Age,mean_IgG)) +
  geom_point() +
  geom_line(aes(y=probsero))

# 7
theta <- 0.07
ages  <- 1:15
mod$probsero <- 1-(exp(-theta*ages))
ggplot(data = mod, aes(Age,mean_IgG)) +
  geom_point() +
  geom_line(aes(y=probsero))

# 8 
thetas  <- seq(0.001,0.1,0.001)
llByTheta <- rep(NaN,length(thetas))
for(i in 1:length(thetas)){
  ll <- dat |> 
    with((1 - as.integer(IgG)) * (-thetas[i] * Age) + 
           as.integer(IgG) * log(1 - exp(-thetas[i] * Age)))
  llByTheta[i] <- sum(ll, na.rm = TRUE)
}
plot(thetas, llByTheta)
llByTheta

# 9 
max(llByTheta)
thetas[which.max(llByTheta)]

theta <- thetas[which.max(llByTheta)]
ages  <- 1:15
mod$probsero <- 1-(exp(-theta*ages))
ggplot(data = mod, aes(Age,mean_IgG)) +
  geom_point() +
  geom_line(aes(y=probsero))

# 10 
print(1-exp(-theta*30))

# 11
ages  <- 1:60
plot(ages, 1 - exp(-theta * ages))

# 12
ages  <- 1:60
prop_seropositive <- exp(-theta * ages)
seropositive_counts <- theta * prop_seropositive
total_seropositive <- seropositive_counts * 10000

# Print result
print(sum(total_seropositive))
     