rm(list = ls())
set.seed(51)

## - Input the data from Figure 18.2 in Gelman, Hill Vehtari
name = c("Audrey", "Anna", "Bob", "Bill", "Caitlin", "Cara", "Dave", "Doug")
female = c(1,1,0,0,1,1,0,0)
age = c(40,40,50,50,60,60,70,70)
treatment = c(0,0,0,0,1,1,1,1)
y0 = c(140,140,150,150,160,160,170,170)
y1 = c(135,135,140,140,155,155,160,160)
y = y0
y[treatment==1] = y1[treatment==1]

print(cbind(name, female, age, treatment, y0, y1, y))

n = length(name)
tau_sate.true = mean(y1) - mean(y0)
tau_sate.true

tau_sate.est = mean(y[treatment==1]) - mean(y[treatment==0])
tau_sate.est

## - Illustration of randomization-based inference under various
#    assignment mechanisms

# ---- Bernoulli Trial
q = 0.5 #probability of each individual being treated
Z= rbinom(n,1,q)
y = y0  # Recalculate values for y based off new assignment vector Z
y[Z==1] = y1[Z==1]
print(cbind(name, female, age, Z, y0, y1, y))
mean(y[Z==1]) - mean(y[Z==0])

# Repeat the Bernoulli trial to get the randomization distribution
#  note: could enumerate every possible Z.bern, but we will simulate
#        more 'reps' than necessary as a general procedure
n_reps = 10000
Z.bern = matrix(NA, n_reps,n) # a matrix to store all of the different Z vectors
for (i in 1:n_reps)
  Z.bern[i,] = rbinom(n,1,q)

# Look at the assignment for a single randomization
print(cbind(name, female, age, Z.bern[551,], y0, y1))
mean(y1[Z.bern[551,]==1]) - mean(y0[Z.bern[551,]==0])

# Look at the distribution of n_t across all randomizations
nt = rowSums(Z.bern) 
hist(nt)

# Calculate the mean difference between treated/control for each 
# randomization
tau_sate.bern = rep(NA, n_reps)
for (i in 1:n_reps)
  tau_sate.bern[i] = mean(y1[Z.bern[i,]==1]) - mean(y0[Z.bern[i,]==0])

# look what happens when the number of treated is either 0 or 8
tau_sate.bern[which(nt==0)]
tau_sate.bern[which(nt==n)]
table(is.na(tau_sate.bern))

m.bern = round(mean(tau_sate.bern, na.rm = TRUE), 2)
sd.bern = round(sd(tau_sate.bern, na.rm = TRUE), 2)
hist(tau_sate.bern, main = paste("Bernoulli: Mean = ", m.bern, " SD = ", sd.bern, sep=""))



# ---- Completely Randomized Design
Z = sample(c(rep(0,4), rep(1,4)), 8)
y = y0  # Recalculate values for y based off new assignment vector Z
y[Z==1] = y1[Z==1]
print(cbind(name, female, age, Z, y0, y1,y))
mean(y[Z==1]) - mean(y[Z==0])

# Repeat the completely randomized trial to get the randomization distribution
#  note: could enumerate every possible Z.bern, but we will simulate
#        more 'reps' than necessary as a general procedure
n_reps = 10000
Z.cr = matrix(NA, n_reps,n) # a matrix to store all of the different Z vectors
for (i in 1:n_reps)
  Z.cr[i,] = sample(c(rep(0,4), rep(1,4)), 8)

# Look at the assignment for a single randomization
print(cbind(name, female, age, Z.cr[551,], y0, y1))
mean(y1[Z.cr[551,]==1]) - mean(y0[Z.cr[551,]==0])

# Calculate the mean difference between treated/control for each randomization
tau_sate.cr = rep(NA, n_reps)
for (i in 1:n_reps)
  tau_sate.cr[i] = mean(y1[Z.cr[i,]==1]) - mean(y0[Z.cr[i,]==0])

m.cr = round(mean(tau_sate.cr), 2)
sd.cr = round(sd(tau_sate.cr), 2)
hist(tau_sate.cr, main = paste("CR: Mean = ", m.cr, " SD = ", sd.cr, sep=""))

