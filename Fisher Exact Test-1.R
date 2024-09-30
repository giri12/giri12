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


# ---- Completely Randomized Design

# Repeat the completely randomized trial to get the randomization distribution
#  note: could enumerate every possible Z.bern, but we will simulate
#        more 'reps' than necessary as a general procedure
n_reps = 10000
Z.cr = matrix(NA, n_reps,n) # a matrix to store all of the different Z vectors
for (i in 1:n_reps)
  Z.cr[i,] = sample(c(rep(0,4), rep(1,4)), 8)


# - Illustrate Fisher's Radomization Test
#   For the completely randomized experiment

# Calculate the test statistic under the observed data
T_obs = mean(y[treatment==1]) - mean(y[treatment==0])
T_obs

# - for each simulated Z.cr, caluclate the test statistic under the sharp null hypothesis
T_sim = rep(NA, n_reps)
for (i in 1:n_reps)
  T_sim[i] = mean(y[Z.cr[i,]==1]) - mean(y[Z.cr[i,]==0])

hist(T_sim)
mean(abs(T_sim) >= T_obs) #This is the Fisher exact P-value
