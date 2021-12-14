##----------------------------------------------------
## EXERCISE 1
## Andrea Álvarez Pérez
##----------------------------------------------------

## Replicate results 
set.seed(123)

## Parameters
days <- 10000 # cambiar a 10000
replicates <- 100000 # cambiar a 100000
N <- c(100,300,1000,3000) # cambiar a c(100, 300, 1000, 3000)
bottleneck <- c(0.25,0.5, 1) # cambiar a 0.25, 0.5, 1

## Results:
# Time fixation matrix
tfix <- matrix(0,length(N),length(bottleneck))

# Probability of fixation matrix
pfix <- matrix(0,length(N),length(bottleneck))

##----------------------------------------------------
## MAIN CODE
##----------------------------------------------------

# Iterate through population size vector
for (i in 1:length(N)){
  
  # Fill matrix by rows
  # fixation time for pop i
  tfix_i <- c()

  # fixation probability for pop i
  pfix_i <- c()
  
  for (j in 1:length(bottleneck)){
    
    cat ("\nPopulation size: ", N[i], "\tBottleneck: ", bottleneck[j])
    cat ("\nReplicates running...")
    
    if ((N[i] == 1000 & bottleneck[j] == 1) | (N[i] == 3000)){
      cat("\nThis may take a while, please wait...")
    }
    
    # number of days the allele takes to get fixed
    days_fix <- c() 
    
    for (rep in 1:replicates){

      # initialize pop: all pop with 0
      pop <- rep(0,N[i])
      # seed my pop with just one mutation in a random individual
      random <- 1:N[i]
      seed <- sample(random, 1)
      pop[seed] <- 1
      indiv_100_rep <- c()
      indiv_1000_rep <- c()

      for (day in 1:days){
        #offspring sampling
        offspring <- sample(1:N[i], N[i]*bottleneck[j], replace=TRUE) 
        
        # update population vector for the offspring
        # I don't sample all the population, but I need a N population size
        # I fill the pop array by concatenating the same array 1/bottleneck[j] times
        pop <- rep(pop[offspring], 1/bottleneck[j])
        
        # with my population in that day, I calculate time and probability of fixation
        # the allele gets fixes when all my pop is composed by 1
        if (all(pop == 0)){ # no fixation
          break
        }else{
          if (all(pop == 1)){ # fixation
            # append days that take for allele fixation under certain conditions 
            days_fix <- append(days_fix, day)
            # Number of replicates the allele gets fixed
            num_fix <- length(days_fix)
            
            # we actualize the probability of fixation 
            # If during all 10000 replicates the allele became fixed X times, the probability
            # of fixation of that allele is X/10000
            prob <- num_fix/replicates
            
            break}
        }

      }
      
    }

    if (exists("prob") == TRUE){
      cat ("\nReplicates done!")
      # mean of days to take to fixation of all replicates + global probability of fixation
      cat("\nAllele gets fixed ", num_fix," times from all ",replicates," replicates")
      cat("\nTime of fixation:", mean(days_fix), "Probability of fixation: ", prob, "\n")
      }
    
    tfix_i<-append(tfix_i,mean(days_fix))
    pfix_i<-append(pfix_i,prob)

  }
  
  # Append results to final matrix
  tfix[i,]<-tfix_i
  pfix[i,]<-pfix_i
  
}

cat ("\nScript runned succesfully, please print results")
rownames(tfix) <- N
colnames(tfix) <- bottleneck

rownames(pfix) <- N
colnames(pfix) <- bottleneck


##----------------------------------------------------
## TABLE REPRESENTATION
##----------------------------------------------------

library(gt)
Population <- c(100, 300, 1000, 3000)
tfix2 <- cbind(Population, tfix)
tfix_data <- gt(as.data.frame(tfix2))
tfix_tbl <- 
  tfix_data %>%
  tab_header(
    title = md("**Time of fixation**")
  ) %>%
  tab_spanner(
    label = md("**Bottleneck**"),
    columns = 2:4
  )
tfix_tbl

Population <- c(100, 300, 1000, 3000)
pfix2 <- cbind(Population, pfix)
pfix_data <- gt(as.data.frame(pfix2))
pfix_tbl <- 
  pfix_data %>%
  tab_header(
    title = md("**Probability of fixation**")
  ) %>%
  tab_spanner(
    label = md("**Bottleneck**"),
    columns = 2:4
  )
pfix_tbl

##----------------------------------------------------
## DATA REPRESENTATION
##----------------------------------------------------

library(ggplot2)
library(reshape2)
library(plyr)
cols = c("red", "darkgreen", "blue")
cols2 = c("darkorange", "darkgreen", "turquoise", "purple")

# Estimate the relationships between tfix and pfix with population and bottleneck size

##----------------------------------------------------
## TFIX VS N

## LINEAR MODEL
# Ref: https://www.datacamp.com/community/tutorials/linear-regression-R
# lm([target variable] ~ [predictor variables], data = [data source])

tfixm <- as.data.frame(tfix)
df<-data.frame(Population=N,"0.25"=tfixm[,1],"0.5"=tfixm[,2],"1"=tfixm[,3])
df.long<-melt(df,id.vars="Population")
colnames(df.long) = c("population", "bottleneck", "Tfix")
head(df.long)
r2 <- c()
# calculate R2 values based on the linear model
for (i in 1:length(bottleneck)){
  linearmodel<-lm(tfix[,i]~N)
  r2 <- c(r2, round(summary(linearmodel)$r.square, 4))}

r2 <- data.frame(bottleneck, r2)

# Graphic representation
ggplot(df.long,aes(population, Tfix, color=factor(bottleneck)))+
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
  ggtitle("Fixation Time VS Population") +
  labs(colour = "Bottleneck") +
  geom_text(data = r2, aes(x = c(2700, 2700, 2300),  y = c(1000, 2700, 4400), 
                           label = paste("R2:", r2)), colour = cols, show.legend = FALSE)

## LOGLINEAR MODEL

dflog<-data.frame(logPopulation=log(N),"0.25"=log(tfixm[,1]),"0.5"=log(tfixm[,2]),"1"=log(tfixm[,3]))
dflog.long<-melt(dflog,id.vars="logPopulation")
colnames(dflog.long) = c("LnPopulation", "bottleneck", "LnTfix")
head(dflog.long)
r2 <- c()

for (i in 1:length(bottleneck)){
  loglinearmodel<-lm(log(tfix[,i])~log(N))
  r2 <- c(r2, round(summary(loglinearmodel)$r.square, 4))}

r2 <- data.frame(bottleneck, r2)

ggplot(dflog.long, aes(LnPopulation, LnTfix, color=factor(bottleneck)))+
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
  ggtitle("Ln(Fixation Time) VS Ln(Population)") +
  labs(colour = "Bottleneck") +
  geom_text(data = r2, aes(x = 6.5,  y = c(5.5, 6.3, 7.5), 
                           label = paste("R2:", r2)), colour = cols, show.legend = FALSE)

##----------------------------------------------------
## PFIX VS POPULATION

## LINEAR MODEL
pfixm <- as.data.frame(pfix)

df2<-data.frame(Population=N,"0.25"=pfixm[,1],"0.5"=pfixm[,2],"1"=pfixm[,3])
df2.long<-melt(df2,id.vars="Population")
colnames(df2.long) = c("population", "bottleneck", "Pfix")
head(df2.long)
r2_pfix <- c()
# calculate R2 values based on the linear model
for (i in 1:length(bottleneck)){
  linearmodel<-lm(pfix[,i]~N)
  r2_pfix <- c(r2_pfix, round(summary(linearmodel)$r.square, 4))}

r2_pfix <- data.frame(bottleneck, r2_pfix)

# Graphic representation
ggplot(df2.long,aes(population, Pfix, color=factor(bottleneck)))+
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
  ggtitle("Fixation Probability VS Population") +
  labs(colour = "Bottleneck") +
  geom_text(data = r2_pfix, aes(x = c(2700, 400, 400),  y = c(0.001, 0.007, 0.004), 
                           label = paste("R2:", r2_pfix)), colour = cols, show.legend = FALSE)

## LOGLINEAR MODEL
dflog2<-data.frame(LnPopulation=log(N),"0.25"=log(pfixm[,1]),"0.5"=log(pfixm[,2]),"1"=log(pfixm[,3]))
dflog2.long<-melt(dflog2,id.vars="LnPopulation")
colnames(dflog2.long) = c("LnPopulation", "bottleneck", "LnPfix")
head(dflog2.long)
r2_pfix <- c()

for (i in 1:length(bottleneck)){
  loglinearmodel<-lm(log(pfix[,i])~log(N))
  r2_pfix <- c(r2_pfix, round(summary(loglinearmodel)$r.square, 4))}

r2_pfix <- data.frame(bottleneck, r2_pfix)

ggplot(dflog2.long, aes(LnPopulation, LnPfix, color=factor(bottleneck)))+
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
  ggtitle("Ln(Fixation Probability) VS Ln(Population)") +
  labs(colour = "Bottleneck") +
  geom_text(data = r2_pfix, aes(x = 7.5,  y = c(-7, -7.5, -8), 
                           label = paste("R2:", r2_pfix)), colour = cols, show.legend = FALSE)

##----------------------------------------------------
## TFIX VS BOTTLENECK

# We have to traspose the tfixm matrix so its easier to operate with it
df_transpose <- data.frame(t(tfixm))
df_transpose

## LINEAR MODEL
df_BN<-data.frame(Bottleneck=bottleneck,"100"=df_transpose[,1],"300"=df_transpose[,2],"1000"=df_transpose[,3],"3000"=df_transpose[,4])
df_BN.long<-melt(df_BN,id.vars="Bottleneck")
colnames(df_BN.long) = c("bottleneck", "population", "Tfix")
head(df_BN.long)
r2 <- c()
# calculate R2 values based on the linear model
for (i in 1:length(N)){
  linearmodel<-lm(df_transpose[,i]~bottleneck)
  r2 <- c(r2, round(summary(linearmodel)$r.square, 4))}

r2 <- data.frame(N, r2)

# Graphic representation
ggplot(df_BN.long,aes(bottleneck, Tfix, color = factor(population)))+
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
  ggtitle("Fixation Time VS Bottleneck") +
  labs(colour = "Population size") +
  geom_text(data = r2, aes(x = 0.9,  y = c(0, 900, 1900, 3700), 
                           label = paste("R2:", r2)), colour = cols2, show.legend = FALSE)

## LOGLINEAR MODEL
dflog_BN<-data.frame(LnBottleneck=log(bottleneck),"100"=log(df_transpose[,1]),"300"=log(df_transpose[,2]),"1000"=log(df_transpose[,3]),"3000"=log(df_transpose[,4]))
dflog_BN.long<-melt(dflog_BN,id.vars="LnBottleneck")
colnames(dflog_BN.long) = c("LnBottleneck", "population", "LnTfix")
head(dflog_BN.long)
r2log <- c()

for (i in 1:length(N)){
  loglinearmodel<-lm(log(df_transpose[,i])~log(bottleneck))
  r2log <- c(r2log, round(summary(loglinearmodel)$r.square, 4))}

r2log <- data.frame(N, r2log)

ggplot(dflog_BN.long, aes(LnBottleneck, LnTfix, color=factor(population)))+
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
  ggtitle("Ln(Fixation Time) VS Ln(Bottleneck)") +
  labs(colour = "Population size") +
  geom_text(data = r2log, aes(x = -0.2,  y = c(4.5, 5.75, 6.75, 8), 
                           label = paste("R2:", r2log)), colour = cols2, show.legend = FALSE)

##----------------------------------------------------
## PFIX VS BOTTLENECK

## LINEAR MODEL
df_transpose_pfix <- data.frame(t(pfixm))
df_transpose_pfix

## LINEAR MODEL
df2_BN<-data.frame(Bottleneck=bottleneck,"100"=df_transpose_pfix[,1],"300"=df_transpose_pfix[,2],"1000"=df_transpose_pfix[,3],"3000"=df_transpose_pfix[,4])
df2_BN.long<-melt(df2_BN,id.vars="Bottleneck")
colnames(df2_BN.long) = c("bottleneck", "population", "Pfix")
head(df2_BN.long)
r2_pfix <- c()
# calculate R2 values based on the linear model
for (i in 1:length(N)){
  linearmodel<-lm(df_transpose_pfix[,i]~bottleneck)
  r2_pfix <- c(r2_pfix, round(summary(linearmodel)$r.square, 4))}

r2_pfix <- data.frame(N, r2_pfix)

# Graphic representation
ggplot(df2_BN.long,aes(bottleneck, Pfix, color = factor(population)))+
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
  ggtitle("Fixation Probability VS Bottleneck") +
  labs(colour = "Population size") +
  geom_text(data = r2, aes(x = 0.9,  y = c(0.008, 0.005, 0.0015, 0.0005), 
                           label = paste("R2:", r2)), colour = cols2, show.legend = FALSE)

## LOGLINEAR MODEL
dflog2_BN<-data.frame(LnBottleneck=log(bottleneck),"100"=log(df_transpose_pfix[,1]),"300"=log(df_transpose_pfix[,2]),"1000"=log(df_transpose_pfix[,3]),"3000"=log(df_transpose_pfix[,4]))
dflog2_BN.long<-melt(dflog2_BN,id.vars="LnBottleneck")
colnames(dflog2_BN.long) = c("LnBottleneck", "population", "LnPfix")
head(dflog2_BN.long)
r2log_pfix <- c()

for (i in 1:length(N)){
  loglinearmodel<-lm(log(df_transpose_pfix[,i])~log(bottleneck))
  print(summary(loglinearmodel))
  r2log_pfix <- c(r2log_pfix, round(summary(loglinearmodel)$r.square, 4))}

r2log_pfix <- data.frame(N, r2log_pfix)
r2log_pfix
ggplot(dflog2_BN.long, aes(LnBottleneck, LnPfix, color=factor(population)))+
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
  ggtitle("Ln(Fixation Probability) VS Ln(Bottleneck)") +
  labs(colour = "Population size") +
  geom_text(data = r2log_pfix, aes(x = -0.25,  y = c(-4.9, -5.8, -7.5, -8.5), 
                              label = paste("R2:", r2log_pfix)), colour = cols2, show.legend = FALSE)

