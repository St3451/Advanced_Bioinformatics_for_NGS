library(R.oo) #package with charToInt function to convert ascii to a value

r <- read.delim("MTnice.pileup",as.is=T,head=F, quote="")
bases <- strsplit(r[,5],"")
quality <- strsplit(r[,6],"")
fun <- function(x){
  y <- R.oo::charToInt(x)-33 #offset
  10^(-y/10)
}
quality <- lapply(quality,fun)
dat <- read.delim("MTnice.pileup",as.is=T,comment.char="",head=F, quote="")
names(dat) <- c("CHR","POS","REF",c("depth","bases","Qscore"))

## bases for individual 1 as a list
bases <- strsplit(toupper(dat$bases),"")
## ascii qualities for individual 1
asciiQ <- strsplit(dat$Qscore,"")
## quality values values
Q <- lapply(asciiQ,function(x) R.oo::charToInt(x) - 33 )


## Implementation 

# Function to compute the probabilities of the base given the genotype
get_prob <- function(genotype, bases, quality){
  ifelse(genotype==bases, 1-quality, quality/3)
}
prob_vec <- Vectorize(get_prob, vectorize.args = c("bases","quality"))

# Expectation maximization for one position
EMfreqStep <- function(p, bases_prob, theta_zero=rep(0.25, 4)){
  # Initialization
  f <- rep(0)
  fTemp <- theta_zero
  c <- 1
  # Iterate until convergence
  while(any(abs(f-fTemp) > 0.00001)){
    f <- fTemp
    # Q-Step (compute the posterior probabilities for each read (q_i))
    qi_den <- c(bases_prob$A[[p]]*f[1] + bases_prob$C[[p]]*f[2] + bases_prob$G[[p]]*f[3] + bases_prob$T[[p]]*f[4])
    qi_A <- bases_prob$A[[p]]*f[1] / qi_den
    qi_C <- bases_prob$C[[p]]*f[2] / qi_den
    qi_G <- bases_prob$G[[p]]*f[3] / qi_den
    qi_T <- bases_prob$T[[p]]*f[4] / qi_den
    # M-Step (compute the maximum likelihood)
    ml_den <- sum(qi_A + qi_C + qi_G + qi_T)
    ml <- c(sum(qi_A), sum(qi_C), sum(qi_G), sum(qi_T)) / ml_den
    fTemp <- ml
    c <- c+1
  }
  return (fTemp)
}

# Base probability for each position
bases_prob <- list(A=prob_vec(genotype="A", bases, quality), 
                   C=prob_vec(genotype="C", bases, quality),
                   G=prob_vec(genotype="G", bases, quality),
                   T=prob_vec(genotype="T", bases, quality))

# Run EM algorithm to estimate allele frequency at each position
thetas <- data.frame(A=numeric(0), C=numeric(0), G=numeric(0), T=numeric(0))
for (i in 1:length(bases)){
  thetas[i,] <- EMfreqStep(i, bases_prob)
  if (i%%1000==0){
    print(paste("Position", i))
  }
}

## Deliverables 

# Allele frequencies data frame
thetas
# Number of sites with most common allele frequency less than 0.9
sum(apply(thetas, 1, max) < 0.9)
# Position of these sites
which(apply(thetas, 1, max) < 0.9)
# Allele frequency across all sites
tot_freq <- apply(thetas, 2, mean)
tot_freq


