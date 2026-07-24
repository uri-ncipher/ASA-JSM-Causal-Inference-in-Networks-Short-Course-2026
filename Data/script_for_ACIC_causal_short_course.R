
# first set working directory to where this script is in

rm(list = ls())
setwd("~/URI-PhD-local/Conferences/2026-ACIC/short_course")

### 10 components
# ncomps <- 10
# csize <- c(185, 9, 6, 3, 3, 2, 2, 2, 2, 2)
# avg_degree <- 3.35

### 20 components
ncomps <- 20
csize <- c(28, 18, 23, 19, 26, 12, 15, 19, 10, 7, 9, 8, 6, 3, 3, 2, 2, 2, 2, 2)
avg_degree <- 3.35



# load R packages
library(lme4)
library(plyr)
library(dplyr)
library(igraph)
library(numDeriv)
library(gtools)


##############################################
############  Define Functions  ##############
##############################################

### function for calculating neighbors' covariates
# avg_nbrs_covs <- function(net0, data, variable){
#   avg_var <- sapply(1:nrow(data), function(i) mean(subset(data, data$id %in% neighbors(net0, i))[[variable]]))
#   return(avg_var)
# }

num_neighbors=function(net, node, data){
  return(length(neighbors(net, node)))
}
trt_neighbors=function(net, node, data, treatment){
  return(sum(subset(data, data$id %in% neighbors(net, node))[[treatment]]))
}
avg_neighbors=function(net, node, data, variable){
  return(mean(subset(data, data$id %in% neighbors(net, node))[[variable]]))
}


### Y.bar(a, alpha) for individual i
y.bar.ind= function(PO, i, a, alpha){
  #######################
  ### Arguments: 
  # PO:    a list of matrix, each element in the list is a matrix, this matrix represents 
  #          all potential outcome for individual i;
  # i:     a numeric number, indicates the index of node/participant;
  # a:     individual's treatment. takes 0 or 1 when treatment is binary;
  # alpha: a numeric value, represents the allocation strategies. e.g. 0.25.
  #######################
  M=PO[[i]]
  A=a+1
  y=0
  for (j in 1:dim(M)[2]){
    y=y+M[A, j]*choose(dim(M)[2]-1, j-1)*alpha^(j-1)*(1-alpha)^(dim(M)[2]-j)
  }
  return(y)
}

### Y.bar(a, alpha) for all participants
y.bar=function(PO, a, alpha){
  #######################
  ### Arguments: 
  # PO:    a list of matrix, each element in the list is a matrix, this matrix represents 
  #          all potential outcome for individual i;
  # a:     individual's treatment. takes 0 or 1 when treatment is binary;
  # alpha: a numeric value, represents the allocation strategies. e.g. 0.25.
  #######################
  y=0
  for (i in 1:n){
    y=y+y.bar.ind(PO, i, a, alpha)
  }
  return(y/n)
}

### Y.bar(alpha) for individual i
y.bar_ind_margin=function(PO, i, alpha){
  #######################
  ### Arguments: 
  # PO:    a list of matrix, each element in the list is a matrix, this matrix represents 
  #          all potential outcome for individual i;
  # i:     a numeric number, indicates the index of node/participant;
  # alpha: a numeric value, represents the allocation strategies. e.g. 0.25.
  #######################
  M=PO[[i]]
  y=0
  for (k in 1:dim(M)[2]){
    y=y+M[1, k]*choose(dim(M)[2]-1, k-1)*alpha^(k-1)*(1-alpha)^(dim(M)[2]-k)*(1-alpha)
  }
  for (k in 1:dim(M)[2]){
    y=y+M[2, k]*choose(dim(M)[2]-1, k-1)*alpha^(k-1)*(1-alpha)^(dim(M)[2]-k)*(alpha)
  }
  return(y)
}

### Y.bar(alpha) for all participants
y.bar_margin=function(PO, alpha){
  #######################
  ### Arguments: 
  # PO:    a list of matrix, each element in the list is a matrix, this matrix represents 
  #          all potential outcome for individual i;
  # alpha: a numeric value, represents the allocation strategies. e.g. 0.25.
  #######################
  y=0
  for (i in 1:n){
    y=y+y.bar_ind_margin(PO, i, alpha)
  }
  return(y/n)
}

# direct effect  (DE = Y(1, alpha) - Y(0, alpha))
DE_true=function(alpha, PO){ # here alpha is a numeric number
  return(y.bar(PO, 1, alpha)-y.bar(PO, 0, alpha))
}

# indirect (spillover) effect   (IE = Y(0, alpha1) - Y(0, alpha0))
IE_true=function(alpha, PO){
  alpha1=alpha[1]
  alpha0=alpha[2]
  return(y.bar(PO, 0, alpha1)-y.bar(PO, 0, alpha0))
}

# total effect    (TE = Y(1, alpha1) - Y(0, alpha0))
TE_true=function(alpha, PO){
  alpha1=alpha[1]
  alpha0=alpha[2]
  return(y.bar(PO, 1, alpha1)-y.bar(PO, 0, alpha0))
}

# overall effect  (OE = Y(alpha1) - Y(alpha0))
OE_true=function(alpha, PO){
  alpha1=alpha[1]
  alpha0=alpha[2]
  return(y.bar_margin(PO, alpha1)-y.bar_margin(PO, alpha0))
}



# ################################################################################
# ############     Fixed coefficients (5 covariates, without age)   ##############
# ################################################################################
# # for generating outcomes
# intercept<- -3.18   # a negative value
# beta.a <- 0.67
# beta.alpha <- 1.53
# beta.a.alpha <- -2.03  # interaction between a and alpha
# 
# beta.var1 <- 0.12
# beta.var2 <- -0.37
# beta.var3 <- 1.62
# beta.var4 <- -0.01
# beta.var5 <- 0.59
# 
# beta.avg_var1 <- 0.29
# beta.avg_var2 <- -0.34
# beta.avg_var3 <- 1.26
# beta.avg_var4 <- 0.51
# beta.avg_var5 <- 0.05
# 
# # for propensity score model
# ps.intercept <- -1.97
# ps.var1 <- 0.02
# ps.var2 <- 0.49
# ps.var3 <- -0.25
# ps.var4 <- -0.53
# ps.var5 <- -0.46
  

################################################################################
############     Fixed coefficients (6 covariates, WITH age)   ##############
################################################################################
# for generating outcomes
intercept<- 0.88   # a negative value
beta.a <- 0.44
beta.alpha <- 1.42
beta.a.alpha <- -2.49  # interaction between a and alpha

beta.var1 <- 0.05
beta.var2 <- -0.63
beta.var3 <- 1.48
beta.var4 <- 0.12
beta.var5 <- 0.72
beta.var6 <- -0.07

# beta.avg_var1 <- 0.18
# beta.avg_var2 <- -0.30
# beta.avg_var3 <- 1.45
# beta.avg_var4 <- 0.46
# beta.avg_var5 <- -0.02
# beta.avg_var6 <- 0.01

# for propensity score model
ps.intercept <- 0.54
ps.var1 <- 0.01  # posneg
# ps.var2 <- 0.45  # date
ps.var2 <- 0.25 
ps.var3 <- -0.45 # share
# ps.var4 <- -0.40 # education 
ps.var4 <- -0.20 # education 
# ps.var5 <- -0.46 # employment
ps.var5 <- -0.11 # employment
ps.var6 <- -0.07 # age




# initial network. m=number of edges.
set.seed(6202)
net0=sample_gnm(n=csize[1], m=ceiling((csize[1] * avg_degree)/2), directed = FALSE) 
while (!is_connected(net0)) {
  net0=sample_gnm(n=csize[1], m=ceiling((csize[1] * avg_degree)/2), directed = FALSE)
}
length(V(net0))
mean(degree(net0))
components(net0)$no

### Add components to net0, in order to generate a network with specific number of components
for (i in 2:length(csize)){
  csize_temp = csize[i]
  print(c(i, csize[i]))
  tryCatch({
    if (csize_temp == 2){
      net_temp=sample_gnm(n=csize_temp, m=1, directed = FALSE) # n: The number of vertices in the graph; 
          # m: The number of edges in the graph.
    } else if (csize_temp == 3) {
      net_temp=sample_gnm(n=csize_temp, m=sample(c(2,3), 1), directed = FALSE)
    } else {
      net_temp=sample_gnm(n=csize_temp, m=ceiling((csize_temp * avg_degree)/2), directed = FALSE)
      while (!is_connected(net_temp)) {
        net_temp=sample_gnm(n=csize_temp, m=ceiling((csize_temp * avg_degree)/2), directed = FALSE)
      }
    }
    net0=net0+net_temp}, error=function(e){cat("ERROR", conditionMessage(e), "\n")})
}
n=length(V(net0))       # count number of nodes/individual/participant.
n
m=components(net0)$no
m
components(net0)
mean(degree(net0))

components_test <- cluster_fast_greedy(net0)$membership
m = max(components_test)
m

data=data.frame(id=1:n, num_nn=unlist(lapply(1:n, num_neighbors, net=net0)),  # num_nn: number of neighbors of each node (is also degree).
                component=components(net0)$membership)  # component: indicate which component each node belongs to.

##### create six baseline covariates: five binary, one continuous 
# set.seed(2026)
data$var1 <- rbinom(n, 1, 0.52) # based on TRIP data covariate `posneg.x`
data$var2 <- rbinom(n, 1, 0.51) # based on TRIP data covariate `date`
data$var3 <- rbinom(n, 1, 0.74) # based on TRIP data covariate `share.x`
data$var4 <- rbinom(n, 1, 0.39) # based on TRIP data covariate `edubin`
data$var5 <- rbinom(n, 1, 0.47) # based on TRIP data covariate `employbin`
data$var6 <- round(rnorm(n, 35.1, 8.19), 0) # based on TRIP data covariate `age` (based on 216 individuals). mean=35.1, sd=8.19.

summary(data$var1)

##### calculate averaged neighbors' covariates
data$avg_var1=unlist(lapply(1:n, avg_neighbors, net=net0, data=data, variable="var1"))
data$avg_var2=unlist(lapply(1:n, avg_neighbors, net=net0, data=data, variable="var2"))
data$avg_var3=unlist(lapply(1:n, avg_neighbors, net=net0, data=data, variable="var3"))
data$avg_var4=unlist(lapply(1:n, avg_neighbors, net=net0, data=data, variable="var4"))
data$avg_var5=unlist(lapply(1:n, avg_neighbors, net=net0, data=data, variable="var5"))
data$avg_var6=unlist(lapply(1:n, avg_neighbors, net=net0, data=data, variable="var6"))
data$avg_var6=round(data$avg_var6, 2)

summary(data$avg_var1)
  
##### generate all potential outcomes #####
outcomes.matrix.list<- list() # all potential outcomes
for (i in 1:n) {
  outcomes.matrix.list[[i]]<- matrix(data = NA, nrow = 2, ncol = data$num_nn[i] + 1) 
  # a list of matrix, each matrix contains 2 rows (realization of i: 0 or 1) and 5 columns (realization of nearest neighbors: {0, 1, 2, 3, 4})
}

# fill in the matrices
set.seed(2026)
for (i in 1:n) {
  for (a in 0:1) {  # 2 rows. indicating participant i receive treatment (a=1) or not (a=0).
    for (k in 0:data$num_nn[i]) { 
      alpha.avg <- k/data$num_nn[i]  # proportion of treated neighbors
      p=plogis(intercept  # baseline outcome status-the intercept; plogis: gives the logistic distribution function
               +(beta.a * a) #treatment effect independent of coverage.    (ai = a-1 = 0 or 1, because a=1 or 2)
               +(beta.alpha*(alpha.avg))  #indirect effect
               +(beta.a.alpha*a*(alpha.avg)) #treatment effect that scales with coverage
               +(beta.var1*data$var1[i] + beta.var2*data$var2[i] + beta.var3*data$var3[i] +
                    beta.var4*data$var4[i] + beta.var5*data$var5[i] + beta.var6*data$var6[i])
               # +(beta.avg_var1*data$avg_var1[i] + beta.avg_var2*data$avg_var2[i] + beta.avg_var3*data$avg_var3[i] +
               #   beta.avg_var4*data$avg_var4[i] + beta.avg_var5*data$avg_var5[i] + beta.avg_var6*data$avg_var6[i])
      )
      outcomes.matrix.list[[i]][a+1,k+1]<- rbinom(1, 1, prob=p)
    }
  }
}

## calculate true causal effects
alpha=c(0.5, 0.4, 0.3, 0.2) 
contrast=t(combn(alpha,2)) 
contrast

true_DE <- cbind(ldply(alpha, DE_true, PO=outcomes.matrix.list), 
                 alpha0=alpha, alpha1=alpha, type="Direct")
true_IE <- cbind(adply(contrast, 1, IE_true, PO=outcomes.matrix.list), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Indirect")[,-1] # use "[,-1]" to remove row index
true_TE <- cbind(adply(contrast, 1, TE_true, PO=outcomes.matrix.list), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Total")[,-1]
true_OE <- cbind(adply(contrast, 1, OE_true, PO=outcomes.matrix.list), 
                 alpha0=contrast[, 1], alpha1=contrast[, 2], type="Overall")[,-1]

## export true causal effects
true_CE_raw <- rbind(true_DE, true_IE, true_TE, true_OE)
names(true_CE_raw) <- c("true_effect", "alpha0", "alpha1", "type")
true_CE_raw2 <- true_CE_raw[order(true_CE_raw$type, true_CE_raw$alpha1, true_CE_raw$alpha0), ]
rownames(true_CE_raw2) <- NULL
true_CE <- true_CE_raw2[c(1:8, 10, 9, 17:20,22,21, 11:14,16,15), ] # reset rows in order to be consistent with row order in Table 7 in AOAS paper
true_CE
write.csv(true_CE, "./TRIPsim_true_causal_effects.csv", row.names = FALSE)


## assign the component-level random effect to data
raneff=rnorm(m, mean=0, sd=0.617)
data$raneff=raneff[data$component]

## assign treatment ("A") to data
data$treatment <- rbinom(
  n,1, plogis(ps.intercept+ ps.var1*data$var1 + ps.var2*data$var2 + ps.var3*data$var3 + 
                + ps.var4*data$var4 + ps.var5*data$var5 + ps.var6*data$var6 + data$raneff))

data$num_tnn= unlist(lapply(1:n, trt_neighbors, data=data, net=net0, treatment = "treatment"))  # number of exposed nearest neighbors.
data$num_utnn=data$num_nn-data$num_tnn    # number of unexposed nearest neighbors.

## assign outcome to data (="observed outcome")
y=c()
for (i in 1:n){
  s=data$treatment[i]+1  # row, s=1 or 2, where treatment is 0 or 1
  t=data$num_tnn[i]+1    # t=any{1, 2, 3 , 4, 5}
  y=c(y, outcomes.matrix.list[[i]][s, t])
}
data$outcome=y  # "y": observed outcome

# rename the nodes covariates 
data_renamed <- data %>%
  rename(
    posneg = var1, date = var2, share = var3, education = var4, employment = var5, age=var6,
    avg_posneg = avg_var1, avg_date = avg_var2, avg_share = avg_var3,
    avg_education = avg_var4, avg_employment = avg_var5, avg_age = avg_var6
  )

head(data_renamed)

# check distribution
table(data_renamed$treatment, data_renamed$share)
table(data_renamed$treatment, data_renamed$date)
table(data_renamed$treatment, data_renamed$posneg)
table(data_renamed$treatment, data_renamed$education)
table(data_renamed$treatment, data_renamed$employment)

## export edgelist
edges <- igraph::as_data_frame(net0, what = "edges")
write.csv(edges, "./TRIPsim_edges.csv", row.names = FALSE)

### export nodelist
write.csv(data_renamed, "./TRIPsim_nodes.csv", row.names = FALSE)
