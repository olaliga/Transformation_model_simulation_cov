packages_main <- c('doSNOW')
lapply(packages_main, require, character.only = T)
source('variance_estimate.R')

print("main")
parameters <- expand.grid(
  r = c(0, 1, 3, 6),
  a_star = c(0.7, 1.2), 
  v_star = c(0.7, 1.2)
)

trial_num = dim(parameters)[1]
final_output_score <- as.list(rep(0,trial_num))
final_output_inform <- as.list(rep(0,trial_num))
final_output_var <- as.list(rep(0,trial_num))

#########
time = proc.time()
#dim(parameters)[1]
trial_num = 2

for(i in c(1:trial_num)){
  cores = detectCores()
  cl <- makeCluster(cores[1]-2)
  registerDoSNOW(cl)
  simulation_times <- 2
  pb <- txtProgressBar(max = simulation_times, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  print("")
  writeLines(paste("Total_Progress : ", i, "Percentage :",round(i/dim(parameters)[1], 2)))
  
  comb <- function(x, ...) {
    lapply(seq_along(x),
           function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
  }
  
  simulation_result <- foreach(j = 1:simulation_times, .combine='comb', .multicombine=TRUE,
                               .init=list(list(), list(), list()), .options.snow = opts) %dopar%
    {
      simulate(sample_size = 100,
                   beta_0 = c(2),
                   inform_t = 5,
                   ph1_real = integrate(integrand, lower = 0, upper = Inf, t = 5,
                                        r_star = parameters[i,]$r, 
                                        a_star = parameters[i,]$a_star, 
                                        v_star = parameters[i,]$v_star, 
                                        beta_0 = c(2))$value,
                   r_star = parameters[i,]$r,
                   a_star = parameters[i,]$a_star,
                   v_star = parameters[i,]$v_star)
    }
  
  score <- simulation_result[[1]]
  inform <- simulation_result[[2]]
  var <- simulation_result[[3]]
  
  final_output_score[[i]] <- score
  final_output_inform[[i]] <- inform
  final_output_var[[i]] <- var
    
  close(pb)
  stopCluster(cl) 
}

total_time = proc.time() - time


#########
# cores = detectCores()
# cl <- makeCluster(cores[1]-2)
# registerDoSNOW(cl)
# 
# comb <- function(x, ...) {
#   lapply(seq_along(x),
#          function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
# }
# 
# oper <- foreach(i=1:10, .combine='comb', .multicombine=TRUE,
#                 .init=list(list(), list(), list())) %dopar% {
#                   list(c(3,4), matrix(c(1,0,0,1), nrow = 2), matrix(c(1,2,2,1), nrow = 2))
#                 }
# 
# oper1 <- oper[[1]]
# oper2 <- oper[[2]]
# oper3 <- oper[[3]]
# stopCluster(cl) 

