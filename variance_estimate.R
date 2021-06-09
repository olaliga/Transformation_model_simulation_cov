package <- c("dplyr", "foreach", "doParallel", "ggplot2", 'rootSolve')
lapply(package, require, character.only = T)

### bug : sigma
parameters <- expand.grid(
  r = c(0, 1, 1.5, 3, 6),
  a_star = c(0.7, 1.2), 
  v_star = c(0.7, 1.2)
)





##### 
simulate <- function(sample_size, beta_0, inform_t, ph1_real, r_star,  a_star, v_star){
  ### data generation
  x_1 <- rnorm(sample_size)
  #x_2 <- rbinom(sample_size, 1, 0.5)
  uni <- runif(sample_size)
  c <- runif(sample_size, 0, 10)
  
  ###
  eps_cox <- cumulate_inverse_cox(x =  uni, r = r_star)
  ht <- - beta_0[1]*x_1 + eps_cox
  t <- h_inver(t = ht, a = a_star, v = v_star)
  
  
  y <- c()
  delta <- c()
  for(i in c(1:sample_size)){
    y[i] <- min(c[i], t[i])
    if(t[i] <= c[i]){ delta[i] = 1}else{delta[i] = 0}
  }
  
  data_d <- data.frame(y = y , delta = delta, x1 = x_1)
  data_d <- data_d[order(data_d$y), ]
  
  
  ### process
  uncensored <- data_d$y[which(data_d$delta == 1)]
  
  sort_time <- sort(data_d$y)
  sort_uncensored <- sort(uncensored)
  
  risk_process_table <- data.frame(matrix(0, nrow = length(data_d$y), ncol = length(uncensored)))
  
  for( i in c(1:dim(risk_process_table)[2])){
    risk_process_table[min(which(sort_time == sort_uncensored [i])):dim(risk_process_table)[1], i] <- 1
  }
  
  
  count_process_table <- data.frame(matrix(1, nrow = length(data_d$y), ncol = length(uncensored)))
  
  for( i in c(1:dim(count_process_table)[2])){
    count_process_table[(max(which(sort_time == sort_uncensored [i])) + 1) :dim(count_process_table)[1], i] <- 0
  }
  
  count_process_table <- count_process_table[- (sample_size+1),] 
  
  for( i in c(1:dim(count_process_table)[2])){
    count_process_table[i] <- (as.matrix(data_d$delta)) * as.matrix(count_process_table[i])   
  }
  
  ############ 以下 check ############
  #### R_hat
  x_df <- data.frame(data_d[,3:dim(data_d)[2]])
  exp_beta_x <- exp(as.matrix(x_df) %*% beta_0)
  
  print("point 1")
  ##### init r 
  find_init_r <- function(t){
    G_array <- G(t = t*exp_beta_x, r = r_star )
    dN <- sum(count_process_table[,2])  - sum(count_process_table[,1]) 
    need <- t(as.matrix(risk_process_table[,1])) %*% G_array - dN
    return(need)
  }
  
  if(r_star == 0){
    lower = -10^5
  }else{
    lower = max(-(1-0.001)/(r_star*exp_beta_x))
  }
  
  rootSearch = uniroot(find_init_r, lower = lower, upper = 10^5)
  init_r <- rootSearch$root
  
  print("point 2")
  #### all r vector
  r <- c(init_r)
  ### i = 2 to dim(count_process_table)[2] 
  
  delta_r_i_vector <- c()
  
  for(i in c(2 : dim(count_process_table)[2])){
    g_array <- g(r[i-1]*exp_beta_x, r = r_star)
    delta_i <- sum(count_process_table[,i])  - sum(count_process_table[,i-1]) 
    frac_under <- t( as.matrix(risk_process_table[,i])*exp_beta_x ) %*% g_array 
    delta_r_i <- delta_i/frac_under
    
    
    delta_r_i_vector[i-1] <- delta_r_i[1,1]
    
    r[i] <- r[i-1] + delta_r_i[1,1]
  }
  
  print("point 3")
  ### score function 
  ### score function 
  
  s_0 <- function(j){
    inside_g <- r[j]*exp_beta_x
    need <- t(risk_process_table[,j]*exp_beta_x) %*% g(t = inside_g, r = r_star )  
    
    return(need)
  }
  
  
  ####new
  partial_r_df <- data.frame(x_1 = rep(0, 1), 
                             x_2 = rep(0, 1))
  
  ##########################
  find_init_partial_r <- function(init){
    dum <- t(t(r[1]*data_d[,3:4])+ init)
    
    need <- t(as.matrix(risk_process_table[,1])*(g(t = r[1]*exp_beta_x, r = r_star )*exp_beta_x))%*%dum
    
    return(sum(t(need)))  
  }
  
  lower = max(-(1-0.001)/(r_star*exp_beta_x))
  rootSearch = uniroot(find_init_r, lower = 10^-5, upper = 10^5)
  partial_r_df[,1] <- rootSearch$root
  print("point 4")
  
  # init_for_root <- c(0.1, 0.1)
  # rootSearch = optim(par = init_for_root, find_init_partial_r, method = 'BFGS', hessian=TRUE)
  # init_partial_r <- rootSearch$par
  # 
  # partial_r_df[,1] <- init_partial_r
  
  part_r_generate <- function(times){
    delta_i <- sum(count_process_table[,times])  - sum(count_process_table[,times-1]) 
    need <- partial_r_df[,times - 1] - delta_i*s_1(times - 1)/(s_0(times - 1))^2
    
    return(need)
  } 
  
  
  psi <- function(i, times){
    inside_g <- r[times]*exp_beta_x[i]
    dum <- r[times]*data_d[i,3] + partial_r_df[, times]
    
    need_1 <- as.numeric(data_d[i, 3] + g_dot(t = inside_g, r = r_star )/g(t = inside_g, r = r_star )*exp_beta_x[i]*dum[1])
    
    need <- c(need_1)
    return(need)
  }
  

  s_1 <- function(j){
    psi_dum <- t(sapply(1:dim(risk_process_table)[1], function(t){psi(t, j)}))
    inside_g <- r[j]*exp_beta_x
    need <- t(risk_process_table[,j]*g(t = inside_g, r = r_star )*exp_beta_x) %*% t(psi_dum)
    
    return(as.vector(need))
  }

  s_2 <- function(j){
    inside_g <- r[j]*exp_beta_x
    need_dum <- risk_process_table[,j]*g(t = inside_g, r = r_star )*exp_beta_x
    
    need <- matrix(rep(0,  length(psi(1,j))^2), nrow =  length(psi(1,j)))
    
    for(i in c(1:dim(risk_process_table)[1])){
      need <- need + psi(i, j)%*%t(psi(i, j))*need_dum[i]
      
    }
    
    return(need)
  }


  #### 久
  for(i in c(2:dim(count_process_table)[2])){
    partial_r_df[,i] <- part_r_generate(i) 
  }
  print("point 5")
  
  
  psi <- Vectorize(psi, vectorize.args = c('i', 'times'))
  
  inside_score_s <- function(j){ return( s_1(j)/s_0(j) ) }
  inside_score_s <- Vectorize(inside_score_s)
  
  inside_score_gamma <- function(j){ return( -s_2(j)/as.numeric(s_0(j))  + (s_1(j)/as.numeric(s_0(j)))%*%t((s_1(j)/as.numeric(s_0(j))))    ) }
  inside_score_gamma <- Vectorize(inside_score_gamma)
  #need_3 = inside_score_gamma(c(1:dim(count_process_table)[2])) ## 久
  #gamma_est = matrix(c(sum(need_3[1,]),sum(need_3[2,]),sum(need_3[3,]),sum(need_3[4,])), nrow = 2)
  
  score <- c(0, 0)
  num <- 0
  uncensored_index <- which(data_d$delta == 1)
  need_1 = inside_score_s(c(1:dim(count_process_table)[2])) ## 久

  print("point 6")
  score <- c(0)
  #dni
  inside_need_mid = function(i, u){
    if(u > 1){
      need = count_process_table[i, u]-count_process_table[i, u-1] 
    }else{
      need = count_process_table[i, u]
    }
    return(need)   
  }
  
  time = 0
  for(i in uncensored_index ){
    time = time + 1
    need_2 = psi(i, c(1:dim(count_process_table)[2]) )
    need_mid = need_2 - need_1
    dni <- sapply(uncensored_index, inside_need_mid, u = time)
    
    need <- c(sum(need_mid[1]*dni))
    score <- score + need
    
    #print(paste("score : ", i))
  }
  
  score = score/dim(count_process_table)[1]  ######## need 1
  #gamma_est = gamma_est/dim(count_process_table)[1]
  
  
  #print(paste("r_star_2:", r_star))
  print("point 7")
  
  ### auxiliary information
  R_step <- stepfun(x = sort_uncensored, y = c(0,r), right = F)
  print("point 8")
  ## k = 1 information X1 > 0, t = 2
  condition <- rep(0, dim(data_d)[1] )
  condition[which(data_d$x1 > 0)] <- 1
  inform <-  (t(condition) %*% (exp(-G(t = R_step(inform_t)*exp_beta_x, r = r_star)) - ph1_real) )/dim(count_process_table)[1]
  print("point 9")
  
  
  ##### var ######
  
  
  d_r <- function(s){
    if(s != 1){
      need <- (r[s]- r[s-1]) 
    }else{
      need <- (r[s]) 
    }
  }
  
  d_r <- Vectorize(d_r)
  
  v_0 <- function(j){
    inside_g <- r[j]*exp_beta_x
    need_dum <- risk_process_table[,j]*g_dot(t = inside_g, r = r_star )*exp_beta_x^2
    need <- sum(need_dum)/sample_size
    return(need)
  }
  
  v_1 <- function(j){
    inside_g <- r[j]*exp_beta_x
    need_dum <- risk_process_table[,j]*g_dot(t = inside_g, r = r_star )*exp_beta_x^2
    
    need <- c(0)
    
    for(i in c(1:dim(risk_process_table)[1])){
      need <- need + psi(i, j)*need_dum[i]
      
    }
    
    return(need)
  }
  
  v_2 <- function(j){
    inside_g <- r[j]*exp_beta_x
    need_dum <- risk_process_table[,j]*g_dot(t = inside_g, r = r_star )*exp_beta_x^2
    
    need <- c(0)
    
    for(i in c(1:dim(risk_process_table)[1])){
      need <- need + psi(i, j)%*%t(psi(i, j))*need_dum[i]
      
    }
    
    return(need)
  }
  
  inside_e_star <- function(j){
    if(j != 1){
      need <- v_0(j)/s_0(j)*(r[j]- r[j-1]) 
    }else{
      need <- v_0(j)/s_0(j)*(r[j]) 
    }
    return(need)
  }
  
  inside_e_star <- Vectorize(inside_e_star)
  
  
  e_star_dum <- function(s){
    inside_exp <- sum(inside_e_star(c(1:s)))
    
    need <- exp(inside_exp)
    
    return(need)
  }
  
  e_star_dum <- Vectorize(e_star_dum)
  
  ## 久 
  e_star_table <- e_star_dum(1:dim(count_process_table)[2])
  print("point 10")
  e_star <- function(s){e_star_table[s]}
  
  
  k <- function(u, s){
    if(u<=s){
      need <- e_star(s)^(-1)*(e_star(u)/s_0(u))
    }else{
      need <- 0
    }
    
    return(need)
  }
  
  k <- Vectorize(k, vectorize.args = "u")
  
  q_star_inside <- function(s){
    inside <- function(u){k(s = s, u = u)*(v_1(u)*s_0(u) - s_1(u)*v_0(u))*d_r(u)}
    inside <- Vectorize(inside)
    need1 <- sum(inside(1:dim(count_process_table)[2]))
    need <- c(need1)
    
    return(need)
  }
  
  init_time = proc.time()
  q_star_table <- sapply(1:dim(count_process_table)[2], q_star_inside)
  q_star <- function(s){q_star_table[s]}
  total_time = proc.time() - init_time
  print("point 11") ## 240 seconds
  dn <- function(i, s){
    if(s == 1){
      need <- count_process_table[i, s]
    }else{
      need <- count_process_table[i, s]-count_process_table[i, s-1]
    }
    return(need)
  }
  
  
  dm <- function(i, s){
    inside_g <- r[s]*exp_beta_x[i]
    need <- dn(i, s) - risk_process_table[i, s]*g(inside_g, r = r_star)*d_r(s)*exp_beta_x[i]
    
    return(need)
  }
  
  
  sigma_inside <- function(i, s){
    need <- (psi(i = i, times = s) - s_1(s)/as.numeric(s_0(s)) - q_star(s) )*dm(i = i, s = s) 
    return(need)
  }
  
  
  sigma_inside <- Vectorize(sigma_inside, vectorize.args = "s")
  
  sigma_inside_int <- function(i){
    need_dum <- sapply(1:dim(count_process_table)[2],function(s){sigma_inside(i, s = s)})
    need1 <- sum(need_dum)
    
    need = c(need1)
    return(need)
  }
  
  ## ksi_i = sigma_inside_int(i)
  ksi_inside <- function(t){
    need <- sigma_inside_int(t)
    print(t)
    return(need)
  } 
  
  ksi <- sapply(1:dim(count_process_table)[1], ksi_inside)
  print("point 12")
  
  ## inform using true r0 = a_star*x^v_star

  r0 <- function(x){
    return(a_star*x^v_star)
  }
  
  informi <- exp(-G(t = R_step(inform_t)*exp_beta_x, r = r_star)) - ph1_real
  j1  <- exp(-G(t = r0(inform_t)*exp_beta_x, r = r_star)) - ph1_real 
  j2_resid <- informi - j1
  
  g_decomp <- data.frame(ksi = ksi, j1 = j1, j2_resid = j2_resid)
  
  simulation_result <- list(score = score, inform = inform, var = var(g_decomp))
  print("point end")
  return(simulation_result)
} 

####### function generate

### generate eps

cumulate_inverse_cox <- function(x, r){
  need <- log(G_inverse(r = r, t = -log(1-x)))
  return(need)
}


### a v specify
h_inver <- function(a, v, t){
  need <- (exp(t)/a)^(1/v)
  return(need)
}

h_inver <-Vectorize(h_inver)


### g function 

G <- function(t, r){
  
  if(r == 0){
    need <- t
  }else{
    need <- log(abs(1+r*t))/r
  }
  
  return(need)
}

g <- function(t, r){
  need = 1/(1+r*t)
  return(need)
}
g <- Vectorize(g)

g_dot <- function(t, r){
  need = -r*(1 + r*t)^(-2)
  return(need)
}
g_dot <- Vectorize(g_dot)


G_inverse <- function(r, t){
  if(r == 0){
    need <- t
  }else{
    need <- (-1 + exp(r*t))/r
  }
  
  return(need)
}

#### real information

integrand <- function(x, t, r_star, a_star, v_star, beta_0) {
  need <- 2*exp(-G( t = a_star*t^(v_star)*exp(beta_0*x), r = r_star) )*1/sqrt(2*pi)*exp(-x^2/2)
  return(need)
}

# i = 1
# time = proc.time()
# try = simulate(sample_size = 100,
#                    beta_0 = c(2),
#                    inform_t = 6,
#                    ph1_real = integrate(integrand, lower = 0, upper = Inf, t = 6,
#                                         r_star = parameters[i,]$r, 
#                                         a_star = parameters[i,]$a_star, 
#                                         v_star = parameters[i,]$v_star, 
#                                         beta_0 = c(2))$value,
#                    r_star = parameters[i,]$r,
#                    a_star = parameters[i,]$a_star,
#                    v_star = parameters[i,]$v_star)
# 
# one_simulate_time = proc.time()-time
