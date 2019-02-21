#' Function to simuulate spatial process for shot chart data for Damian Lillard
#'
#' This function allows you to explore differnt "true" parameter values assocciated with the simulated spatial process.
#' @param sigma.sq.T "True" value for sigma.sq
#' @param tau.sq.T "True" value for tau.sq
#' @param phi.T "True" value for phi
#' @param grid Length of one side of the grid assocciated with the simulated spatial process
#' @keywords Spatial
#' @export
#' @examples
#' Spatial.Function()

Dame0.Function <- function(sigma.sq.T , tau.sq.T, phi.T , grid) {
  
  X <- c(24,  31,  34,  39,  41,  47,  48,  58,  62,  70  ,77,  81, 100, 103, 125, 138, 155, 165, 170, 172, 179, 185, 193, 199)
  x1 <- c(-24, -114,  228,   20,    6, -135,  165,  -15, -169,   63,   52,    3,   62,  135, -164,   40,  -20,   82 , -42,   48
          ,-10,  -36,  -25,   97)
  
  y1 <- c(35 ,235,   5, 244, -15, 198, 189,  -6, 193, 256, 242, 234, 242, 205, 209, 236,  -7, 229,  25, 234,  31 ,253,  28 ,245)
  
  Prop.Made <- c(0.2075472 ,0.4070796, 0.2307692, 0.3877551, 0.6498856, 0.3676471, 0.4090909, 0.6062603, 0.3956044, 0.2857143,
                 0.3157895, 0.4477612, 0.2692308, 0.4111111, 0.3552632, 0.3731343, 0.5510204, 0.2361111, 0.3214286, 0.3750000
                 ,0.3846154, 0.3448276, 0.3402062, 0.2800000)
  
  Distance.From.Basket.FT <- c(4.243819, 26.119150, 22.805482, 24.481830,  1.615549, 23.964349, 25.089041,  1.615549, 25.653460, 26.363801,
                               24.752374, 23.401923, 24.981593, 24.545875, 26.566332, 23.936583,  2.118962, 24.323857,  4.887740, 23.887235,
                               3.257299, 25.554843,  3.753665, 26.350332)
  
  SD <- data.frame(X, x1 , y1 , Prop.Made , Distance.From.Basket.FT)
  
  #set.seed(seed.number)
  
  library(mnormt)
  library(tidyr)
  library(dplyr)
  library(ggplot2)
  library(raster)
  
  sigma.sq <- sigma.sq.T
  tau.sq <- tau.sq.T
  phi <- phi.T
  
  #SD <- read.csv("Project 1 Final Data.csv" , header = T)
  
  ## Fit Model ##
  
  lm1 <- lm(Prop.Made ~ Distance.From.Basket.FT , data = SD)
  cov.matrix <- model.matrix.lm(lm1)
  
  beta <- lm1$coefficients #Assume estimates from the model are true
  
  num.pts <- length(SD$Prop.Made)
  
  #Observed Points
  X2 <- cbind(SD$x1 , SD$y1)
  
  #grid - Creating a 50 x 50 grid = 2500 points on the grid
  grid.side <- grid
  
  x.g1 <- seq(min(X2[,1]) , max(X2[,1]) , length.out = grid.side)
  x.g2 <- seq(min(X2[,2]) , max(X2[,2]) , length.out = grid.side)
  
  grid <- expand.grid(x.g1 , x.g2)
  
  mu2 <- c(cov.matrix %*% beta)
  
  d2 <- dist(X2 , upper = T , diag = T) %>% as.matrix()
  
  Omega22 <- (sigma.sq *exp(-d2 * phi)) + tau.sq*diag(num.pts)
  
  resp <- model.frame(lm1)[,1]
  
  dat <- data.frame(x1 = X2[,1] , x2 = X2[,2] , y2 = resp)
  
  Basket.Location <- t(as.matrix(c(0, 0)))
  
  
  Basket <- function(j){
    B <- pointDistance(Basket.Location, grid[j,], lonlat = F, allpairs=TRUE)
    
    return(B)
  }
  
  Distance.From.Basket.Grid <- NULL
  for(i in 1:length(grid[,1])){
    
    Distance.From.Basket.Grid[i] <- Basket(i)
    
    
  }
  
  
  D1 <- as.vector(Distance.From.Basket.Grid)
  D1 <- as.matrix(D1)
  Int <- rep(1 , length(D1))
  
  X1 <- cbind(Int , D1/10)
  
  colnames(X1) <- c("Int" , "Dist.2.Bask")
  
  mu1 <- c(X1 %*% beta)
  
  
  #distance matrix of the grid
  
  d11 <- dist(grid, upper = T , diag  = T) %>% as.matrix()
  
  
  Omega11 <- (sigma.sq*exp(-d11 * phi)) + tau.sq * diag(grid.side^2)
  
  X2 <- as.matrix(X2)
  grid <- as.matrix(grid)
  
  d.big <- dist(rbind(grid, X2) , upper = T, diag = T) %>% as.matrix()
  
  d12 <- d.big[(1:grid.side^2) , ((grid.side^2+1) : (num.pts + grid.side^2))]
  
  
  Omega12 <- sigma.sq * exp(-d12 *phi)
  
  cond.exp <- mu1 + Omega12 %*% solve(Omega22) %*% (resp - mu2)
  
  cond.cov <- (sigma.sq + tau.sq) * diag(grid.side^2) - (Omega12) %*% solve(Omega22) %*% t(Omega12)
  
  exp.df <- data.frame(x1 = grid[,1] , x2 = grid[,2] , cond.exp)
  
  
  try <- ggplot() +
    geom_point(data = exp.df, aes(x=x1, y=x2 , color = cond.exp) ) +
    geom_point(data = dat , aes(x = x1 , y = x2 , color = y2) , size = 4) +
    geom_point(aes(x=0, y=-25), colour="red" , size = 5)
  
  return(try + scale_colour_gradientn(colours = topo.colors(10)))
  
  
}



  
  
