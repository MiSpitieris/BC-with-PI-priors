# create data in the [0,1] scale
create_data = function(flow, time, Rtrue, Ctrue, Ztrue=NULL, nc=3, nP=20, nI=17, n_pred=50, Pnoise=3, Inoise=4, seed = 0){
  # 1. simulate data from WK2 or WK3
  if(is.null(Ztrue)){
    Psim = WK2_simulate(flow = flow, time = time, R = Rtrue, C = Ctrue)
  }else{
    Psim = WK3_simulate(flow = flow, time = time, R = Rtrue, C = Ctrue, Z=Ztrue)
  }
  # 2. selected indices according to nP and nI
  nflow = length(flow)
  indP = round(seq(1, nflow, length.out = nP)); indI = round(seq(1, nflow, length.out = nI))
  tP = matrix(d$time[indP], ncol = 1); tI = matrix(d$time[indI], ncol = 1)
  yP_real = Psim[indP]; yI_real = flow[indI]
  # 3. Add noise
  set.seed(seed)
  Pnoise = rnorm(nP*nc, 0, Pnoise)
  Inoise =rnorm(nI*nc, 0, Inoise)
  yP_real = rep(yP_real,each=nc)
  yI_real = rep(yI_real,each=nc)
  
  yP_temp = yP_real + Pnoise
  yI_temp = yI_real + Inoise
  y_temp = c(yP_temp, yI_temp)
  # tranform y in [0,1]
  rPI = range(c(yP_temp, yI_temp))
  yP = (yP_temp - rPI[1])/diff(rPI)
  yI = (yI_temp - rPI[1])/diff(rPI)
  
  ind_pred = round(seq(1,101, length.out = n_pred))
  data_noisy_pred = list(nP = nc*nP, nI = nc*nI, tP = matrix(rep(tP,each=nc),ncol=1)
                         , tI = matrix(rep(tI,each=nc),ncol=1), yP = yP, yI = yI
                         , tP_pred = matrix(d$time[ind_pred],ncol = 1), nP_pred = n_pred
                         , tI_pred = matrix(d$time[ind_pred],ncol = 1), nI_pred = n_pred
                         )
  
  # true data for plotting
  data_mod_true = data.frame(time = time, I = flow, P = Psim)
  return(list(data_noisy_pred = data_noisy_pred, data_mod_true =data_mod_true, y=y_temp))
}

# transform to the real scale
transform_post = function(y, fit){
  post_df = as.data.frame(extract(fit))
  ry = diff(range(y))
  my = min(y)
  cn = colnames(post_df)
  
  post_df[,c("R", "C")] = post_df[,c("R", "C")]*2.5 +0.5
  indPI = c(grep("_P", cn), grep("_I", cn))
  post_df[,indPI] = post_df[,indPI]*ry + my
  # ind_delta = grep("delta",cn)
  ind_delta = c(grep("sigma",cn), grep("mu_wk2",cn))
  md = apply(post_df[,ind_delta],2,min)
  post_df[,ind_delta] = post_df[,ind_delta]*ry + md
  return(post_df)
}