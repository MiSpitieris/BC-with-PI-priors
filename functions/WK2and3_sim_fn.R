WK2_simulate = function(flow, time, R, C){
  Q = flow
  impedance_WK2 = function(w,R,C){
    Z=R/(1.0+1i*w*R*C)
    return(Z)
  }
  d = time[2]-time[1]
  n = length(time)
  if((n %% 2) == 0) {
    f = c(seq(0,n/2-1, by = 1), seq(-n/2, -1,  by=1))/(d*n)
  } else {
    f = c(seq(0,(n-1)/2, by = 1), seq(-(n-1)/2, -1,  by=1))/(d*n)
  }
  w = 2*pi*f
  Zmod = impedance_WK2(w,R,C)
  Qfft = fft(Q)
  P = Re(fft(Zmod*Qfft, inverse = TRUE))/n
  return(P)
}

WK3_simulate = function(flow, time, R, C, Z = NULL){
  Q = flow
  Z_ao = Z
  impedance_WK3 = function(w,R,C,Z_ao){
    Z=(Z_ao + R + 1i*w*R*C*Z_ao)/(1.0+1i*w*R*C)
    return(Z)
  }
  d = time[2]-time[1]
  n = length(time)
  if((n %% 2) == 0) {
    f = c(seq(0,n/2-1, by = 1), seq(-n/2, -1,  by=1))/(d*n)
  } else {
    f = c(seq(0,(n-1)/2, by = 1), seq(-(n-1)/2, -1,  by=1))/(d*n)
  }
  w = 2*pi*f
  Qfft = fft(Q)
  Zmod = impedance_WK3(w,R,C,Z_ao)
  P = Re(fft(Zmod*Qfft, inverse = TRUE))/n
  return(P)
}
