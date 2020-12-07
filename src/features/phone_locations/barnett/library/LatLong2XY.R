LatLong2XY <-
function(lat_v,lon_v,R=6.371*10^6){
  th0 = min(lon_v)
  th1 = max(lon_v)
  ph0 = min(lat_v)
  ph1 = max(lat_v)
  d1 = 2*pi*R*((ph1-ph0)*2*pi/360)/(2*pi)
  d2 = 2*pi*(R*sin(pi/2-ph1*2*pi/360))*((th1-th0)*2*pi/360)/(2*pi)
  d3 = 2*pi*(R*sin(pi/2-ph0*2*pi/360))*((th1-th0)*2*pi/360)/(2*pi)
  x_v=rep(0,length(lon_v))
  y_v=rep(0,length(lat_v))
  for(i in 1:length(lat_v)){
    w1=(lat_v[i]-ph0)/(ph1-ph0)
    w2=(lon_v[i]-th0)/(th1-th0)
    x_v[i]=w1*abs(d3-d2)/2+w2*(d3*(1-w1)+d2*w1)
    y_v[i]=w1*d1*sin(acos(abs((d3-d2)/(2*d1))))
  }
  return(list("x_v"=x_v,"y_v"=y_v))
}
