upp <- function(#define coeficients and anatomical traits
  rings = 6,
  intact = 12,
  PLC = 0.5,              #Percentage (in decimal) loss of conductance
  Hcc = 1.83E-02,	       #Hcc Henry's law constant (weighted average for N2 & O2)
  Cg_air = 4.09E+01,	     #mM or mol m-3,Cg,air at standard atmospheric pressure
  Dgw = 5.00E-11,    	   #Dgw, m2 s-1 diffusion coef of gas in wet wood: Justification Sorz J, Hietz P (2006) Gas diffusion through wood: implications for oxygen supply. Trees 20: 34-41
  Dg = 1.00E-09,	         #Dg, m2 s-1 diffusion coef of gas in wet pit membrane
  dv = 4.00E-05,	         #m, dv mean vessel diameter, diameter shared when overlapping
  Vf = 0.05,	             #% (in decimal), fraction of vessel wall surface in common between vessels
  Vlf = 0.1,          	   #% (in decimal) vessel lumen in stem cross section, ax
  Pmf = 0.5,         	   #fraction vessel pit field that is pit membrane
  Pmt = 5.00E-07,    	   #dm m, pit membrane thickness
  Vt = 1.50E-06,	         #m3 volume of external tubing, Vt
  dw = 0.015,	           #m, diam of wood in stem,dw
  dt = 0.05,              #dt, simulation
  total_time = 150        # duration, s
  ) {
  
  Vl = 200000*dv^1.48    #m, mean vessel length
  Vlc = Vl/2             #m, cut open vessel length
  Vv = pi*Vlc*(dv/2)^2   #m3, volume of cut open vessel, Vv
  Vone = pi*Vl*(dv/2)^2  #m3, volume of one vessel, Vone
  Ap = pi*dv*Vl*Pmf*Vf   #Ap m2, total pit membrane surface area
  Nv = Vlf*((dw/dv)^2)   #Nv, Number of vessels= ax(0.015/dv)2
  Nve = Nv*PLC           #Nve, Number of vessels embolized to right
  Vvt = Vv*Nv            #vol all cut open vessels
  Vt_v = Vt+Vvt          #Vt+v, volume of cut open vessel plus tubing in model
  ka = dt*Ap*Dg*Hcc/Pmt    #calculate axial transport rate (ka= dt PM_a Dg Hcc/ddm)
  
  #radial diffusion - cut open vessel
  
  # rmax, m, max radius of UnitPipe
  rmax <- vector()
  rmax[1] <- dw/(2*sqrt(Nv)) #rmax1
  rmax[2] <- dw/(2*sqrt(Nve)) #rmax
  
  dr <- vector()
  rmid <- data.frame() # shell radius, at middle or center of ring
  ki <- data.frame() # radial conductance
  vsh <- data.frame() # shell volume, from edge boundaries
  
  # rings radius for rmax1
  dr[1] <- (rmax[1] - dv/2)/rings
  rmid[1, 1] <- dv/2 + dr[1]/2
  ki[1, 1] <- 2*pi*Dgw*Vlc*dt/log(rmid[1, 1]/(rmid[1, 1] - dr[1]/2))
  
  for (i in seq(2, rings)) {
    rmid[i, 1] <- dr[1] + rmid[i - 1, 1]
    ki[i, 1] <- 2*pi*Dgw*Vlc*dt/log(rmid[i, 1]/(rmid[i, 1] - dr[1]))
  }
  
  for (i in seq(rings)) {
    vsh[i, 1] <- pi*((rmid[i, 1] + dr[1]/2)^2-(rmid[i, 1] - dr[1]/2)^2)*Vlc
  }
  
  # rings radius for rmax
  dr[2] <- (rmax[2] - dv/2)/rings
  rmid[1, 2] <- dv/2 + dr[2]/2
  ki[1, 2] <- 2*pi*Dgw*Vl*dt/log(rmid[1, 2]/(rmid[1, 2] - dr[2]/2))
  
  for (i in seq(2, rings)) {
    rmid[i, 2] <- dr[2] + rmid[i - 1, 2]
    ki[i, 2] <- 2*pi*Dgw*Vl*dt/log(rmid[i, 2]/(rmid[i, 2] - dr[2]))
  }
  
  for (i in seq(rings)) {
    vsh[i, 2] <- Vl*pi*((rmid[i, 2] + dr[2]/2)^2-(rmid[i, 2] - dr[2]/2)^2)
  }
  
  
  # start from atmospheric gas concentration (Cg)
  cg <- rep(Cg_air, intact)
  
  # intial value of radial gas concentration (equlibrium gas/liquid phases = Cg_air*Hcc)
  v <- matrix(Cg_air*Hcc, intact-1, rings)
  
  #initial count (x) and time (dt1)
  x = 1
  dt1 = 0
  
  # create dataframe to colect data
  data <- data.frame()
  data1 <- c(dt1, cg, v[1, ])
  data <-rbind(data, data1) 
  data_names <- vector()
  data_names[1] <- "dt1"
  for (i in seq(intact)) {
    data_names[i+1] <- paste0("Cg", i)
  }
  for (i in seq(rings)) {
    data_names[i+1+intact] <- paste0("v1", i)
  }
  names(data) <- data_names
  
  dn <- vector()
  dnr <- vector()
  #start calculations
  repeat {
    if (dt1 <= 1){
      dp <- Cg_air*(1 - 0.35)*dt          #starting vacuum ~35 kPa in one second
    } else {
      dp = 0
    }
    
    #cut-open vessel
    
    dn[1] <- ka*(cg[2] - cg[1])                                     #axial diffusion, mol
    dnr[1] <- (v[1,1] - Hcc*cg[1])*ki[1,1]                              #radial diffusion, mol
    cg[1] <- cg[1] - dp + (Nve*dn[1] + Nv*dnr[1])/Vt_v                  #gas concentration, mol m^-3
    
    v[1,1] <- v[1,1] + (ki[2,1]*(v[1,2] - v[1,1]) - ki[1,1]*(v[1,1] - Hcc*cg[1]))/vsh[1,1]
    for (i in seq(2, rings - 1)) {
      v[1,i] <- v[1,i]+(ki[i+1,1]*(v[1,i+1]-v[1,i])-ki[i,1]*(v[1,i]-v[1,i-1]))/vsh[i,1]
    }
    v[1,rings] <- v[1,rings]-ki[rings,1]*(v[1,rings]-v[1,rings-1])/vsh[rings,1]
    
    #intact vessel
    for (j in seq(2, intact-1)) {
      dn[j] <- ka*(cg[j+1] - cg[j])
      dnr[j] <- (v[j,1] - Hcc*cg[j])*ki[1,2]
      cg[j] <- cg[j] + (dnr[j] + dn[j] - dn[j-1])/Vone
      
      v[j,1] <- v[j,1] + (ki[2,2]*(v[j,2] - v[j,1]) - ki[1,2]*(v[j,1] - Hcc*cg[j]))/vsh[1,2]
      for (i in seq(2, rings - 1)) {
        v[j,i] <- v[j,i]+(ki[i+1,2]*(v[j,i+1]-v[j,i])-ki[i,2]*(v[j,i]-v[j,i-1]))/vsh[i,2]
      }
      v[j,rings] <- v[j,rings]-ki[rings,2]*(v[j,rings]-v[j,rings-1])/vsh[rings,2]
    }
    #last intact vessel
    cg[intact] <- cg[intact-1]
    
    #save data in dataframe data1
    data1 <- c(dt1, cg, v[1, ])
    data <-rbind(data, data1) 
    
    #count time and repeatition
    dt1 <- dt1 + dt
    x <- x + 1
    if (x == total_time/dt){ # repeat it 3000 times (150 s considering dt = 0.05)
      break
    }
  }
  return(data)
}