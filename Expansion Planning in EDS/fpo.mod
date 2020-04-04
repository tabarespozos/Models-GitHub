set N within {1..4000};
set L within {1..4000} cross N cross N;

set T;				# set of load levels
param dT     {T};
param cese   {T};
param fdem   {T};

set A;                       # set of types of conductors
param RT    {A};             # resistance (ohm/km)
param XT    {A};             # reactance (ohm/km)
param ZT2   {A};             # reactance (ohm/km)
param IT    {A};             # max current (A)
param CFT   {A};             # cost of conductor ($/km)


param D     {L};             # length of the branch (km)
param TO    {L};             # initial conductor type of branch
param Tini  {L};             # optional conductor type of branch
param Xini  {L};             # optional state of branch
param custo {L,A};           # branch maximum current (kA)

set NSE within N;
param Sgo   {NSE};
param Sgp   {NSE};
param Sgc   {NSE};


param vb;
param CLOSS;                 # power losses cost
param C_VARSE;               # substations operation cost
param VMAX;
param VMIN;
param VFK;
param pf;
param K;  # years into one state
param tj; # anual interest rate
param LF; # load factor
param g ; # interest rate for lines
param cl; # cost of Energy (US$/kWh)
param gs; # interest rate for substation
param cs; # operation cost of substation
param alfa; #number of hours in one year

param r3 := 1.73205081;
var volmin;
var xbmax;
param ex;
let ex := 0;

var V2      {T,N}     >= 0;
var Sg2     {T,N}     >= 0;
var Pg      {T,N}     >= 0;
var Qg      {T,N}     >= 0;
param Sd    {N};
param Pd    {N};
param Qd    {N};
param vini  {N,T};
var Ps      {T,L};
var Qs      {T,L};
var I2s     {T,L}     >= 0;
var P       {T,L,A};
var Q       {T,L,A};
var I2      {T,L,A}   >= 0;

var xg      {NSE} binary;

param bmax;
var b       {T,L};
var yp      {L}  binary;
var yn      {L}  binary;


var z       {L,A} binary;

param PWmax >= 1;
set PW := 1..PWmax ;
param mS    {PW};
param DS;
var Dip     {T,L,PW} >= 0;
var Diq     {T,L,PW} >= 0;
var Qp      {T,L}    >= 0;
var Qn      {T,L}    >= 0;
var Pp      {T,L}    >= 0;
var Pn      {T,L}    >= 0;


param mG    {NSE,PW};
param DG    {NSE};
var Dgp     {T,NSE,PW} >= 0;
var Dgq     {T,NSE,PW} >= 0;


subject to SE_capacity1 {n in NSE, t in T}:
  Sg2[t,n] >= Pg[t,n]^2 + Qg[t,n]^2;

# subject to SE_capacity1 {n in NSE, t in T}:
  # Sg2[t,n] = sum {p in PW} mG[n,p] * (Dgp[t,n,p] + Dgq[t,n,p]);

# subject to SE_capacity2 {n in NSE, t in T}:
  # Pg[t,n] = sum {p in PW} Dgp[t,n,p];

# subject to SE_capacity3 {n in NSE, t in T}:
  # Qg[t,n] = sum {p in PW} Dgq[t,n,p];

# subject to SE_capacity4 {n in NSE, p in PW, t in T}:
  # Dgp[t,n,p] <= DG[n];

# subject to SE_capacity5 {n in NSE, p in PW, t in T}:
  # Dgq[t,n,p] <= DG[n];


  
subject to I_soma {(l,m,n) in L, t in T}:
  I2s[t,l,m,n] = sum{a in A} I2[t,l,m,n,a];

subject to P_soma {(l,m,n) in L, t in T}:
  Ps[t,l,m,n] = sum{a in A} P[t,l,m,n,a];  

subject to Q_soma {(l,m,n) in L, t in T}:
  Qs[t,l,m,n] = sum{a in A} Q[t,l,m,n,a];  

subject to CURRENT_FLOW_SQUARE {(l,m,n) in L, t in T}:
  V2[t,n] * I2s[t,l,m,n] >= Ps[t,l,m,n]^2 + Qs[t,l,m,n]^2;

# subject to CURRENT_FLOW_SQUARE {(l,m,n) in L, t in T}:
  # vini[n,t]^2 * I2s[t,l,m,n] = sum {p in PW} mS[p] * fdem[t]^2 * (Dip[t,l,m,n,p] + Diq[t,l,m,n,p]);

# subject to potencia_ativa_soma1 {(l,m,n) in L, t in T}:
  # Ps[t,l,m,n] / fdem[t] = Pp[t,l,m,n] - Pn[t,l,m,n];

# subject to potencia_reativa_soma {(l,m,n) in L, t in T}:
  # Qs[t,l,m,n] / fdem[t] = Qp[t,l,m,n] - Qn[t,l,m,n];

# subject to DELTA_CURRENT_FLOW_re {(l,m,n) in L, t in T}:
  # Pp[t,l,m,n] + Pn[t,l,m,n] = sum {p in PW} Dip[t,l,m,n,p];

# subject to DELTA_CURRENT_FLOW_im {(l,m,n) in L, t in T}:
  # Qp[t,l,m,n] + Qn[t,l,m,n] = sum {p in PW} Diq[t,l,m,n,p];

# subject to MAXIMUN_DELTA_CURRENT_FLOW_re {(l,m,n) in L, p in PW, t in T}:
  # Dip[t,l,m,n,p] <= DS;

# subject to MAXIMUN_DELTA_CURRENT_FLOW_im {(l,m,n) in L, p in PW, t in T}:
  # Diq[t,l,m,n,p] <= DS;

# END OF DISCRETIZATION

# BEGIN POWER FLOW EQUATIONS

subject to ACTIVE_POWER_BALANCE {n in N, t in T}:
  - sum{(l,n,m) in L, a in A} (P[t,l,n,m,a] + I2[t,l,n,m,a] * RT[a] * D[l,n,m]) +
    sum{(l,m,n) in L, a in A} P[t,l,m,n,a] + Pg[t,n] = Pd[n] * fdem[t];

subject to REACTIVE_POWER_BALANCE {n in N, t in T}:
  - sum{(l,n,m) in L, a in A} (Q[t,l,n,m,a] + I2[t,l,n,m,a] * XT[a] * D[l,n,m]) +
    sum{(l,m,n) in L, a in A} Q[t,l,m,n,a] + Qg[t,n] = Qd[n] * fdem[t];

subject to VOLTAGE_DROP {(l,m,n) in L, t in T}:
  V2[t,m] - V2[t,n] = b[t,l,m,n] + sum{a in A} (2 * ( P[t,l,m,n,a] * RT[a] +
       Q[t,l,m,n,a] * XT[a] ) * D[l,m,n] + I2[t,l,m,n,a] * ZT2[a] * D[l,m,n]^2);

# END POWER FLOW EQUATIONS

# OPERATIONAL LIMITS

subject to VM_limits {n in N, t in T}:
  VMIN^2 <= V2[t,n] <= VMAX^2;

subject to IM_max {(l,m,n) in L, a in A, t in T}:
  I2[t,l,m,n,a] <= IT[a]^2 * z[l,m,n,a];

subject to r1:
  sum{(l,m,n) in L} (yp[l,m,n] + yn[l,m,n]) = card(N) - card(NSE);
  
subject to r2_limit_binaria {(l,m,n) in L }:
  yp[l,m,n] + yn[l,m,n] <= 1;

subject to r3_limit_b_neg {(l,m,n) in L, t in T}:
  b[t,l,m,n] >= -bmax * (1 - yp[l,m,n] - yn[l,m,n]);
 
subject to r4_limit_b_pos {(l,m,n) in L, t in T}:
  b[t,l,m,n] <= bmax * (1 - yp[l,m,n] - yn[l,m,n]);
 
subject to r5_potencia_ativa_pos {(l,m,n) in L, a in A, t in T}:
  P[t,l,m,n,a] <= VMAX * IT[a] * yp[l,m,n];

subject to r6_potencia_ativa_neg {(l,m,n) in L, a in A, t in T}:
  P[t,l,m,n,a] >= -VMAX * IT[a] * yn[l,m,n];

subject to r5_potencia_reativa_pos {(l,m,n) in L, a in A, t in T}:
  Q[t,l,m,n,a] <= VMAX * IT[a] * (yp[l,m,n] + yn[l,m,n]);

subject to r6_potencia_reativa_neg {(l,m,n) in L, a in A, t in T}:
  Q[t,l,m,n,a] >= -VMAX * IT[a] * (yp[l,m,n] + yn[l,m,n]);

subject to r11_current_limit {(l,m,n) in L, a in A, t in T}:
  I2[t,l,m,n,a] <= IT[a]^2 * (yp[l,m,n] + yn[l,m,n]);

subject to r11_active_power_limit_z1 {(l,m,n) in L, a in A, t in T}:
  P[t,l,m,n,a] <= VMAX * IT[a] * z[l,m,n,a];

subject to r11_active_power_limit_z2 {(l,m,n) in L, a in A, t in T}:
  P[t,l,m,n,a] >= -VMAX * IT[a] * z[l,m,n,a];

subject to r11_reactive_power_limit_z1 {(l,m,n) in L, a in A, t in T}:
  Q[t,l,m,n,a] <= VMAX * IT[a] * z[l,m,n,a];

subject to r11_reactive_power_limit_z {(l,m,n) in L, a in A, t in T}:
  Q[t,l,m,n,a] >= -VMAX * IT[a] * z[l,m,n,a];

  
subject to r10_reconfig_barras {n in N: Pd[n] > 0}:
  sum{(l,n,m) in L} yn[l,n,m] + sum{(l,m,n) in L} yp[l,m,n] >= 1;
  
subject to UNIQUENESS_type {(l,m,n) in L}:
  sum{a in A} z[l,m,n,a] = yp[l,m,n] + yn[l,m,n];

subject to SE_capacity_total {n in NSE, t in T}:
  Sg2[t,n] <= Sgo[n]^2 + 2 * Sgo[n] * Sgp[n] * xg[n] + Sgp[n]^2 * xg[n];

# conductor's investments
var CF_se = sum{n in NSE} Sgc[n] * xg[n];

var CF_lin = sum{(l,n,m) in L, a in A} z[l,n,m,a] * custo[l,n,m,a] * D[l,n,m];

# power losses
var LOSS = sum{(l,m,n) in L, a in A, t in T} dT[t] * cese[t]* I2[t,l,m,n,a] * RT[a] * D[l,m,n];

var CSE = sum{n in NSE, t in T} Sg2[t,n];

#---- OBJECTIVE FUNCTION.
# Total invesment and operation cost 
minimize COST : CF_se + CF_lin + CLOSS * LOSS + CSE * C_VARSE;
