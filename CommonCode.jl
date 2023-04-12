################################################################
################################################################
################################################################
################################################################
# Copyright (C) 2023 Robert Gowers and Magnus Richardson
#
# This program is free software; you can redistribute it and/or 
# modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation; 
# either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program; if not, see <https://www.gnu.org/licenses>.
################################################################
################################################################
################################################################
################################################################

################################################################
################################################################
################################################################
################################################################
# Common parameters and functions
################################################################
################################################################
################################################################
################################################################

################################################################
################################################################
# Common parameters
################################################################
################################################################

# Useful constants
const K=1000.0
const r2d=360.0/(2.0*pi)

# Voltage parameters
const EL=-60.0 #[mV] resting potential in absence of input
const Ee=0.0 #[mV]
const Ei=-80.0 #[mV] reversal potentials
const Vre=-60.0
const Vth=-50.0 #[mV] threshold reset for spikes

# Time constants
const tauL=40.0 #[ms]
const taue=3.0 #[ms]
const taui=10.0 #[ms]

# Mean synaptic conductances - single synpase conductances typically in range 0.1-1nS 
const game=0.1
const gami=0.1

###############################################################
# Other dendrite parameters
###############################################################

c=0.01 #[pF/um^2] from C=1.0 uF/cm^2 (0.9 in Gentet et al Biophy. J., 79 (2000) 314–320)
a=0.5 #[um] 1um reasonable for a relatively thick apical dendrite
ra = 0.02 #[GOhm um]

const lamL=sqrt(a*tauL/(2*ra*c)); # [um]
const lame=2*game*taue/(2*pi*a*c) # [um] NB exponentially distributed synaptic conductances
const lami=2*gami*taui/(2*pi*a*c) # [um] NB exponentially distributed synaptic conductances

const D=lamL^2/tauL #[um^2/ms] this is numerically large!

const kappa=0.6 # Conductance state of the neuron ratio tau0/tauL
const tau0=kappa*tauL #[ms] effective membrane time constant. 

###############################################################
# Point-neuron parameters for matched variance
###############################################################

const C=2*pi*a*lamL*c  #[pF] 
const kappae=(lame/lamL)*((tau0+taue)/tau0)*(1-sqrt(taue/(tau0+taue)))/(2*sqrt(kappa))
const kappai=(lami/lamL)*((tau0+taui)/tau0)*(1-sqrt(taui/(tau0+taui)))/(2*sqrt(kappa))

###############################################################
# A quantity related to the state of the neuro during firing
###############################################################

const rstarHz=5.0
const rstar=5.0/1000.0
const L=2000.00

###############################################################
# A quantity related to the state of the neuro during firing
###############################################################

println("\n------------------------------------------")
println("EL=$EL \tEe=$Ee  \tEi=$Ei")
println("tauL=$tauL \ttaue=$taue \ttaui=$taui")

lamLS,lameS,lamiS=(round.((lamL,lame,lami),digits=0))
println("lamL=$lamLS \tlame=$lameS \tlami=$lamiS")

println("\t\tgame=$game \tgami=$gami")
println("\nkappa=$kappa tau0=$tau0")
println("D=$D")
println("Vth=$Vth Vre=$Vre")
println("L=$L")
println("------------------------------------------\n")

################################################################
################################################################
# General functions
################################################################
################################################################

function FixIt(X) # tidies up NaNs and Infs 
    sx=findall(isinf.(X))
    X[sx].=1
    sx=findall(isnan.(X))
    X[sx].=1
    return X 
end

function Smooth(x,c,eps)
    n=length(x)
    for k=1:c
        x=(1-2*eps)*x .+ eps*([x[1];x[1:end-1]] .+ [x[2:end];x[end]])
    end
    return x
end

function MyCosFitSim(t,x,w)
    
    mx=mean(x)
    T=t[end]
    dt=t[2]-t[1]
    X1=(2/T)*sum(dt*(x.-mx).*cos.(t*w))
    X2=(2/T)*sum(dt*(x.-mx).*sin.(t*w))
    asim=sqrt(X1^2 + X2^2)
    psim=angle(X1 -1im*X2)
    
    return mx,asim,psim
    
end

################################################################
################################################################
# Steady-state theory
################################################################
################################################################

function PointSteadyTheory(ae0::Float64,ai0::Float64)
    
    # First-order moment
    tau0=1/(1.0/tauL+ae0+ai0)
    V0=(EL/tauL + ae0*Ee + ai0*Ei)*tau0
    
    # Useful quantities for second-order moments
    Ee0=Ee-V0 
    Ei0=Ei-V0
    vth=Vth-V0
    
    # Second-order moments
    vhe0=Ee0*(ae0*kappae/2taue)/(1/taue+1/tau0)
    vhi0=Ei0*(ai0*kappai/2taui)/(1/taui+1/tau0)
    
    vve0=Ee0^2*(ae0*kappae/2taue)*tau0/(1/tau0+1/taue)
    vvi0=Ei0^2*(ai0*kappai/2taui)*tau0/(1/tau0+1/taui)
    vv0=vve0 + vvi0

    dvdve0=Ee0^2*(ae0*kappae/(2*taue^2))/(1/tau0+1/taue)
    dvdvi0=Ei0^2*(ai0*kappai/(2*taui^2))/(1/tau0+1/taui)
    dvdv0=dvdve0+dvdvi0
    
    # Upcrossing rate    
    ru0=(1/(2pi))*(sqrt(dvdv0)/sqrt(vv0))*exp(-vth^2/(2*vv0))
    
    return V0,vhe0,vhi0,vv0,dvdv0,ru0

end

function DendrSteadyTheory(ae0::Float64,ai0::Float64)
    
    # First-order moment
    tau0=1/(1.0/tauL+ae0+ai0)
    V0=(EL/tauL + ae0*Ee + ai0*Ei)*tau0
    
    # Useful quantities for second-order moments
    Ee0=Ee-V0 
    Ei0=Ei-V0
    vth=Vth-V0
    lam0=lamL*sqrt(tau0/tauL)
    
    # Second-order moments
    vhe0=(Ee0/(4*taue))*(tau0*ae0)*(lame/lam0)*sqrt(taue/(tau0+taue))
    vhi0=(Ei0/(4*taui))*(tau0*ai0)*(lami/lam0)*sqrt(taui/(tau0+taui))

    vve0=(Ee0^2/4)*(tau0*ae0)*(lame/lam0)*(1 - sqrt(taue/(tau0+taue)))
    vvi0=(Ei0^2/4)*(tau0*ai0)*(lami/lam0)*(1 - sqrt(taui/(tau0+taui)))
    vv0=vve0+vvi0
    
    # variance of dvdt
    dvdve0=(Ee0^2/(2*taue)^2)*(tau0*ae0)*(lame/lam0)*sqrt(taue/(tau0+taue))
    dvdvi0=(Ei0^2/(2*taui)^2)*(tau0*ai0)*(lami/lam0)*sqrt(taui/(tau0+taui))
    dvdv0=dvdve0+dvdvi0
    
    # Upcrossing rate
    ru0=(1/(2pi))*(sqrt(dvdv0)/sqrt(vv0))*exp(-vth^2/(2*vv0))
     
    return V0,vhe0,vhi0,vv0,dvdv0,ru0
    
end;

function DendrCovarianceSteadyTheory(x,ae0::Float64,ai0::Float64)

    # Get the x=0 covariance results
    V0,vhe0,vhi0,vv0,dvdv0,ru0=DendrSteadyTheory(ae0,ai0)

    # Useful quantities for second-order moments
    Ee0=Ee-V0 
    Ei0=Ei-V0
    vth=Vth-V0
    tau0=1/(1.0/tauL+ae0+ai0)
    lam0=lamL*sqrt(tau0/tauL)
    
    ke0=sqrt((taue+tau0)/taue)/lam0
    ki0=sqrt((taui+tau0)/taui)/lam0
    k0=1/lam0

    # covariance between synaptic inputs and voltage
    vhe0x=vhe0*exp.(-abs.(x)*ke0)
    vhi0x=vhi0*exp.(-abs.(x)*ki0)

    # voltage autocovariance
    ce=(Ee0^2)*(tau0*ae0/4.0)*(lame/lam0)
    ci=(Ei0^2)*(tau0*ai0/4.0)*(lami/lam0)
    vv0xe=ce*( exp.(-abs.(x)*k0) .- (sqrt(taue/(tau0+taue)))*exp.(-abs.(x)*ke0))
    vv0xi=ci*( exp.(-abs.(x)*k0) .- (sqrt(taui/(tau0+taui)))*exp.(-abs.(x)*ki0))
    vv0x=vv0xe.+vv0xi

    # rate of change of voltage covariance
    dvdv0x=Ee0*vhe0x./taue .+ Ei0*vhi0x./taui

    return vhe0x,vhi0x,vv0x,dvdv0x
    
end

################################################################
################################################################
# Theory for general dynamics
################################################################
################################################################

function PointDynamicTheory(t::Vector{Float64},aet::Vector{Float64},ait::Vector{Float64})
    
    dt=t[2]-t[1]; 
    nt=length(t)
    
    mHL=1/tauL
    
    ##########################################################
    # Deterministic quantities
    ##########################################################
    mHe=zeros(nt)
    mHi=zeros(nt)
    mU=EL*ones(nt); # at rest at the beginning
    
    # the deterministic evolution
    for k=1:nt-1
        mHe[k+1]=mHe[k] + (dt/taue)*(aet[k] - mHe[k])
        mHi[k+1]=mHi[k] + (dt/taui)*(ait[k] - mHi[k])
        mU[k+1]=mU[k]   + dt*(mHL*(EL-mU[k]) + mHe[k]*(Ee-mU[k]) + mHi[k]*(Ei-mU[k]) )
    end

    # useful quantities
    Ee0=(Ee .-mU)
    Ei0=(Ei .-mU)
    H=mHL .+ mHe .+ mHi;
    mdUdt=(mHL*(EL .-mU) .+ mHe.*(Ee .-mU) .+ mHi.*(Ei .-mU))

    ##########################################################
    # Second-order evolution
    ##########################################################

    hehe,hihi,vhe,vhi,vv,dvdv=(zeros(nt) for k=1:6)
    
    # second-order evolution
    for k=1:nt-1
        
        hehe[k+1]=hehe[k] + (2*dt/taue)*(kappae*aet[k]/(2*taue) - hehe[k])
        hihi[k+1]=hihi[k] + (2*dt/taui)*(kappai*ait[k]/(2*taui) - hihi[k])
        
        vhe[k+1]=vhe[k] + dt*(Ee0[k]*hehe[k] - (H[k] + 1/taue)*vhe[k])
        vhi[k+1]=vhi[k] + dt*(Ei0[k]*hihi[k] - (H[k] + 1/taui)*vhi[k])
        
        vv[k+1]=vv[k] + (2*dt)*(Ee0[k]*vhe[k] + Ei0[k]*vhi[k] - H[k]*vv[k])
        
    end
    
    dvhe=Ee0.*hehe .- H.*vhe
    dvhi=Ei0.*hihi .- H.*vhi
    vdv=Ee0.*vhe .+ Ei0.*vhi .- H.*vv
    dvdv=Ee0.*dvhe .+ Ei0.*dvhi .- H.*vdv

    ##########################################################
    # Upcrossing rate
    ##########################################################
    
    vth=Vth.-mU
    
    kappa=vdv./vv; 
    kappa=FixIt(kappa)
    
    s2=dvdv .-(vdv.^2)./vv
    s2=FixIt(s2)
    s2[findall(s2.<0)].=1
    
    bbeta=(mdUdt .+ kappa.*(Vth .-mU))./sqrt.(2*s2)
    bbeta=FixIt(bbeta);
    
    Q1=(1/(2*pi))*sqrt.(s2./vv).*exp.(-(vth.^2)./(2*vv))
    Q2=(exp.(-bbeta.^2) .+ (bbeta.*sqrt(pi)).*(1 .+ erf.(bbeta))   )
    rU=Q1.*Q2;
    rU[isnan.(rU)].=0;
    
    return mU,mHe,mHi,mdUdt,vhe,vhi,vv,vdv,dvdv,rU
    
end

function DendrDynamicTheory(x::Vector{Float64},t::Vector{Float64},aet::Vector{Float64},ait::Vector{Float64})
    
    nt=length(t); dt=t[2]-t[1]
    nx=length(x); dx=x[2]-x[1]

    # deterministic quantities
    mHe=zeros(nt)
    mHi=zeros(nt)
    mHL=1/tauL
    mV=EL*ones(nt); # at rest at the beginning
    
    # the deterministic evolution
    for k=1:nt-1
        mHe[k+1]=mHe[k] + (dt/taue)*(aet[k] - mHe[k])
        mHi[k+1]=mHi[k] + (dt/taui)*(ait[k] - mHi[k])
        mV[k+1] =mV[k]  + dt*(mHL*(EL -mV[k]) + mHe[k]*(Ee-mV[k]) + mHi[k]*(Ei-mV[k]))
    end
    mdVdt=mHL*(EL .-mV) .+ mHe.*(Ee .-mV) .+ mHi.*(Ei .-mV) # deterministic dVdt
    
    # set up the noisy part
    fe=zeros(nt)
    fi=zeros(nt)
    vhe=zeros(nx,nt)
    vhi=zeros(nx,nt)
    vv=zeros(nx,nt)
    vdv=zeros(nx,nt)
    dvdv=zeros(nx,nt)
    
    # useful quantities
    EEe,EEi=(Ee .-mV),(Ei .-mV)
    H=(mHL .+ mHe .+ mHi)
    
    dd0=zeros(nx); dd0[1]=1.0./dx # Dirac delta at 0
    
    # second-order evolution. 
    # Variables on half sites, gradients on site because gradients at 0 are known
    
    for k=1:nt-1
        
        fe[k+1]=fe[k] + dt*(aet[k]*lame/taue^2 - fe[k]*2/taue)
        fi[k+1]=fi[k] + dt*(ait[k]*lami/taui^2 - fi[k]*2/taui)
        
        # hev and hiv

        Dd2vhedx2=D*diff(diff([vhe[1,k];vhe[:,k];vhe[end,k]]))/dx^2
        Dd2vhidx2=D*diff(diff([vhi[1,k];vhi[:,k];vhi[end,k]]))/dx^2
        Dd2vhedx2[1]=( D*(vhe[2,k]-vhe[1,k])/dx + EEe[k]*fe[k]/2 )/dx
        Dd2vhidx2[1]=( D*(vhi[2,k]-vhi[1,k])/dx + EEi[k]*fi[k]/2 )/dx
        vhe[:,k+1]=vhe[:,k] .+ dt*(Dd2vhedx2 .- (H[k].+1/taue).*vhe[:,k]) 
        vhi[:,k+1]=vhi[:,k] .+ dt*(Dd2vhidx2 .- (H[k].+1/taui).*vhi[:,k]) 
        
        # v2 & vdvdt
        Dd2vvdx2=D*diff(diff([vv[1,k];vv[:,k];vv[end,k]]))/dx^2
        Dd2vvdx2[1]=D*(vv[2,k]-vv[1,k])/ dx^2
        vv[:,k+1]=vv[:,k] .+ 2*dt*( EEe[k]*vhe[:,k] .+EEi[k]*vhi[:,k] .-H[k]*vv[:,k] .+Dd2vvdx2 ) 
        vdv[:,k]=EEe[k]*vhe[:,k] .+EEi[k]*vhi[:,k] .-H[k]*vv[:,k] .+Dd2vvdx2
        
        # stuff for dvdt2
        dvhe=(vhe[:,k+1]-vhe[:,k])/dt .+vhe[:,k]/taue
        dvhi=(vhi[:,k+1]-vhi[:,k])/dt .+vhi[:,k]/taui
        Dd2vdvdx2=D*diff(diff([vdv[1,k];vdv[:,k];vdv[end,k]]))/dx^2
        Dd2vdvdx2[1]=D*(vdv[2,k]-vdv[1,k])/ dx^2
        dvdv[:,k]=EEe[k]*dvhe .+ EEi[k]*dvhi - H[k]*vdv[:,k] .+ Dd2vdvdx2
    
    end
    
    ##########################################################
    # Upcrossing rate 
    ##########################################################

    # make correction to get quanatities at the origin
    vvo=(9vv[1,:] .- vv[2,:])/8 # using zero gradient at origin
    vdvo=(9vdv[1,:] .- vdv[2,:])/8 # using zero gradient at origin
    dvdvo=(3dvdv[1,:] .- dvdv[2,:])/2 # linear extrapolation

    vth=Vth.-mV
    kappa=vdvo./vvo; 
    kappa=FixIt(kappa);
    s2=dvdvo .-(vdvo.^2)./vvo
    s2=FixIt(s2); s2[findall(s2.<0)].=1
    Q1=(1/(2*pi))*sqrt.(s2./vvo).*exp.(-(vth.^2)./(2*vvo))
    Q1[findall(isnan.(Q1))].=0
    bbeta=(mdVdt .+ kappa.*(Vth .-mV))./sqrt.(2*s2)
    bbeta=FixIt(bbeta);
    Q2=(exp.(-bbeta.^2) .+ (bbeta.*sqrt(pi)).*(1 .+ erf.(bbeta))   )
    
    ru=Q1.*Q2;
    
    return mV,mHe,mHi,mdVdt,vhe,vhi,vv,vdv,dvdv,ru
    
end

################################################################
################################################################
# Theory for oscillations
################################################################
################################################################

function PointOscTheory(ae0::Float64,ai0::Float64,fkHz::Vector{Float64})
    
    # unpack
    w=2*pi*fkHz
    aL=1/tauL
    
    ##################################################
    # Steady-state results
    ##################################################
    V0,vhe0,vhi0,vv0,dvdv0,ru0=PointSteadyTheory(ae0,ai0)
    tau0=1/(aL+ae0+ai0)
    Ee0=(Ee-V0)
    Ei0=(Ei-V0)
    H0=1/tau0
    
    # Other steady-state deterministic
    hehe0=kappae*ae0/2taue
    hihi0=kappai*ai0/2taui
    dvhe0=Ee0*hehe0 .- vhe0/tau0
    dvhi0=Ei0*hihi0 .- vhi0/tau0
    
    # first-order moments for Re modulation
    V1=Ee0./((1.0 .+1im*w*taue).*(1im*w .+H0))
    dVdt1=1im*w.*V1
    He1=1.0./(1.0 .+1im*w*taue)
    Hi1=0  # No inhib osc here
    H1=He1
    
    # second-order moments
    hehe1=(kappae/(2taue))./(1.0 .+1im*w*taue/2)
    hihi1=0.0 # No inhib osc here
    
    vhe1=(Ee0*hehe1 .- V1.*hehe0 .- H1.*vhe0)./(1im*w .+ 1/tau0 .+ 1/taue)
    vhi1=-(V1.*hihi0 .+ H1.*vhi0)./(1im*w .+ 1/tau0 .+ 1/taui)
    
    vv1=((Ee0*vhe1 .+ Ei0*vhi1) .- V1.*(vhe0 .+ vhi0) .- H1.*vv0)./(1im*w/2 .+ 1/tau0)
    vdv1=vv1.*1im.*w/2
    
    dvhe1=(Ee0*hehe1 .-V1.*hehe0) .- H1.*vhe0 .- vhe1/tau0
    dvhi1=-V1.*hihi0 .- H1.*vhi0 .- vhi1/tau0
    dvdv1=Ee0*dvhe1 .+ Ei0*dvhi1 .- V1.*(dvhe0 .+ dvhi0) .-H0*vdv1
    
    # upcrossing rate
    
    theta0=Vth.-V0
    K1=0.5*dvdv1./dvdv0
    K2=-0.5*vv1./vv0
    K3=sqrt(pi/(2*dvdv0))*dVdt1
    K4=sqrt(pi/(2*dvdv0))*theta0*(vdv1/vv0)
    K5=(theta0^2/(2vv0))*(2*V1/theta0)
    K6=(theta0^2/(2vv0))*(vv1/vv0)
    ru1=ru0.*(K1.+K2.+K3.+K4.+K5.+K6)

    return V0,vv0,dvdv0,V1,dVdt1,vhe1,vhi1,vv1,vdv1,dvdv1,ru1

end

function PointOscAsym(ae0::Float64,ai0::Float64,fkHz::Vector{Float64})
    
    w=2*pi*fkHz # for the high-freq case only
    nw=length(w)
    
    # Get steady-state objects
    V0,vhe0,vhi0,vv0,dvdv0,ru0=PointSteadyTheory(ae0,ai0)
    HL,He0,Hi0=1/tauL,ae0,ai0
    H0=HL+He0+Hi0
    V0=(HL*EL+He0*Ee+Hi0*Ei)/H0
    tau0=1/H0
    Ee0=(Ee-V0)
    Ei0=(Ei-V0)
    
    #########################################################
    # Low-frequency limit
    #########################################################
    
    # V1
    V1L=Ee0*tau0*ones(nw)
    
    # vhe and vhi
    Qe=(1 - ae0*tau0 - ae0*tau0*taue/(tau0+taue))
    vhe1L=(Ee0*kappae/2)*(tau0/(tau0+taue))*Qe*ones(nw)
    Qi=Ee0 + Ei0*taui/(tau0+taui)
    vhi1L=-(kappai/2)*(tau0*ai0)*(tau0/(tau0+taui))*Qi*ones(nw)
    
    # vv
    Qe=1-3*ae0*tau0-ae0*taue*tau0/(tau0+taue)
    vv1Le=tau0*(kappae*Ee0^2/2)*(tau0/(tau0+taue))*Qe*ones(nw)
    Qi=2Ee0*Ei0*tau0*ai0+Ei0^2*tau0*ai0+Ei0^2*ai0*taui*tau0/(taui+tau0)
    vv1Li=-tau0*(kappai/2)*(tau0/(tau0+taui))*Qi*ones(nw)
    vv1L=vv1Le.+vv1Li
    
    # dvdv
    Qe=(1 - 2*tau0*ae0 - ae0*tau0*taue/(tau0+taue))
    dvdv1Le=(kappae*Ee0^2/(2taue))*(tau0/(tau0+taue))*Qe*ones(nw)
    Qi=(2*Ei0*Ee0*tau0*ai0+Ei0^2*ai0*taui*tau0/(taui+tau0))
    dvdv1Li=-(kappai/(2taui))*(tau0/(tau0+taui))*Qi*ones(nw)
    dvdv1L=dvdv1Le .+ dvdv1Li
    
    # upcrossing rate
    theta0=Vth.-V0
    Qru1=theta0*V1L/vv0
    Qru2=0.5*(vv1L/vv0)*(theta0^2/vv0 - 1.0)
    Qru3=0.5*dvdv1L/dvdv0
    rU1L=ru0*(Qru1.+Qru2.+Qru3)
    
    #########################################################
    # High-frequency limit
    #########################################################
    
    # V1 asymptotics
    dVdt1H=Ee0./(1im*w*taue)
    H1=1.0./(1im*w*taue)
    
    # hehe1 & H1
    hehe1H=kappae./(1im*w*taue*taue)
    
    # vdv
    vdv1H=-vv0./(1im*w*taue)
    
    # vhe and vhi
    vhe1H=(hehe1H*Ee0 .-H1.*vhe0)./(1im*w)
    vhi1H=.-H1.*vhi0./(1im*w)
    
    # dvdv
    dvhe1H=Ee0*kappae./(1im*w*taue*taue) .- vhe0./(1im*w*taue)
    dvhi1H= - vhi0./(1im*w*taue)
    dvdv1H=Ee0*dvhe1H .+ Ei0*dvhi1H .-vdv1H/tau0
    
    # the high-freq asymptote divided by Re1
    theta0=Vth-V0
    A1=sqrt(pi/(2*dvdv0))*dVdt1H
    A2=sqrt(pi/(2*dvdv0))*theta0*vdv1H/vv0
    A3=(1/2)*dvdv1H/dvdv0
    rU1H=ru0*(A1 .+A2 .+A3)
    
    return vhe1L,vhi1L,vv1L,dvdv1L,rU1L,dVdt1H,vhe1H,vhi1H,vdv1H,dvdv1H,rU1H
    
end

function DendrOscTheoryALT(ae0::Float64,ai0::Float64,fkHz::Vector{Float64})
    
    w=2*pi*fkHz
    aL=1/tauL

    ##################################################
    # Steady-state results
    ##################################################
    V0,vhe0,vhi0,vv0,dvdv0,ru0=DendrSteadyTheory(ae0,ai0)
    tau0=1/(aL+ae0+ai0)
    lam0=lamL*sqrt(tau0/tauL)
    Ee0=(Ee-V0)
    Ei0=(Ei-V0)

    ##################################################
    # Modulation, deterministic
    # All quantities stripped of alpha_1
    # Only the excitatory case included. 
    ##################################################
    He1=1.0./(1.0 .+1im*w*taue)
    Hi1=0.0
    H1=He1
    V1=Ee0.*He1./(1im*w .+ 1/tau0)
    dVdt1=1im*w.*V1
    Ee1=-V1
    Ei1=-V1

    ##################################################
    # Modulation, deterministic
    # All quantities stripped of alpha_1
    # Only the excitatory case included. 
    ##################################################
    
    qe0=sqrt.((1im*w .+ 1/tau0 .+1/taue)/D)
    ke0=sqrt.((1/tau0 .+1/taue)/D)
    qi0=sqrt.((1im*w .+ 1/tau0 .+1/taui)/D)
    ki0=sqrt.((1/tau0 .+1/taui)/D)
    q0=sqrt.((1im*w/2 .+ 1/tau0)/D)
    k0=sqrt.((1/tau0)/D)
  
    ##################################################
    # <vhe>_1 
    ##################################################
    
    A1=Ee0*tau0./(1 .+ 1im*w*taue/2)
    A2=Ee1*ae0*tau0 .+ H1*Ee0*ae0*tau0./(1im.*w)
    A3=vhe0*H1./(1im.*w)
    vhe1=(1.0 ./(4.0 .*taue .*qe0 .*lam0)).*(lame/lam0).*(A1 .+A2) .- A3

    ##################################################
    # <vhi>_1 
    ##################################################

    A1=0
    A2=Ei1*ai0*tau0 .+ H1*Ei0*ai0*tau0./(1im*w)
    A3=vhi0*H1./(1im*w)
    vhi1=(1.0 ./(4*taui*qi0*lam0)).*(lami/lam0).*(A1 .+A2) .- A3 

    ##################################################
    # vv1 
    ##################################################

    a=-2*H1 ./(1im*w)
    be=-Ee0*taue ./(1 .+ 1im*w*taue/2)
    bi=-Ei0*taui ./(1 .+ 1im*w*taui/2)
    ce=taue*(a*Ee0 - Ee1 - be*H1) ./ (1 .- 1im*w*taue/2)
    ci=taui*(a*Ei0 - Ei1 - bi*H1) ./ (1 .- 1im*w*taui/2)
    
    Ae=be.*( Ee0*tau0./(1 .+ 1im*w*taue/2) .+ Ee1 .* ae0 * tau0) .+ ce.*Ee0*ae0.*tau0
    psie=(1/(4*taue))*(lame/lam0)*(1 ./(lam0.*qv)).*Ae

    Ai=bi.*( Ei0*tau0./(1 .+ 1im*w*taui/2) .+ Ei1 .* ai0 * tau0) .+ ci.*Ei0*ai0.*tau0
    psii=(1/(4*taui))*(lami/lam0)*(1 ./(lam0.*qv)).*Ai

    vve1=ae.*vv0 .+ be.*vhe1 .+ ce.*vhe0 .+ psie
    vvi1=ai.*vv0 .+ bi.*vhi1 .+ ci.*vhi0 .+ psii
    vv1=vve1 .+ vvi1

    vdv1=vv1
    dvdv1=vv1
    ru1=vv1

    return V0,vv0,dvdv0,V1,dVdt1,vhe1,vhi1,vv1,vdv1,dvdv1,ru1

end

function DendrOscTheory(ae0::Float64,ai0::Float64,fkHz::Vector{Float64})
    
    w=2*pi*fkHz
    aL=1/tauL

    ##################################################
    # Steady-state results
    ##################################################
    V0,vhe0,vhi0,vv0,dvdv0,ru0=DendrSteadyTheory(ae0,ai0)
    tau0=1/(aL+ae0+ai0)
    lam0=lamL*sqrt(tau0/tauL)
    Ee0=(Ee-V0)
    Ei0=(Ei-V0)

    ##################################################
    # Modulation, deterministic
    # All quantities stripped of alpha_1
    # Only the excitatory case included. 
    ##################################################
    He1=1.0./(1.0 .+1im*w*taue)
    Hi1=0.0
    H1=He1.+Hi1
    V1=Ee0.*He1./(1im*w .+ 1/tau0)
    dVdt1=1im*w.*V1
    H1=He1
    
    ##################################################
    # Some useful functions for second-order moments
    ##################################################
    Q(k1,k2)=1.0./((lam0^2).*(k1.^2 .- k2.^2))
    X1(k1)=1.0./(2.0.*lam0.*k1)
    X2(k1,k2)=Q(k2,k1).*X1(k1) .+Q(k1,k2).*X1(k2)
    X3(k1,k2,k3)=Q(k2,k1).*Q(k3,k1).*X1(k1) .+ Q(k1,k2).*Q(k3,k2).*X1(k2) .+ Q(k1,k3).*Q(k2,k3).*X1(k3)
   
    Y1(k1)=-lam0*k1.*(1im.*w*tau0)/(4(tau0^2))
    Y2(k1,k2)=Q(k2,k1).*Y1(k1) .+Q(k1,k2).*Y1(k2)
    Y3(k1,k2,k3)=Q(k2,k1).*Q(k3,k1).*Y1(k1) .+ Q(k1,k2).*Q(k3,k2).*Y1(k2) .+ Q(k1,k3).*Q(k2,k3).*Y1(k3)
    
    ##################################################
    # Modulation, deterministic
    # All quantities stripped of alpha_1
    # Only the excitatory case included. 
    ##################################################
    
    qe0=sqrt.((1im*w .+ 1/tau0 .+1/taue)/D)
    ke0=sqrt.((1/tau0 .+1/taue)/D)
    qi0=sqrt.((1im*w .+ 1/tau0 .+1/taui)/D)
    ki0=sqrt.((1/tau0 .+1/taui)/D)
    q0=sqrt.((1im*w/2 .+ 1/tau0)/D)
    k0=sqrt.((1/tau0)/D)
    
    ##################################################
    # <hev>_1 
    ##################################################
    A1=tau0*(Ee0/(2taue))*(lame/lam0)*X1(qe0)./(1 .+1im*w*taue/2)
    A2=tau0*(Ee0/(2taue))*(lame/lam0)*(ae0*tau0)*X1(qe0)./((1 .+1im*w*taue).*(1 .+1im*w*tau0));
    A3=tau0*(Ee0/(2taue))*(lame/lam0)*(ae0*tau0)*X2(qe0,ke0)./(1 .+1im*w*taue);
    vhe1=A1.-A2.-A3
    
    ##################################################
    # <hiv>_1
    ##################################################
    B2=tau0*(Ee0/(2taui))*(lami/lam0)*(ai0*tau0)*X1(qi0)./((1 .+1im*w*taue).*(1 .+1im*w*tau0));
    B3=tau0*(Ei0/(2taui))*(lami/lam0)*(ai0*tau0)*X2(qi0,ki0)./(1 .+1im*w*taue);
    vhi1= .-B2 .-B3

    ##################################################
    # <v2>_1
    ##################################################
    c1=tau0*(Ee0^2/2)*(tau0/taue)*(lame/lam0)
    c2=tau0*(Ee0^2/2)*(tau0/taue)*(lame/lam0)*(ae0*tau0)
    c3=tau0*(Ee0*Ei0/2)*(tau0/taui)*(lami/lam0)*(ai0*tau0)
    c4=tau0*(Ei0^2/2)*(tau0/taui)*(lami/lam0)*(ai0*tau0)
    
    C1=c1*X2(qe0,q0)./(1 .+ 1im*w*taue/2);
    C2=c2*X2(qe0,q0)./((1 .+ 1im*w*taue).*(1 .+ 1im*w*tau0));
    C3=c2*X3(qe0,ke0,q0)./((1 .+ 1im*w*taue))
    C4=c3*X2(qi0,q0)./((1 .+ 1im*w*taue).*(1 .+ 1im*w*tau0))
    C5=c4*X3(qi0,ki0,q0)./(1 .+ 1im*w*taue)
    C6=c2*X2(ke0,q0)./((1 .+ 1im*w*taue).*(1 .+ 1im*w*tau0))
    C7=c3*X2(ki0,q0)./((1 .+ 1im*w*taue).*(1 .+ 1im*w*tau0))
    C8=c2*X3(q0,k0,ke0)./(1 .+ 1im*w*taue)
    C9=c4*X3(q0,k0,ki0)./(1 .+ 1im*w*taue)
    
    vv1=C1.-C2.-C3.-C4.-C5.-C6.-C7.-C8.-C9
    
    ##################################################
    # <vdv>_1. Straightforward
    ##################################################
    vdv1=1im*w.*vv1/2
    
    ##################################################
    # <dvdv>_1
    ##################################################
    D1=Ee0*vhe1.*(1/taue .+1im.*w).+Ei0*vhi1.*(1/taui .+1im.*w)
    D2=V1.*(vhe0.*(1/taue .+1im.*w) .+ vhi0.*(1/taui .+1im.*w))
    D3=vdv1/tau0
    #NB bit below is the same as for <v^2>_1 but with X->Y to give Dk^2<vdv>_1
    E1=c1*Y2(qe0,q0)./(1 .+ 1im*w*taue/2);
    E2=c2*Y2(qe0,q0)./((1 .+ 1im*w*taue).*(1 .+ 1im*w*tau0));
    E3=c2*Y3(qe0,ke0,q0)./((1 .+ 1im*w*taue))
    E4=c3*Y2(qi0,q0)./((1 .+ 1im*w*taue).*(1 .+ 1im*w*tau0))
    E5=c4*Y3(qi0,ki0,q0)./(1 .+ 1im*w*taue)
    E6=c2*Y2(ke0,q0)./((1 .+ 1im*w*taue).*(1 .+ 1im*w*tau0))
    E7=c3*Y2(ki0,q0)./((1 .+ 1im*w*taue).*(1 .+ 1im*w*tau0))
    E8=c2*Y3(q0,k0,ke0)./(1 .+ 1im*w*taue)
    E9=c4*Y3(q0,k0,ki0)./(1 .+ 1im*w*taue)
    E=E1.-E2.-E3.-E4.-E5.-E6.-E7.-E8.-E9
    dvdv1=D1 .-D2 .- D3 .-E
    
    ##################################################
    # finally the upcrossing theory
    ##################################################
    theta0=Vth.-V0
    
    # oscillatory factors
    K1=0.5*dvdv1./dvdv0
    K2=-0.5*vv1./vv0
    K3=sqrt(pi/(2*dvdv0))*dVdt1
    K4=sqrt(pi/(2*dvdv0))*theta0*(vdv1/vv0)
    K5=(theta0^2/(2vv0))*(2*V1/theta0)
    K6=(theta0^2/(2vv0))*(vv1/vv0)
    ru1=ru0.*(K1.+K2.+K3.+K4.+K5.+K6)
    
    return V0,vv0,dvdv0,V1,dVdt1,vhe1,vhi1,vv1,vdv1,dvdv1,ru1

end

function DendrOscAsym(ae0::Float64,ai0::Float64,fkHz::Vector{Float64})
    
    w=2*pi*fkHz # for the high-freq case only
    nw=length(w)
    
    V0,vhe0,vhi0,vv0,dvdv0,ru0=DendrSteadyTheory(ae0,ai0)
    HL,He0,Hi0=1/tauL,ae0,ai0
    H0=HL+He0+Hi0
    V0=(HL*EL+He0*Ee+Hi0*Ei)/H0
    tau0=1/H0
    lam0=lamL*sqrt(tau0/tauL)
    Ee0=(Ee-V0)
    Ei0=(Ei-V0)
    
    #########################################################
    # Low-frequency limit
    #########################################################
    
    # V1
    V1L=tau0*Ee0*ones(nw)
    
    # vhe and vhi
    
    Qvhe1L=1-tau0*ae0*(1 + taue/(2*(taue+tau0)))
    vhe1L=tau0*(Ee0/4taue)*(lame/lam0)*sqrt(taue/(taue+tau0))*Qvhe1L*ones(nw)
    
    Qvhi1L=(Ee0 + Ei0*taui/(2*(taui+tau0)))
    vhi1L=-tau0*(1/4taui)*(lami/lam0)*tau0*ai0*sqrt(taui/(taui+tau0))*Qvhi1L*ones(nw)
    
    # vv
    ce=sqrt(taue/(taue+tau0))
    Qvve1L=(1-2*ae0*tau0)*(1-ce)-(ae0*tau0/2)*(1-ce*ce^2)
    vve1L=tau0*(Ee0^2/4)*(lame/lam0)*Qvve1L*ones(nw)
    ci=sqrt(taui/(taui+tau0))
    Qvvi1L=2*Ee0*Ei0*(1-ci)+(Ei0^2/2)*(1-ci*ci^2)
    vvi1L=-tau0*(tau0*ai0/4)*(lami/lam0)*Qvvi1L*ones(nw)
    vv1L=vve1L+vvi1L
        
    # dvdv
    Qdvdve1L=1-tau0*ae0*(2 + taue/(2*(taue+tau0)))
    dvdve1L=tau0*(Ee0^2/(4*taue^2))*(lame/lam0)*sqrt(taue/(taue+tau0))*Qdvdve1L*ones(nw)
    Qdvdvi1L=(2*Ee0*Ei0 + Ei0^2*taui/(2*(taui+tau0)))
    dvdvi1L=-tau0*(1/(4*taui^2))*(lami/lam0)*tau0*ai0*sqrt(taui/(taui+tau0))*Qdvdvi1L*ones(nw)
    dvdv1L=dvdve1L + dvdvi1L
    
    # now the upcrossing rate asymptotics
    theta0=Vth.-V0
    Qru1=theta0*V1L/vv0
    Qru2=0.5*(vv1L/vv0)*(theta0^2/vv0 - 1.0)
    Qru3=0.5*dvdv1L/dvdv0
    ru1L=ru0*(Qru1.+Qru2.+Qru3)
    
    #########################################################
    # High-frequency asymptotics
    #########################################################
    
    # variance of dvdt
    dvdv0e=(Ee0^2/(4taue^2))*(tau0*ae0)*(lame/lam0)*sqrt(taue/(tau0+taue))
    dvdv0i=(Ei0^2/(4taui^2))*(tau0*ai0)*(lami/lam0)*sqrt(taui/(tau0+taui))
    dvdv0=dvdv0e+dvdv0i
    
    # V1 and dV1dt asymptotics
    V1H=tau0*Ee0./(-w.^2*taue*tau0)
    dVdt1H=1im*w.*V1H
    
    # hev1 asymptotics
    vhe1H=tau0*(Ee0/(4taue))*(lame/lam0).*(1.0./((1im*w*taue/2).*sqrt.(1im*w*tau0)))
    
    # hiv1 asymptotics
    vhi1H=tau0*(Ei0/(4taui))*(ai0*tau0)*(lami/lam0).*(1.0./((1im*w*taue).*(1im*w*tau0)))*sqrt(taui/(taui+tau0))
    
    # vv1 asymptotics
    vv1He=-tau0*(Ee0^2/2)*(ae0*tau0)*(lame/lam0)*(1.0./((1im*w*tau0).*(1im*w*taue)))*(1-sqrt(taue/(tau0+taue)))
    vv1Hi=-tau0*(Ei0^2/2)*(ai0*tau0)*(lami/lam0)*(1.0./((1im*w*tau0).*(1im*w*taue)))*(1-sqrt(taui/(tau0+taui)))
    vv1H=vv1He.+vv1Hi
    
    # vdv1 asymptotics
    vdv1H=vv1H.*(1im*w/2)
    
    # dvdv1 asymptotics 
    dvdv1H=tau0*(Ee0^2/(2*taue^2))*(lame/lam0)./sqrt.(2*1im.*w*tau0)
    
    # finally the upcrossing rate asymptotics
    ru1H=ru0*0.5*dvdv1H/dvdv0
    
    return vhe1L,vhi1L,vv1L,dvdv1L,ru1L,dVdt1H,vhe1H,vhi1H,vdv1H,dvdv1H,ru1H
    
end

################################################################
################################################################
# Gaussian simulation functions
################################################################
################################################################

# Single run for the point-neuron in the Gaussian approximation
function PointGaussApproxSim(t::Vector{Float64},aet::Vector{Float64},ait::Vector{Float64},statestart)
    
    # get various parameters
    (Hestart,Histart,bUstart,bHestart,bHistart,ustart)=statestart
    dt=t[2]-t[1]
    nt=length(t); 
    ntp=nt+1; # extend to get starting point for next run
    
    # initialise
    He=zeros(ntp);  He[1]=Hestart
    Hi=zeros(ntp);  Hi[1]=Histart
    bU=zeros(ntp);    bU[1]=bUstart
    bHe=zeros(ntp);   bHe[1]=bHestart
    bHi=zeros(ntp);   bHi[1]=bHistart
    u=zeros(ntp);   u[1]=ustart
    ru=zeros(nt)
    
    rV=zeros(nt) # with reset  
    rU=zeros(nt) # upcrossing but has cross terms like v*he
    ru=zeros(nt) # upcrossing in guassian approximation
    
    # precalculate noise
    flucte=randn(nt).*sqrt.(aet.*kappae)/sqrt(dt)
    flucti=randn(nt).*sqrt.(ait.*kappai)/sqrt(dt)
    
    for k=1:nt
        
        ##########################################
        # Deterministic part
        ##########################################
        bHe[k+1]=bHe[k] .+ (dt/taue)*(aet[k] .- bHe[k])
        bHi[k+1]=bHi[k] .+ (dt/taui)*(ait[k] .- bHi[k])
        bHeU=bHe[k].*(Ee .-bU[k])
        bHiU=bHi[k].*(Ei .-bU[k])
        bHLU=(EL .-bU[k])/tauL
        bU[k+1]=bU[k] .+ dt*( bHLU .+ bHeU .+ bHiU)
        
        ##########################################
        # synaptic dynamics are used by all cases
        ##########################################
        He[k+1]=He[k]+(dt/taue)*(flucte[k] .+ aet[k] .-He[k])
        Hi[k+1]=Hi[k]+(dt/taui)*(flucti[k] .+ ait[k] .-Hi[k])
        
        ##########################################
        # upcrossing in gaussian approx
        ##########################################
        uHeHi=-u[k].*(bHe[k] .+ bHi[k] .+ 1/tauL)
        EeHeEiHi=(Ee .- bU[k]).*(He[k] .-bHe[k]) .+(Ei .- bU[k]).*(Hi[k] .-bHi[k])
        u[k+1]=u[k] .+ dt*(uHeHi .+ EeHeEiHi)
        if (((u[k]+bU[k])<Vth))&((u[k+1]+bU[k+1])>Vth)
            ru[k]=1.0/dt
        end
        
    end
    
    # tidy up arrays 
    statestart=(He[end],Hi[end],bU[end],bHe[end],bHi[end],u[end])
    
    return bU[1:end-1],u[1:end-1],He[1:end-1],Hi[1:end-1],ru,statestart

end

# Runs multiple simulations and returns various averages
function PointSweeps(t::Vector{Float64},aet::Vector{Float64},ait::Vector{Float64},ns::Int64)
    
    dt=t[2]-t[1]
    nt=length(t)

    # Thermalise for 1sec
    statestart=(0,0,EL,0,0,0)
    ~,~,~,~,~,statestart=PointGaussApproxSim(t,aet,ait,statestart)
    
    #################################################
    # Main run
    #################################################

    rus=zeros(nt,ns) 
    Hes=zeros(nt,ns)
    His=zeros(nt,ns)
    us=zeros(nt,ns)
    dus=zeros(nt,ns)

    for s=1:ns
        ~,us[:,s],Hes[:,s],His[:,s],rus[:,s],statestart=PointGaussApproxSim(t,aet,ait,statestart);
        dus[:,s]=diff([us[:,s] ; statestart[6]])/dt
    end
    
    mru=mean(rus,dims=2)
    hes,his=Hes.-aet,His.-ait
    muu=mean(us.*us,dims=2)
    muhe=mean(us.*hes,dims=2)
    muhi=mean(us.*his,dims=2)
    mdudu=mean(dus.*dus,dims=2)

    return mru,muu,muhe,muhi,mdudu
    
end

function DendrGaussApproxSim(x::Vector{Float64},t::Vector{Float64},aet::Vector{Float64},ait::Vector{Float64},statestart)
    
    dt=t[2]-t[1]; nt=length(t); ntp=nt+1
    dx=x[2]-x[1]; nx=length(x)

    (Hestart,Histart,bUstart,bHestart,bHistart,ustart)=statestart
    
    # Initialise
    u=zeros(ntp,nx);   u[1,:]=ustart
    He=zeros(ntp,nx); He[1,:]=Hestart
    Hi=zeros(ntp,nx); Hi[1,:]=Histart
    bU=zeros(ntp);    bU[1]=bUstart
    bHe=zeros(ntp);   bHe[1]=bHestart
    bHi=zeros(ntp);   bHi[1]=bHistart

    ru=zeros(nt)
    
    # useful constants
    HL=1/tauL
    D=lamL^2/tauL
    
    flucte=sqrt.(aet.*lame)/sqrt(dx*dt)
    flucti=sqrt.(ait.*lami)/sqrt(dx*dt)
    
    for k=1:nt
        
        ##########################################
        # Deterministic part
        ##########################################
        bHe[k+1]=bHe[k] .+ (dt/taue)*(aet[k] .- bHe[k])
        bHi[k+1]=bHi[k] .+ (dt/taui)*(ait[k] .- bHi[k])
        bHeU=bHe[k].*(Ee .-bU[k])
        bHiU=bHi[k].*(Ei .-bU[k])
        bHLU=(EL .-bU[k])/tauL
        bU[k+1]=bU[k] .+ dt*( bHLU .+ bHeU .+ bHiU)
        
        ##########################################
        # Gaussian synaptic input 
        ##########################################
        He[k+1,:]=He[k,:] .+ (dt/taue)*(aet[k] .- He[k,:] .+ flucte[k]*randn(nx))
        Hi[k+1,:]=Hi[k,:] .+ (dt/taui)*(ait[k] .- Hi[k,:] .+ flucti[k]*randn(nx))
       
        ##########################################
        # Upcrossing in gaussian approx
        ##########################################
        uHeHi=-u[k,:].*(bHe[k] .+ bHi[k] .+ 1/tauL)
        EeHeEiHi=(Ee .- bU[k]).*(He[k,:] .-bHe[k]) .+(Ei .- bU[k]).*(Hi[k,:] .-bHi[k])
        Dd2udx2=D*diff(diff([u[k,end];u[k,:];u[k,1]]))/dx^2 # wrap around
        u[k+1,:]=u[k,:] .+ dt*(uHeHi .+ EeHeEiHi .+ Dd2udx2)

        spikes=length(findall(((u[k,:] .+bU[k]).<Vth).&((u[k+1,:] .+bU[k+1]).>=Vth)))
        ru[k]=spikes/(dt*nx)

    end
    
    # Tidy up arrays and return
    statestart=(He[end,:],Hi[end,:],bU[end],bHe[end],bHi[end],u[end,:])
    
    return bU[1:end-1],u[1:end-1,:],He[1:end-1,:],Hi[1:end-1,:],ru,statestart

end


# Runs multiple simulations and returns the average
function DendrSweeps(x::Vector{Float64},t::Vector{Float64},aet::Vector{Float64},ait::Vector{Float64},ns::Int64)

    dt=t[2]-t[1]; nt=length(t); ntp=nt+1
    dx=x[2]-x[1]; nx=length(x)
    sh=Int(round(nx/2,digits=0))

    # Thermalise for 1sec
    statestart=(zeros(nx),zeros(nx),EL,0,0,zeros(nx))
    ~,~,~,~,~,statestart=DendrGaussApproxSim(x,t,aet,ait,statestart)
    
    #################################################
    # Main run
    #################################################
    
    ru=zeros(nt,ns) 
    uus=zeros(nt,ns)
    hehes=zeros(nt,ns)
    hihis=zeros(nt,ns)
    uhes=zeros(nt,ns)
    uhis=zeros(nt,ns)
    dudus=zeros(nt,ns)
    
    # The full run
    for s=1:ns
        bU,u,He,Hi,ru[:,s],statestart=DendrGaussApproxSim(x,t,aet,ait,statestart)
        uhes[:,s]=mean(u.*(He.- aet),dims=2)
        uhis[:,s]=mean(u.*(Hi.- ait),dims=2)
        uus[:,s]=mean(u.^2,dims=2)
        dudus[:,s]=mean( (diff([u ; statestart[6]'],dims=1)/dt).^2,dims=2)
    end
    
    mru=mean(ru,dims=2)
    muu=mean(uus,dims=2)
    muhe=mean(uhes,dims=2)
    muhi=mean(uhis,dims=2)
    mdudu=mean(dudus,dims=2)

    return mru,muu,muhe,muhi,mdudu
    
end

