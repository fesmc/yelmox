!................................................................................!
!                                                                                ! 
! Adapting the box model of Zickfeld/Rahsmtorf clim. dyn. 2004 to fortran 90     !
!                                                                                !
!              Jorge Alvarez-Solas | UCM 2019                                    ! 
!                                                                                !
! And "liberating" it from calibration to PD as in the 2004 paper                !                                                                           ! 
!................................................................................!

!............................................................................
!
! script boxmod shows time-dependent solution of the 
! Zickfeld/Rahmstorf box model. 
! The model equations are (see Ocean Dyn. paper):
! (only valid for positive m!)
!
!          dT1/dt = lam1*(Tr1-T1) + m/V1*(T4-T1)
!          dT2/dt = lam2*(Tr2-T2) + m/V2*(T3-T2)
!          dT3/dt = lam3*(Tr3-T3) + m/V3*(T1-T3)
!          dT4/dt =                 m/V4*(T2-T4)
!
!          dS1/dt =  S0*f1/V1       + m/V1*(S4-S1)
!          dS2/dt = -S0*f2/V2       + m/V2*(S3-S2)
!          dS3/dt =  S0*(f2-f1)/V3  + m/V3*(S1-S3)
!          dS4/dt =                   m/V4*(S2-S4)
!
!          m = k * (beta*(S2-S1) - alpha*(T2-T1))  ;  m = flow rate
!
! Box model parameters:
!
!         f1,f2 = freshwater fluxes 
!         k = empirical flow constant
!         alpha, beta = expansion coefficients
!         Tri = restoring temperatures
!         lami = thermal relaxation constants
!         Vi = volume ratio of box i
!.............................................................................


program zira04

    implicit none

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the working precision (wp) of the library (sp,dp)
    integer,  parameter :: wp = sp


    real(wp) :: time_init, time_end, dt, time
    integer  :: n, ntot, ncid

    !           physical constants
    !---------------------------------------------------
    real(wp) :: alpha                    ! 1/K
    real(wp) :: beta                     ! 1/psu
    real(wp) :: csp                      ! specific heat capacity, J/K/Kg
    real(wp) :: ro                       ! density, Kg/m/m/m
    real(wp) :: secy                     ! seconds per year
    real(wp) :: PI=3.141592

    !           model parameters
    !--------------------------------------------------
    real(wp) :: gamma             ! thermal coupling constant, W/m^2/K
    real(wp) :: S0                   ! reference salinity, psu
    real(wp) :: F2sv                ! freshwater transport, Sv
    
    real(wp) :: boxvol              ! reference box volume, m^3
    real(wp) :: V1                   ! volume ratio of box i w.r.t boxvol
    real(wp) :: V2 
    real(wp) :: V3 
    real(wp) :: V4 
    
    real(wp) :: A1               ! box surface area, m^2
    real(wp) :: A2 
    real(wp) :: A3 
    
    real(wp) :: Dx1  !V1*boxvol/A1   ! depth of box i, m
    real(wp) :: Dx2  !V2*boxvol/A2 
    real(wp) :: Dx3 
    
    real(wp) :: lam1   ! thermal relaxation constants 
    real(wp) :: lam2 
    real(wp) :: lam3 
    
    real(wp) :: tau1                   ! restoring times, yr
    real(wp) :: tau2 
    real(wp) :: tau3 

    real(wp) :: tauref1                   ! jalv: restoring times for salinities, yr
    real(wp) :: tauref2
    real(wp) :: tauref3
    real(wp) :: tauref4


    ! Tri and k are derived from a fit to the CLIMBER hysteresis
    real(wp) :: Tr1              ! degC
    real(wp) :: Tr2              ! degC
    real(wp) :: Tr3              ! degC
    real(wp) :: k                ! hydraulic constant, 1/yr

    ! Srefs | jalv ad hoc values to ensure a kind of salinity conservation at the equilibrium and not stupid salinity values
    real(wp) :: Sref1              ! psu
    real(wp) :: Sref2              ! psu
    real(wp) :: Sref3              ! psu
    real(wp) :: Sref4              ! psu


    ! conversion factor of m and F from 1/yr  into Sv    
    real(wp) :: fact 
    
    !                climate parameters
    !---------------------------------------------------
    ! hydrological sensitivities derived from CLIMBER [Sv/K]
    ! hysens = DF_ideltaDTreg_i according to Eq. 26-27
    real(wp) :: hysens1                ! DF1/DT3
    real(wp) :: hysens2                ! DF2/DT2
    real(wp) :: hysens3                ! DF3/DT2 (F3 = meltwater)     
    real(wp) :: hysens4         

    ! regional "patterns" of temperature increase derived from CLIMBER
    ! pati = DTreg_i/DTglobal according to Eq. 25
    real(wp) :: pat2     ! North Atlantic
    real(wp) :: pat1     ! South Atlantic
    real(wp) :: pat3     ! Tropical Atlantic
    real(wp) :: patSH    ! Southern Hemisphere
    real(wp) :: patNH    ! Northern Hemisphere

    !                prescribed climatic changes
    !---------------------------------------------------
    real(wp) :: wyrs               ! years of sustained warming
    real(wp) :: totw               ! total amount of warming
    real(wp) :: trend              ! degrees per century 
        
    ! -----------------------EQUILIBRIUM CONDITIONS ---------------------        
    ! assign the equilibrium flow strength:
    ! msv = Equilibrium overturning in Sv   
    real(wp) ::    msve 
    real(wp) ::    me
    
    ! compute the corresponding steady-state freshwater flux f1:
    real(wp) ::  deltaT
    real(wp) ::  f1
    real(wp) ::  F1sv               ! In Sv            
    real(wp) ::  f2   

    ! compute the other equilibrium values
    real(wp) ::  T1e
    real(wp) ::  T2e
    real(wp) ::  T3e
    real(wp) ::  T4e
    real(wp) ::  S1e
    real(wp) ::  S2e
    real(wp) ::  S3e
    real(wp) ::  S4e
    
    real(wp) ::  DS
    real(wp) ::  DSN
    real(wp) ::  DSS
    
    ! -------------------------- INITIAL CONDITIONS ------------------------
    !
    ! now add a perturbation deltam to the equilibrium flow
    ! to initialize the model away from equilibrium
    !
    real(wp) :: deltam
    real(wp) :: mi
    real(wp) :: deltaTi
    real(wp) :: fi

    ! variables entering the differential equations..........................
    real(wp) :: h1
    real(wp) :: h2
    real(wp) :: h3
    real(wp) :: h4

    real(wp) :: FWF1         ! Additional freshwater fluxes (e.g. from ice melting) 
    real(wp) :: FWF2         ! in Sv
    real(wp) :: FWF3         !
    real(wp) :: FWF4         !
    
    real(wp) :: sf1         ! Salt fluxes from additional freshwater fluxes (e.g. from ice melting) 
    real(wp) :: sf2         ! in salt flux units
    real(wp) :: sf3         !
    real(wp) :: sf4 

    real(wp) ::  W1
    real(wp) ::  W2
    real(wp) ::  W3
    real(wp) ::  Wp1
    real(wp) ::  Wp2
    real(wp) ::  Wp3
    
    ! prognsotic variables ..................................................
    real(wp) ::  T1
    real(wp) ::  T2
    real(wp) ::  T3
    real(wp) ::  T4
    real(wp) ::  S1
    real(wp) ::  S2
    real(wp) ::  S3
    real(wp) ::  S4
    real(wp) ::  msv             ! In Sv       Msv = m*fact
    real(wp) ::  m               ! In 1/yr 

    !________________________________________________________________________
    !                                                                        !
    !                         END of DECLARING                               !
    !________________________________________________________________________!

    !........................................................................!

    !________________________________________________________________________
    !                                                                        !
    !                     INIT of assigning values                           !
    !________________________________________________________________________!
 
    !           physical constants
    !---------------------------------------------------
    alpha= 0.00017
    beta = 0.0008
    csp = 4000
    ro = 1025
    secy = 86400*365

    !           model parameters
    !--------------------------------------------------   
    gamma = 23.1             
    S0 = 35                  
    F2sv = 0.065               
    
    boxvol=1e17              

    V1=1.1                  
    V2=0.4
    V3=0.68
    V4=0.04*1.36

    A1= 4.6e13              
    A2= 1e13
    A3= 6.8e13

    Dx1=3000 
    Dx2=3000  
    Dx3=V3*boxvol/A3

    lam1= gamma*secy/csp/ro/Dx1   
    lam2= gamma*secy/csp/ro/Dx2
    lam3= gamma*secy/csp/ro/Dx3

    ! jalv . test for increasing the time reference for the southern ocean to change temperatures 
    lam1 = lam1 * 0.01 
    !lam1 = lam1 * 10.0 
    !lam2 = lam2 * 10.0 
    !lam3 = lam3 * 10.0

    tau1=1/lam1                   
    tau2=1/lam2
    tau3=1/lam3

    ! Tri and k are derived from a fit to the CLIMBER hysteresis
    Tr1 = 6.6             
    Tr2 = 2.7             
    Tr3 = 11.7           
    ! jalv: testing a change in Trefs
    Tr1 = 10.0
    Tr2 = 5.0
    Tr3 = 20.0

    k = 25.4           

    ! Srefs | jalv
    Sref1 = 35.0
    Sref2 = 35.0
    Sref3 = 35.0
    Sref4 = 35.0

    Sref1 = 33.0
    Sref2 = 32.0
    Sref3 = 37.0
    Sref4 = 35.0

    ! taurefs for salinity | jalv
    tauref1 = 2000.0         ! years
    tauref2 = 1500.0
    tauref3 = 500.0
    tauref4 = 500.0


    ! conversion factor of m and F from 1/yr  into Sv
    fact=boxvol/1e6/secy
        
    !                climate parameters
    !---------------------------------------------------

    ! hydrological sensitivities derived from CLIMBER [Sv/K]
    ! hysens = DF_i/DTreg_i according to Eq. 26-27    
    hysens1 = -0.005               
    hysens2 = 0.013               
    hysens3 = 0.00                    
    hysens4 = 0.0
    
    ! regional "patterns" of temperature increase derived from CLIMBER
    ! pati = DTreg_i/DTglobal according to Eq. 25
    pat2 = 1.07     ! North Atlantic
    pat1 = 0.86     ! South Atlantic
    pat3 = 0.79     ! Tropical Atlantic
    patSH = 0.93    ! Southern Hemisphere
    patNH = 1.07    ! Northern Hemisphere

    !                prescribed climatic changes
    !---------------------------------------------------
    wyrs = 150               ! years of sustained warming
    totw = 4.5               ! total amount of warming
    trend = totw*100/wyrs    ! degrees per century    
        
    ! -----------------------EQUILIBRIUM CONDITIONS ---------------------        
    ! assign the equilibrium flow strength:
    ! Me = Equilibrium overturning in Sv      
    msve  = 22.6     ! in sverdrups
    me = Msve/fact    ! jalv: not in Sv 

    ! compute the corresponding steady-state freshwater flux f1:                     
    deltaT=(me*lam1*(lam3*(Tr3-Tr1)/V2+lam2*(Tr2-Tr1)/V3)+lam1*lam2*lam3*(Tr2-Tr1))/((me/V1+lam1)*(me/V2+lam2)*(me/V3+lam3)-me*me*me/V1/V2/V3)

    f1 = -me* (me + k*alpha*deltaT)/k/beta/S0 
    F1sv = f1*fact
    f2   = F2sv/fact

    ! compute the other equilibrium values           
    T1e = (me*deltaT/V1+lam1*Tr1)/lam1           
    T3e = (lam3*Tr3+me/V3*T1e)/(lam3+me/V3)  
    T2e = (lam2*Tr2+me/V2*T3e)/(lam2+me/V2)  
    T4e = T2e

    S1e = S0*f1/me 
    S2e = 0 
    S3e = S0*f2/me  
    S4e = 0 

    DS  = S2e-S1e 
    DSN = S3e-S2e
    DSS = S3e-S1e

    ! -------------------------- INITIAL CONDITIONS ------------------------
    !
    ! now add a perturbation deltam to the equilibrium flow
    ! to initialize the model away from equilibrium
    !
       
    deltam = 0          ! initial deviation of flow from equilibrium in Sv
    mi = me + deltam/fact 
    deltaTi =(mi*lam1*(lam3*(Tr3-Tr1)/V2+lam2*(Tr2-Tr1)/V3)+lam1*lam2*lam3*(Tr2-Tr1))/((mi/V1+lam1)*(mi/V2+lam2)*(mi/V3+lam3)-mi*mi*mi/V1/V2/V3) 

    fi = -mi* (mi + k*alpha*deltaTi)/k/beta/S0  
        
    ! Initialization of solution vector:
    m =0.0 
    T1=0.0
    T2=0.0 
    T3=0.0 
    T4=0.0 
    S1=0.0 
    S2=0.0 
    S3=0.0 
    S4=0.0 

    ! compute initial values of T,S (assuming initial state is a steady state)
    m  = mi 
    T1  = (mi*deltaTi/V1+lam1*Tr1)/lam1              ! Eq. 8
    T3  = (lam3*Tr3+mi/V3*T1 )/(lam3+mi/V3)    ! Eq. 9
    T2  = (lam2*Tr2+mi/V2*T3 )/(lam2+mi/V2)    ! Eq. 10
    T4  = T2  
        
    S1  = S0*fi/mi 
    S2  = 0 
    S3  = S0*f2/mi 
    S4  = 0 

    ! jalv: test with initial Si = S0:
    S1  = S0
    S2  = S0
    S3  = S0
    S4  = S0

   
    ! open output  writing files
    OPEN(11,FILE='/home/jalvarez/stommel/output/sal1.dat',STATUS='unknown')
    OPEN(12,FILE='/home/jalvarez/stommel/output/sal2.dat',STATUS='unknown')
    OPEN(13,FILE='/home/jalvarez/stommel/output/sal3.dat',STATUS='unknown')
    OPEN(14,FILE='/home/jalvarez/stommel/output/temp1.dat',STATUS='unknown')
    OPEN(15,FILE='/home/jalvarez/stommel/output/temp2.dat',STATUS='unknown')
    OPEN(16,FILE='/home/jalvarez/stommel/output/temp3.dat',STATUS='unknown')
    OPEN(17,FILE='/home/jalvarez/stommel/output/THC.dat',STATUS='unknown')
    OPEN(18,FILE='/home/jalvarez/stommel/output/FWF.dat',STATUS='unknown')
    OPEN(19,FILE='/home/jalvarez/stommel/output/FWFs.dat',STATUS='unknown')
    OPEN(20,FILE='/home/jalvarez/stommel/output/hysteresis.dat',STATUS='unknown')
    OPEN(21,FILE='/home/jalvarez/stommel/output/hysteresis_south.dat',STATUS='unknown')


    ! ========================================================
    ! User parameters 

    ! Define the time steps 
    time_init =   0000.0
    time_end  =   500000.0
    dt        =       1.0
    ! ========================================================

    ! Determine total timesteps
    ntot = ceiling((time_end-time_init)/dt)

    ! Advance timesteps
    do n = 1, ntot
        ! Get current time 
        time = time_init + n*dt
       
        ! ---------------------------------------------------------------------------!
        !               * * * Compute THC Response * * *
        ! ---------------------------------------------------------------------------!
    
        
        !               define forcing time series |  jalv: left behind from former matlab program
        ! ---------------------------------------------------------------------------!
        
        ! add linear warming trend to Tr2
        !W2= [0:nt]*pat2*trend*dt/100 
        !W2(wyrs+1:ns)=ones(1,ns-wyrs)*W2(wyrs+1)   ! stop warming after some time
        W2 = 0.0  
       
        ! add linear warming trend to Tr1
        !W1= [0:nt]*pat1*trend*dt/100 
        !W1(wyrs+1:ns)=ones(1,ns-wyrs)*W1(wyrs+1)   ! stop warming after some time
        W1 = 0.0
            
        ! add linear warming trend to Tr3
        !W3=[0:nt]*pat3*trend*dt/100 
        !W3(wyrs+1:ns)=ones(1,ns-wyrs)*W3(wyrs+1)   ! stop warming after some time
        W3 = 0.0
            
        ! add linear warming trend to T_SH
        !W1p= [0:nt]*patSH*trend*dt/100 
        !W1p(wyrs+1:ns)=ones(1,ns-wyrs)*W1p(wyrs+1) 
    
        ! add linear warming trend to T_NH
        !W2p= [0:nt]*patNH*trend*dt/100 
        !W2p(wyrs+1:ns)=ones(1,ns-wyrs)*W2p(wyrs+1) 

        !h1=-ones(1,ns)*f1*S0      ! salt flux
        !h2=-ones(1,ns)*f2*S0      ! salt flux
    
        h1=-f1*S0      ! salt flux
        h2=-f2*S0      ! salt flux
    
        !add freshwater perturbation to f1
        !h1=h1(1:ns)-hysens1*W1p/fact*S0             

        !add freshwater perturbation to f2
        !h2=h2(1:ns)-hysens2*W2p/fact*S0             
        
        ! add meltwater perturbation (f3)
        !h3=-hysens3*W2p/fact*S0 
        h3 = 0.0
        
        ! "Latif" effect
        !h4=-hysens4*W3/fact*S0 
        h4 = 0.0

        ! Additional FWFs (in Sv) and then converted in salt fluxes:
        
        ! a) constant values
        FWF1 = 0.0
        FWF2 = 0.0
        FWF3 = 0.0
        FWF4 = 0.0

        !FWF1 = -0.12
        !FWF2 = 0.0
        !FWF3 = 0.12
        !FWF4 = 0.0

        !h1 = 0.0
        !h2 = 0.0 
 
        !b) evolving values
        h1 = 0.0
        h2 = 0.0
        !FWF1=0.3*sin(2*PI*time/3000) ! FWF SOUTH centered in 0 ==> the minimum amplitude for stopping THC seems to be in PD (in glacial) ?? (??)     Sv
        FWF2=0.25*sin(2*PI*time/3000) ! FWF NORTH centered in 0 ==> the minimum amplitude for stopping THC seems to be in PD (in glacial) 0.2 (0.015) Sv
        !FWF3 = -FWF2
        !   c) Hysteresis experiment :
        ! North
        if (.false.) then
            FWF2=0.
            h1=0.0
            h2=0.0
            if (time.ge.100000) then
                FWF2=(1./200000.)*(time-100000.)
            endif
            if (time.ge.200000) then
                FWF2=0.5-(1./200000.)*(time-200000.)
            endif
            if (time.ge.400000) then
                FWF2=-0.5+(1./200000.)*(time-400000.)
            endif
            if (time.ge.500000.) then
                FWF2=0.
            endif
            !FWF3= -FWF
        endif

        ! South
        if (.false.) then
            FWF1=0.
            h1=0.0
            h2=0.0
            if (time.ge.100000) then
                FWF1=(-1./100000.)*(time-100000.)
            endif
            if (time.ge.200000) then
                FWF1=-1.0+(1./100000.)*(time-200000.)
            endif
            if (time.ge.400000) then
                FWF1=1.0-(1./100000.)*(time-400000.)
            endif
            if (time.ge.500000.) then
                FWF1=0.
            endif
        endif


        sf1 = -FWF1*S0/fact
        sf2 = -FWF2*S0/fact
        sf3 = -FWF3*S0/fact
        sf4 = -FWF4*S0/fact
       
        ! jalv mimicking glacial a lo cutre (without calibrating the equilibrium quantities for the given change in the Trefs first)
        ! jalv: does not seem to work, if I cool the trefs no shut down of the AMOC is simulated even for 1 Sv
        Tr3 = Tr3 - 0.0
        Tr2 = Tr2 - 0.0
        Tr1 = Tr1 - 0.0

        ! jalv: trial to make the hysteresis loop to FWF2 make sense:
        ! by default, as soon as M=0 (if M limited to postive) the AMOC does noet recover 
        ! and if M is allowed to become negative, the model crashes (solution diverges)
        !m = abs(m)

        ! jalv: trial for having Srefs that account for the fact that the FWF have changed the total salinity of the ocean
        !       in this case the Srefs wuld be the mean salinity of all the boxes in order to adjust to the fact that has beeen a net FWF input (even if negative)
        !       for simplicity I assume (as I do by default that all teh Srefs are the same) (if inititally the Srefs are not the same. i.e more saline in some places ...
        !        ... then this averaging has to be rethough to account for a potential spatial variability of the Srefs
        !sref1 = (S1 + S2 + S3 + S4) / 4.0
        !sref2 = sref1
        !sref3 = sref1
        !sref4 = sref1


        !jalv: be careful, by default it does not work, cause if all the Srefs are the mean, if one salinity is reduced through FWF input , then the mean is lowered and all...
        !      ... the salinities want to go to that lower value, eventaully driving all to 0.
        !      So, nbew trial defining the Sref of one box as the mean of the others  
        !sref1 = (     S2 + S3 + S4) / 6.0 + S0/2.0
        !sref2 = (S1 +    + S3 + S4) / 6.0 + S0/2.0
        !sref3 = (S1 + S2 +    + S4) / 6.0 + S0/2.0
        !sref4 = (S1 + S2 + S3     ) / 6.0 + S0/2.0
        ! --------------------------compute m(t) ---------------------------------------!
        
        ! solve DEs: jalv: the following is the way it was in the former model
        ! jalv: note that h3 and h4 do not appear in the original equations
        ! so I assume h3 and h4 appear because prepared to be used with the perturbation experiments
        ! if h3 and h4 = 0, the original equations are recovered    
        ! the same applies for the temperature perturbations W1, W2 and W3
        if (.false.) then
            T1 = T1 + (lam1*(Tr1+W1-T1) + m/V1*(T4-T1)) * dt
            T2 = T2 + (lam2*(Tr2+W2-T2) + m/V2*(T3-T2)) * dt 
            T3 = T3 + (lam3*(Tr3+W3-T3) + m/V3*(T1-T3)) * dt 
            T4 = T4 + (                   m/V4*(T2-T4)) * dt 
            
            S1 = S1 + (   -h1/V1           + m/V1*(S4-S1)) * dt 
            S2 = S2 + (   (h2+h3)/V2       + m/V2*(S3-S2)) * dt 
            S3 = S3 + (   (h1-h2-h4)/V3    + m/V3*(S1-S3)) * dt 
            S4 = S4 + (                      m/V4*(S2-S4)) * dt 
        end if

        ! solve DEs: jalv adding additional FWFs and restoring salinities
        if (.true.) then
            T1 = T1 + (lam1*(Tr1+W1-T1) + m/V1*(T4-T1)) * dt    ! jalv : ojo con el 0.01 !!!
            T2 = T2 + (lam2*(Tr2+W2-T2) + m/V2*(T3-T2)) * dt
            T3 = T3 + (lam3*(Tr3+W3-T3) + m/V3*(T1-T3)) * dt
            T4 = T4 + (                   m/V4*(T2-T4)) * dt
    
            S1 = S1 + ((sref1 - S1)/tauref1    +  sf1/V1   -h1/V1            + m/V1*(S4-S1)) * dt
            S2 = S2 + ((sref2 - S2)/tauref2    +  sf2/V2   +(h2+h3)/V2       + m/V2*(S3-S2)) * dt
            S3 = S3 + ((sref3 - S3)/tauref3    +  sf3/V3   +(h1-h2-h4)/V3    + m/V3*(S1-S3)) * dt
            S4 = S4 + ((sref4 - S4)/tauref4    +  sf4/V4                     + m/V4*(S2-S4)) * dt
        end if

        !T1 = max(T1, 0.0)
        !T2 = max(T2, 0.0)
        !T3 = max(T3, 0.0)
        !T4 = max(T4, 0.0)

        S1 = max(S1, 0.0)
        S2 = max(S2, 0.0)
        S3 = max(S3, 0.0)
        S4 = max(S4, 0.0)


        m = k * (beta*(S2-S1) - alpha*(T2-T1))   
        ! jalv testing the AMOC as the tropics and north density gradient (seems to help for AMOC reactivity when forcing cyclically the south, but not for Antarctic temperatures)
        !m = k * (beta*(S3-S1) - alpha*(T3-T1))
 
        m = max(m, 0.0)
 
        msv = m * fact

        !msv = max(msv,-1.5)

        !m = msv/fact

        !write(*,*) msv
        !write(*,*) FWF1

        !write(*,*) T1 
        !write(*,*) FWF2, msv
        !write(*,*) FWF1, msv

        WRITE(11,*)S1
        WRITE(12,*)S2
        WRITE(13,*)S3
        WRITE(14,*)T1
        WRITE(15,*)T2
        WRITE(16,*)T3
        WRITE(17,*)msv
        WRITE(18,*)FWF2
        WRITE(19,*)FWF1
        WRITE(20,*)FWF2,msv
        WRITE(21,*)FWF1,msv
         
    end do


   CLOSE(11)
   CLOSE(12)
   CLOSE(13)
   CLOSE(14)
   CLOSE(15)
   CLOSE(16)
   CLOSE(17)
   CLOSE(18)
   CLOSE(19)
   CLOSE(20)
   CLOSE(21)


end program


