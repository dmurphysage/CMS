      program CMS_Haz45

c     This program will compute the conditional mean spectra for conditioning periods
c     from PSHA runs. This version is compatible with the output files from Haz45.2
c     which only outputs the mean hazard code combined over all attenuation models. 

c     Version 45.2, last modified 04/2017

      implicit none
      include 'cms.h'
      
      real*8 haz(MAX_PROB, MAX_ATTENTYPE, MAX_ATTEN, MAX_INTEN)
      real*8 hazmean1(MAX_ATTENTYPE,MAX_INTEN)      

      integer iPer, nAttenType, HWFlag, vs30_class, forearc, nScenario, 
     1        intflag(4,MAX_PROB), iScen, iInten, nProb, jProb, kType,
     2        nGM_model(MAX_PROB,MAX_ATTENTYPE), jType, iAtten, iPer2, 
     3        attenType(MAX_FLT), jCalc(MAX_PROB,4,MAX_ATTEN), nINten1, 
     4        nInten(MAX_PROB)      
      real version, varAdd(MAX_PROB,4,MAX_ATTEN), period(MAX_PROB), ratio,
     1     lgInten, sigmaY, sum, sum1, sum2, wt_deagg(MAX_ATTENTYPE,MAX_ATTEN),
     2     hazLevel, haz_GMPE(MAX_INTEN), haz10, period1(4,MAX_PROB), tau, phi,
     3     UHS1(MAX_PROB), epsilonstar, CMS(MAX_PROB), rRup, rJB, rSeismo, 
     4     rHypo, Rx, Ry0, mag, ftype, vs30, hypoDep, AR, dip, Z1, Z15, Z25, 
     5     zTOR, theta_site, RupWidth, testInten(MAX_PROB, MAX_INTEN),
     6     gmScale(MAX_PROB,4,MAX_ATTEN), gm_wt(MAX_PROB,4,MAX_ATTEN),
     7     cfcoefrrup(MAX_Atten,11), cfcoefrjb(MAX_Atten,11), med(MAX_PROB), 
     8     sig(MAX_PROB)      
      character*80 filein, file1, dummy, attenName(4,MAX_ATTEN)
                  
      write (*,*) '******************************'
      write (*,*) '*      CMS Code for GMC      *'
      write (*,*) '*  compatible with Haz45.2   *'
      write (*,*) '*   Tagged April 10, 2017    *'
      write (*,*) '******************************'

      write (*,*) 'Enter the input filename.'
      
      read (*,'(a80)') filein
      open (31,file=filein,status='old',ERR=2000)

      read (31,*,ERR=2001) Hazlevel
      read (31,*,ERR=2003) jProb

c     Read the hazard run file to get the logic tree weights for the GMPEs 
c     and the testInten values
      call RdInput (nInten, testInten, nGM_model, nattentype,attenType,
     1               nProb, GM_wt, period, jCalc, gmscale, VarAdd, version)
       
c     Read the out6 file
      write (*,'( 2x,''reading logic tree file out6'')')
      nInten1 = nInten(1)
      call read_logichaz_out6 ( haz, nProb, nAttenType, nGM_model, nInten1)
      write (*,'( 2x,''out of logichaz_out6'')')
      
c     Read the out3 file and compute the UHS at the desired haz level
      call RdOut3 ( UHS1, HazLevel ) 
      write (*,'( 2x,''out of out3'')')
      
c     Open the output file
      read (31,'( a80)') file1
      open (50,file=file1,status='new')
      write (50,*) ' *** Output file from program CMS  *** '
      write (50,*) '        *** Version 45.2 ***           '    
      write (50,*) 
      write (50,'(a17,2x,a80)') ' Input filename: ', filein 
      write (50,*) 
      write (50,'( ''    period    AttenType  iAtten  HazRatio  Wt_LT     Wt_D'')')
       
c      Compute the mean hazard for each AttenType
       do jType=1,nAttenType
         do iInten=1,nInten(jProb)
           sum = 0.
           do iAtten=1,nGM_Model(jProb,jType)
             sum = sum + haz(jProb,jType,iAtten,iInten)*gm_wt(jProb,jType,iAtten)
           enddo
           hazMean1(jType,iInten) = sum
         enddo
       enddo
       
c      Compute the hazard for each GMPE, adding the mean hazard from other AttenTypes
       do jType=1,nAttenType
        do iAtten=1,nGM_Model(jProb,jType)
         do iInten=1,nInten(jProb)
           haz_GMPE(iInten) =  haz(jProb,jType,iAtten,iInten)
c          Now add the mean from the other Atten Types
           do kType=1,nAttenType
            if ( kType .ne. jType) then
             haz_GMPE(iInten) = haz_GMPE(iInten) + hazMean1(kType,iInten)
            endif
           enddo
         enddo
 
c        Interpolate to find the hazard at the UHS level
         iPer = jProb
         do iInten=2,nInten(iPer)
          if ( UHS1(iPer) .ge. testInten(iPer,iInten-1) .and. 
     1       UHS1(iPer) .le. testInten(iPer,iInten) ) then
      
           haz10 = exp( alog(UHS1(iPer) / testInten(iPer,iInten-1)) / 
     1                  alog( testInten(iPer,iInten)/ testInten(iPer,iInten-1))
     2                  * alog( haz_GMPE(iInten)/haz_GMPE(iInten-1) ) 
     3                  + alog(haz_GMPE(iInten-1)) )
           goto 20
          endif
         enddo
 20      continue

c        compute the deagg weights for this GMPE
         ratio = haz10 / hazLevel
         wt_deagg(jType,iAtten) = haz10 / hazLevel * gm_wt(iPer,jType,iAtten)

c        Write out the deagg wts
         write (50,'( f10.3,2i10,3f10.5)') period(iPer), jType,iAtten, 
     1          ratio, gm_wt(iPer,jType,iAtten), wt_deagg(jType,iAtten)
        enddo
       enddo
       
c      Compute the median and sigma for each GMPE for each scenario
       read (31,*,err=3000) nScenario
       read (31,'( a1)') dummy
       do iScen=1,nScenario

c        Write header for the scenario Med and Sigma         
         write (50,'( /,''   Scen    Mag     Rrup      Period    Med(g)    Sigma '' )') 

        read (31,*) jType, mag, rRup, rJB, Rx, rSeismo, Ry0, rHypo, dip,
     1        zTOR, theta_site, RupWidth, hypoDep, vs30, Z1, Z15, Z25,  
     2        vs30_class, HWFlag, ftype, forearc
     
       do iPer2=1,nProb    
         sum1 = 0.
         sum2 = 0.
         do iAtten=1,nGM_Model(iPer2,jType)

          AR = 1.
          
c         Compute the median and sigma for this GMPE
          call meanInten ( rRup, rJB, rSeismo,
     1           HWFlag, mag, jcalc(iPer2,jType,iAtten), period(iPer2),  
     2           lgInten,sigmaY, ftype, attenName, period1, 
     3           iAtten, iPer2, jType, vs30, hypoDep, intflag, AR, dip,
     4           rhypo, z1, z15, z25, tau,
     5           zTOR, theta_site, RupWidth, vs30_class, forearc, Rx, phi,
     6           cfcoefrrup, cfcoefrjb, Ry0 )


c         Add epistemic uncertainty term (constant shift) to median
          lgInten = lgInten + gmScale(iPer2,jType,iAtten)
          
c         Adjust the sigma
          sigmaY = sqrt(sigmaY*sigmaY + varadd(iPer2,jType,iAtten) )
          
c         Sum to get weighted median and sigma          
          sum1 = sum1 + lgInten * wt_deagg(jType,iAtten) 
          sum2 = sum2 + sigmaY * wt_deagg(jType,iAtten) 
         enddo

c        Write the deagg-weighted median and sigma for the scenario
         write (50,'( i5,5f10.4)') iScen, mag, rrup, period(iPer2), exp(sum1), sum2
         med(iPer2) = sum1
         sig(iPer2) = sum2
        enddo
        
c      Write header for the CMS_part1  
       write (50,'( /,''    Period    UHS       Med       Sigma     Tamp15    Eps*      Rho       Eps_bar   CMS '' )') 

c      Compute the CMS   
          call calcCMS ( med, sig, period, UHS1, epsilonstar, CMS,
     1                   iPer2, iPer, iScen, nScenario, nProb )
       enddo
       write (*,'( 2x,''CMS calc complete'')')
      

      close (50)

      write (*,*) 
      write (*,*) '*** CMS input Code (45) Completed with Normal Termination ***'
      stop

 2000 write (*,'( 2x,''bad CMS run file:'',a80)') filein
      stop
 2001 write (*,'( 2x,''Error in the hazard level in CMS run file'')')
      stop
 2002 write (*,'( 2x,''Error in the number of T* in CMS run file'')')
      stop
 2003 write (*,'( 2x,''Error in the T* indexes in CMS run file'')')
      stop
 3000 write (*,'( 2x,''Error in nScenario of input file'')')
      stop
 3001 write (*,'( 2x,''Error in scenario inputs for scenario number: '',i5)') iScen
      stop
      
      end

