
      subroutine calcCMS ( med, sig, period, UHS1, epsilonstar, CMS, 
     1                     iPer2, iPer, iScen, nScenario, nProb )

      implicit none
      include 'cms.h'
      
      integer iScen, nScenario, iPer2, nProb, iFlag, iPer, i
      real period(MAX_PROB), med(MAX_PROB), sig(MAX_PROB), Tstarprime, 
     1     Tprime(MAX_PROB), T_max, T_min, c1, c2, c3, c4, rho(MAX_PROB),
     2     epsilon_bar(MAX_PROB), UHS2(MAX_PROB), period5(MAX_PROB), Tstar, 
     3     temp, epsilonstar, cms(MAX_PROB), pi, Tamp15, t1, Trockratio, 
     4     shape(MAX_PROB), UHS1(MAX_PROB)

      pi =3.1415926

c     Copy UHS to new array (this new one will be shifted)    
      do i=1,nProb
        UHS2(i) = UHS1(i)
        period5(i) = period(i)
      enddo

      do iScen=1,nScenario
        do iPer2 = 1,nProb
          shape(iPer2) = med(iPer2)-med(1)
        enddo

        t1 = alog(1.5)
          
        do iPer2 = 1,nProb-1
          if (shape(iPer2) .le. t1 .and. shape(iPer2+1) .ge. t1) then
            if (period5(iPer2) .eq. 0) then
               write (*,*) 'Re-run hazard with a PGA that is non-zero, interpolation error occurred'
	       stop 99
            endif
          call interp (exp(shape(iPer2)),exp(shape(iPer2+1)),
     1       alog(period5(iPer2)),alog(period5(iPer2+1)),exp(t1),Tamp15,iFlag)
          Tamp15=exp(Tamp15)
          endif
        enddo        
      enddo      
          
c     iPer is the index of tstar
c     Compute number of standard deviations (epsilon) to reach UHS at T*
      epsilonstar =  (log(UHS2(iPer))-med(iPer))/sig(iPer)

      do iPer2=1,nProb     
          
c         Calculate correlation coefficients based on Baker and Jayaram (2008)          
          Tstar = period5(iPer)
          Trockratio = Tamp15/0.1
          Tstarprime = Tstar/Trockratio

          Tprime(iPer2) = period5(iPer2)/Trockratio
          T_min = amin1(Tstarprime,Tprime(iper2))
          T_max = amax1(Tstarprime,Tprime(iper2))
          
          temp = 0.109
          if (T_min .gt. 0.109) temp = T_min
            C1 = 1-(cos((pi/2) - 0.366*alog(T_max/temp)))
                 
          if (T_max .lt. 0.2) then
            C2 = 1 - 0.105*(1 - 1./(1+exp(100*T_max-5)))*(T_max-T_min)/(T_max-0.0099)
          else
            C2 = 0
          endif
                   
          if (T_max .lt. 0.109) then
            C3 = C2
          else
            C3 = C1
          endif
                   
          C4 = C1 + 0.5 * (sqrt(C3) - C3) * (1 + cos(pi*(T_min)/(0.109)))
                   
          if (T_max .le. 0.109) then 
            rho(iPer2) = C2
          elseif (T_min .gt. 0.109) then
            rho(iPer2) = C1
          elseif (T_max .lt. 0.2) then
            rho(iPer2) = amin1(C2,C4)
          else
            rho(iPer2) = C4
          endif
          
          Epsilon_bar(iPer2) = rho(iPer2) * epsilonstar
          
          CMS(iPer2) = exp(med(iPer2))*exp(sig(iPer2)*epsilon_bar(iPer2))

      enddo
                      
       do iPer2=1,nProb     
        write (50,'(9f10.4)') Period5(iPer2), UHS2(iPer2), exp(med(iPer2)), 
     1                      sig(iPer2), Tamp15, epsilonstar, rho(iPer2), epsilon_bar(iPer2), CMS(iPer2)
       enddo
      return
      end
