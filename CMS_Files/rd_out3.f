      subroutine RdOut3 ( UHS1, hazLevel ) 
c      implicit none
      include 'cms.h'

      character*80 file1, dummy, fname(MAX_FLT)
      integer nRd, nFlt, nProb, nAmp(MAX_PROB)
      integer iProb, iAtten, iFlt, iAmp, k
      real lat, long, version
      real segModelwt(MAX_FLT),al_segwt(MAX_FLT),
     1     mindist(MAX_FLT),Haz(MAX_PROB,MAX_INTEN,MAX_FLT)
      real HazTotal(MAX_PROB,MAX_INTEN), amp(MAX_PROB,MAX_INTEN)
      real mBar(MAX_PROB,MAX_INTEN), dBar(MAX_PROB,MAX_INTEN),
     1  eBar(MAX_PROB,MAX_INTEN)
      real testHaz, x, UHS1(MAX_PROB)

c     open the out3 file
      read (31,'( a80)') file1 
      write (*,'( a80)') file1
c      pause 'in out3' 
      nRd = 23
      open (nRd,file=file1,status='old')

C     Check for version compatibility with hazard code
        read (nRd,*) version
         if (version .ne. 45.2) then
           write (*,*) 'out3 from incompatible version of Haz45, use Haz45.2'
           stop 99
         endif 
c     Read the output 3 file, keeping the total hazard and the background 
c     zone hazard 

c     Read in the hazard curves for all periods in given Haz45 ouptut file.
         read (nRD,*) nflt, testnum 

         do j=1,testnum
            read (nRD,'(a1)') dummy
         enddo

         do j=1,4
            read (nRD,'( a1)') dummy
         enddo

         read (nRD,'( 21x,2f9.3)') long, lat
         read (nRD,*) nProb

         do iProb=1,nProb

            read (nRD,'( 15x,i5)')  iAtten
            read (nRD,*) nAmp(iProb)
            read (nRD,'( 61x,30f12.4)') (amp(iProb,k),k=1,nAmp(iProb))
            
c            write(*,*) amp(iProb,1)
c            pause 'amp'

            do l=1,nFlt
               read (nRD,'( 2x,a38,2f6.3,f8.1,1x,30e12.4)') fname(l),
     1              segModelwt(l),al_segwt(l),
     1              mindist(l),(haz(iProb,k,l),k=1,nAmp(iProb))
            enddo

            read (nRD,'( 61x,50e12.4)') (haztotal(iProb,k),k=1,nAmp(iProb))
            read (nRD,'(a1)') dummy

            read (nRD,'( 61x,50e12.3)') (mBar(iProb,k),k=1,nAmp(iProb))
            read (nRD,'( 61x,50e12.3)') (dBar(iProb,k),k=1,nAmp(iProb))
            read (nRD,'( 61x,50e12.3)') (eBar(iProb,k),k=1,nAmp(iProb))

            do idum=1,3
               read (nRD,'( a80)') dummy
            enddo
            

         enddo
         close (nRD)

c     Interpolate to desired return period for the total hazard
      testHaz = hazLevel
      do iProb=1,nProb
	do iAmp=2,nAmp(iProb)

c         Check for zero values in hazard curve.
          if (hazTotal(iProb,iAmp) .eq. 0. ) then
            write (*,*) 'warning: Zero Values for hazard curve at desired haz level.'
            write (*,*) 'Setting UHS to last nonzero value'
            UHS1(iProb) = exp(x)
          endif
          
c         Interpolate the hazard curve.
          if ( hazTotal(iProb,iAmp) .lt. testHaz ) then
            x = ( alog(testHaz) - alog(hazTotal(iProb,iAmp-1)) )/
     1            ( alog(hazTotal(iProb,iAmp))-alog(hazTotal(iProb,iAmp-1))) 
     2          * (alog(amp(iProb,iAmp))-alog(amp(iProb,iAmp-1))) 
     3          + alog(amp(iProb,iAmp-1))
            UHS1(iProb) = exp(x)
            goto 10
          endif
        enddo
 10     continue
      enddo

      return
      end

    
