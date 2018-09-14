c -------------------------------------------------------------------------

      subroutine RdInput (nInten,  testInten, nGM_model, nattentype,
     1           attenType, nProb, gm_wt, period, jcalc, gm_scale, varadd, 
     2           version, scalc, sigfix)

      implicit none
      include 'cms.h'

      integer nProb, nattentype, nGM_Model(MAX_PROB,MAX_ATTENTYPE),
     1        nInten(MAX_PROB), attentype(MAX_FLT), iprob, dirflag, j, 
     2        jj, iMix(MAX_PROB,4,MAX_ATTEN), jcalc(MAX_PROB,4,MAX_ATTEN),
     3        scalc(MAX_PROB,4,MAX_ATTEN) 
      real version, specT, sigtrunc, testInten(MAX_PROB,MAX_INTEN),
     1     minlat, maxlat, minlong, maxlong, maxdist, 
     2     gm_scale(MAX_PROB,4,MAX_ATTEN), varadd(MAX_PROB,4,MAX_ATTEN),
     3     checkwt, c1, c2, gm_wt(MAX_PROB,4,MAX_ATTEN), period(MAX_PROB),
     4     sigfix(MAX_PROB,4,MAX_ATTEN)
      character*80 filein, title

c     For CMS, all we need is the logic tree weights for the GMPEs
c     This will read past the inputs until it gets to the GMPE weights

c     Open PSHA Run Input File
      read (31,'( a80)') filein
      write (*,'( a80)') filein
      open (20,file=filein,status='old',ERR=2000)

c     Open Input PSHA Source/Fault file
      read (20,'( a80)') filein

c     Check for version compatibility with hazard code
        read (20,*) version
         if (version .ne. 45.3 .and. version .ne. 45.2 .and. version .ne. 45.1) then
         write (*,*) version
         write (*,*) 'Hazard faultfile format incompatible version of Haz45, use Haz45.3 or Haz45.2 or Haz45.1'
         stop 99
        endif

c     Read in parameters for background grid.
      read (20,*,ERR=2001) minlat, maxlat, minlong, maxlong

c     Added back in read of single maxdist
      read (20,*,ERR=2002) maxdist

c     Input Title (not used) 
      read(20,'( a80)') title

c     Number of Spectral Periods and Number of attenuation relations types
      read(20,*,ERR=2003) nProb, nattentype
      
      do iprob=1,nProb

c       Read period, maxeps dir flag and gm intensities
        read (20,*,ERR=2004) specT, sigtrunc, dirflag 
        read (20,*,ERR=2005) nInten(iProb), (testInten(iProb,j), j=1,nInten(iProb))
        call CheckDim ( nInten(iProb), MAX_INTEN, 'MAX_INTEN' )
        period(iProb) = specT

c       Read in the suite of attenuation models and wts for each attentype
        do j=1,nattentype
          checkwt = 0.0
          read (20,*,ERR=2006) nGM_model(iProb,j)
          
c         Check for Max number of attenuation model
          call checkDim ( nGM_model(iProb,j), MAX_ATTEN, 'MAX_ATTEN' )

          do jj=1,nGM_model(iProb,j)
            read (20,*,ERR=2007) jcalc(iProb,j,jj), c1, c2, gm_wt(iProb,j,jj), varadd(iProb,j,jj), iMix(iProb,j,jj)
            if ( jcalc(iProb,j,jj) .lt. 0 ) then
              backspace (20)
              read (20,*,ERR=2007) jcalc(iProb,j,jj), c1, c2, gm_wt(iProb,j,jj), varadd(iProb,j,jj), iMix(iProb,j,jj),
     1              sCalc(iProb,j,jj), sigfix(iProb,j,jj)
            endif
            gm_scale(iProb,j,jj) = c1 + c2
          enddo
        enddo

      enddo
      
      close (20)
       
      return
 2000 write (*,'( 2x,''Error: PHSA run file not found:'',a80)') filein
      stop
 2001 write (*,'( 2x,''Error in PHSA run file: min max lat long line '')')
      stop
 2002 write (*,'( 2x,''Error in PHSA run file: maxdist line '')')
      stop
 2003 write (*,'( 2x,''Error in PHSA run file: nProb,nattentype line'')')
      stop
 2004 write (*,'( 2x,''Error in PHSA run file: specT, sigtrunc, dirflag line'')')
      stop
 2005 write (*,'( 2x,''Error in PHSA run file: nInten, testInten line'')')
      stop
 2006 write (*,'( 2x,''Error in PHSA run file: nGM_model line'')')
      stop
 2007 write (*,'( 2x,''Error in PHSA run file: jcalc line'')')
      stop
     
      
      end
