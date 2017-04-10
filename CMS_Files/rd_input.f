c -------------------------------------------------------------------------

      subroutine RdInput (nInten,  testInten, nGM_model, nattentype,
     1 attenType, nProb, gm_wt, period, jcalc, gm_scale, VarAdd, version)

      include 'cms.h'

      real version      
      real testInten(MAX_PROB,MAX_INTEN)
      real segModelWt(1), probAct(1), minlat,maxlat,minlong,maxlong,maxdist
      integer nInten(MAX_PROB), ntotal, attentype(MAX_FLT)
      integer iii
      integer jCalc(MAX_PROB,4,MAX_ATTEN)
      real gm_Scale(MAX_PROB,4,MAX_ATTEN), varAdd(MAX_PROB,4,MAX_ATTEN)

      character*80 filein, title, fname(1), dummy
      integer nfiles, ix(MAX_FILES)
      integer nWidth(MAX_FLT)
      real period(MAX_PROB)

      integer nMagBins, nDistBins, nEpsBins, nXcostBins, soilampflag
      real magBins(MAX_MAG), distBins(MAX_DIST), epsBins(MAX_EPS)
      real XcostBins(MAX_XCOST)
      integer nProb, nattentype, nGM_Model(MAX_PROB,MAX_ATTENTYPE)
      real testwt, checkwt, c1, c2, gm_WT(MAX_PROB,4,MAX_ATTEN)

c     For CMS, all we need is the logic tree weights for the GMPEs
c     This will read past the inputs until it gets to the GMPE weights

c     Set Data file units
      nwr = 11

      ntotal = 0

c     Read in the number of data files.
c      read (5,*) nfiles
c     Program no longer allowed to read from multiple files.
      nFiles = 1

c     Loop over the number of files.
c      do 111 iii=1,nfiles

c     Open PSHA Run Input File
      read (31,'( a80)') filein
      write (*,'( a80)') filein
      open (20,file=filein,status='old',ERR=2000)

c     Open Input PSHA Source/Fault file
      read (20,'( a80)') filein
c      open (10,file=filein,status='old')

c     Check for version compatibility with hazard code
        read (20,*) version
         if (version .ne. 45.2 .and. version .ne. 45.1) then
         write (*,*) version
         write (*,*) 'Hazard faultfile format incompatible version of Haz45, use Haz45.2 or Haz45.1'
         stop 99
        endif

c     Read in parameters for background grid.
      read (20,*,ERR=2001) minlat,maxlat,minlong,maxlong

C     Added back in read of single maxdist
      read (20,*,ERR=2002) maxdist

c     Input Title (not used) 
      read(20,'( a80)') title

c     Number of Spectral Periods and Number of attenuation relations types
      read(20,*,ERR=2003) nProb, nattentype
c      write (*,'( 2i5)') nProb, nattentype
      
      do iprob=1,nProb

C       Read period, maxeps dir flag and gm intensities
        read (20,*,ERR=2004) specT, sigtrunc, dirflag 
        read (20,*,ERR=2005) nInten(iProb), (testInten(iProb,j), j=1,nInten(iProb))
        call CheckDim ( nInten(iProb), MAX_INTEN, 'MAX_INTEN' )
        period(iProb) = specT
x        write (*,'( i5,f10.3)') iProb, period(iProb)

C       Read in the suite of attenution models and wts for each attentype
        do j=1,nattentype
          checkwt = 0.0
          read (20,*,ERR=2006) nGM_model(iProb,j)
          
c         Check for Max number of attenuation model
          call checkDim ( nGM_model(iProb,j), MAX_ATTEN, 'MAX_ATTEN' )

          do jj=1,nGM_model(iProb,j)
            read (20,*,ERR=2007) jcalc(iProb,j,jj), c1, c2, gm_wt(iProb,j,jj), Varadd1, iMix
            gm_scale(iProb,j,jj) = c1 + c2
            varAdd(iProb,j,jj) = varAdd1
c            if ( jCalc(j,jj) .lt. 0 ) then
c               backspace (20)
c               read (20,*) jcalc(j,jj), c1, c2, wtgm(j,jj), Varadd, sCalc(j,jj), sigfix(j,jj), sssCalc(j,jj)
c            endif
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
