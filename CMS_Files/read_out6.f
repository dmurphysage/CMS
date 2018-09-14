c -------------------------------------------------------------------

      subroutine read_logichaz_out6 ( haz, nProb, nattenType, nAtten, nInten)
     
      implicit none
      include 'cms.h'

      real*8 haz(MAX_PROB, MAX_ATTENTYPE, MAX_ATTEN, MAX_INTEN)
      integer nInten, nProb, nattenType, nAtten(MAX_PROB,MAX_ATTENTYPE),
     1        iProb, jType, iAtten, iProb1, jType1, iAtten1, j, nwr
      real version
      character*80 file1
      
      nwr = 12
                     
c     Open output file
      read (31,'( a80)') file1
      write (*,'( 2x,a80)') file1

      write (*,*) 'Opening out6 file'
      write (*,*) file1
      open (nwr,file=file1,status='old',ERR=2000)

c     Check for version compatibility with hazard code
        read (nwr,*) version
         if (version .ne. 45.2) then
           write (*,*) 'out6 from incompatible version of Haz45, use Haz45.2'
           stop 99
         endif      

c     Currently, only works if nInten is the same for all periods
c     Fix this later

      do iProb=1,nProb
        do jType=1,nattenType
           do iAtten = 1,nAtten(iProb,jType)
               read (nwr,'( 3i5, 100e12.4 )',ERR=2001)  iProb1, jType1, iAtten1, 
     1             (Haz(iProb, jType, iAtten, j ), j=1,nInten)
           enddo
         enddo
      enddo
      close (nwr)

      return
 2000 write (*,'( 2x,''out6 file not found:'',a80)') file1
      stop
 2001 write (*,'( 2x,''Error in out6 file'')')
      write (*,'( 2x,''iProb,jType, iatten'',3i5)') iProb,jType, iatten
      stop
      end

c -------------------------------------------------------------------------
