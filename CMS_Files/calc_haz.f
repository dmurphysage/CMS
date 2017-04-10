      subroutine calcHaz ( haz, haz1, nBR_GMC, nFlt, iAtten, attenType, nInten, wt1 )

      implicit none
      include 'cms.h'
      integer iAtten(4,MAX_MED,MAX_epistemic,MAX_SIGMA)
      integer nBR_GMC(3), nFlt, attenTYpe(MAX_FLT), nInten, jtype
      integer iBR1, iBR2, iBR3, kFlt, i, jAtten
      real haz(MAX_INTEN,MAX_ATTEN, MAX_FLT)
      real*8 haz1(MAX_INTEN) 
      real wt1(4,MAX_ATTEN)
      
c     FIX LATER!
      jtype=1
      
      do i=1,nInten           
           haz1(i)= 0.
      enddo
      
      do iBR1=1,nBR_GMC(1)
       do iBR2=1,nBR_GMC(2)
        do iBR3=1,nBR_GMC(3)
         jAtten = iAtten(jType,iBR1,iBR2,iBR3)
         
c         write (*,*) jAtten, jType, iBR1, iBR2, iBR3
c         pause 'jatten, jType, iBR1, iBr2, iBR2'

         do kFlt=1,nFlt
          do i=1,nInten           
           haz1(i)= haz1(i) + haz(i,jAtten,kFlt)*wt1(AttenType(kFlt),jAtten)
C           if (i .eq. 1) then
C           write (*,*) iBr1, iBr2, iBr3, kFlt, jatten,i, haz(i,jAtten,kFlt), wt1(AttenType(kFlt),jAtten), haz1(i)
C           endif
          enddo
         enddo


        enddo
       enddo
      enddo 

C        pause 'hazcalc2'
      
      return
      end
