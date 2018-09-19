
      subroutine Rd_Fault_Data (nFltTotal, nFlt0, f_start, f_num, AttenType, 
     1           n_Dip, n_bValue, nActRate, nSR, nRecInt, nMoRate, nMagRecur,  
     2           nThick1, nRefMag, nFtypeModels, dipWt, bValueWt, actRateWt, 
     3           wt_sr, wt_recInt, wt_MoRate, magRecurWt, faultThickWt, 
     4           refMagWt, ftmodelwt, nFtype, ftype_wt1, wt_RateMethod, al_Segwt, 
     5           nRate, rateType, nBR_SSC, nSegModel1, segwt1, segFlag, indexRate)

      implicit none
      include 'cms.h'
      
      integer nSegModel1(MAX_FLT), segFlag(MAX_FLT, MAX_SEG), iSeg, nFlt0, 
     1        nFltTotal, f_start(MAX_FLT), f_num(MAX_FLT), nSegModel(MAX_FLT), 
     2        faultflag(MAX_FLT,MAX_SEG,MAX_FLT), nBR_SSC(MAX_FLT,MAX_BR),
     3        n_Dip(MAX_FLT), n_bvalue(MAX_FLT), nActRate(MAX_FLT), iRecur, 
     4        nSR(MAX_FLT), nRecInt(MAX_FLT), nMoRate(MAX_FLT), Nrate(MAX_FLT),
     5        nMagRecur(MAX_FLT), nThick1(MAX_FLT), nRefMag(MAX_FLT,MAX_WIDTH), 
     6        nFtypeModels(MAX_FLT), Attentype(MAX_FLT), indexrate(MAX_FLT,4),
     7        nFtype(MAX_FLT,MAXPARAM), rateType(MAX_FLT,MAXPARAM), ii, iThick, 
     8        iCoor, nFLT2, iFlt, iFlt2, nsyn, nfp, iDepthModel, iOverRideMag,
     9        iFLt0, k, i, isourceType, insyn, synflag, directflag, ipt, nDownDip, 
     1        iThick1, iFM
      real segwt1 (MAX_FLT,MAX_SEG), al_segWt(MAX_FLT), segwt(MAX_FLT,MAX_FLT), 
     1     dipWt(MAX_FLT,MAXPARAM), bValueWt(MAX_FLT,MAXPARAM), actRateWt(MAX_FLT,MAXPARAM), 
     2     wt_SR(MAX_FLT,MAXPARAM), wt_RecInt(MAX_FLT,MAXPARAM), wt_MoRate(MAX_FLT,MAXPARAM), 
     3     magRecurWt(MAX_FLT,MAXPARAM), faultThickWt(MAX_FLT,MAXPARAM), probAct0, 
     4     refMagWt(MAX_FLT,MAX_Width,MAXPARAM), ftmodelwt(MAX_FLT,MAXPARAM),
     5     wt_rateMethod(MAX_FLT,4), ftype_wt1(MAX_FLt,MAXPARAM,5), ftype1(10,10),
     6     flat, flong, fz, sampleStep, minmag, sigArea, sigWidth, magsyn, dip1, 
     7     top, x(100), x2(100,100), bValue2, actRate 
      character*80 fName1, fName 

c     Input Fault Parameters
      read (10,*) iCoor
      read (10,*) NFLT0

      call CheckDim ( NFLT0, MAX_FLT, 'MAX_FLT' )
  
      iflt = 0

c     Loop over each fault in the source file....
      DO iFlt0=1,NFLT0
        read (10,'( a80)') fName1
        read (10,*) probAct0

c       Read number of segmentation models for this fault system       
        read (10,*) nSegModel(iFlt0)
        read (10,*) (segWt(iFlt0,k),k=1,nSegModel(iFlt0))

c       Read total number of fault segments for this fault system
        read (10,*) nFlt2
        do i=1, nSegModel(iFlt0)
         read (10,*) (faultFlag(iFlt0,i,k),k=1,nFlt2)
        enddo

c       Set the index for the first fault in this fault system and the number
c       This allows us to find the right fault from the large list
        f_start(iFlt0) = iFlt + 1
        f_num(iFlt0) = nFlt2 

c       Loop over number of individual fault segments....                        
        do iflt2=1,nflt2

          iFlt = iFlt + 1
          call CheckDim ( iflt, MAX_FLT, 'MAX_FLT   ' )

c         Set up segmentation cases and weights by individual faults
          nSegModel1(iFlt) = nSegModel(iFlt0)
          do iSeg=1,nSegModel1(iFlt)
            segwt1(iFlt,iSeg) = segWt(iFlt0,iSeg)
            segFlag(iFlt,iSeg) = faultFlag(iFlt0,iSeg,iflt2)
          enddo

c         Read past name of this segment
          read(10,'( a80)') fname

          read (10,*) isourceType, attenType(iFlt), sampleStep, directflag, synflag

c         Read past the synchronous Rupture parameters
          if (synflag .gt. 0) then
            read (10,*) nsyn
            do insyn=1,nsyn
              read (10,*) magsyn
            enddo
          endif

c         Read aleatory segmentation wts
          read (10,*) al_segWt(iFlt)

c         Check for standard fault source or areal source
          if ( isourceType .eq. 1 .or. isourceType .eq. 2) then
            read (10,*) dip1, top
            read(10,*) nfp     
            do ipt=1,nfp
              read (10,*) fLong, fLat
            enddo
          endif

c         Check for grid source (w/o depth)
          if ( isourceType .eq. 3 .or. isourceType .eq. 7 ) then
            read (10,*)  dip1, top
c           Read past grid filename...
            read (10,*)
          endif

c         Check for grid source (w/ depth)
          if ( isourceType .eq. 4 ) then
            read (10,*)  dip1
c           Read over grid filename...
            read (10,*)
          endif

c         Check for custom fault source
          if ( isourceType .eq. 5) then
            read(10,*) nDownDip, nfp
c           Only read in the first downdip case - rest is not needed...
            do ipt=1,nfp
              read (10,*) fLong, fLat, fZ
            enddo
          endif

c         Read dip Variation
          if ( isourceType .ne. 5 ) then
            read (10,*) n_Dip(iflt)
            nBR_SSC(iFLt,1) = n_Dip(iflt)
            read (10,*) (x(i),i=1,n_Dip(iflt))
            read (10,*) (dipWt(iFlt,i),i=1,n_Dip(iflt))
          else
            n_Dip(iflt) = 1
            dipWt(iFlt,1) = 1.
          endif

c         Read b-values (not for activity rate cases)
          read (10,*) n_bValue(iFlt)
          nBR_SSC(iFLt,11) = n_bValue(iFlt)
          if ( n_bValue(iFlt) .gt. 0 ) then
            read (10,*) (x(i),i=1,n_bValue(iFlt))
            read (10,*) (bValueWt(iFlt,i),i=1,n_bValue(iFlt))
          endif
                
c         Read activity rate - b-value pairs
          read (10,*) nActRate(iFlt)
          nBR_SSC(iFLt,8) = nActRate(iFlt)
          if ( nActRate(iFlt) .ne. 0 ) then
            do ii=1,nActRate(iFlt)
              read (10,*) bValue2, actRate, actRateWt(iFlt,ii)
            enddo
          endif

c         Read weights for rate methods
          read (10,*) (wt_RateMethod(iFlt,i), i=1,4)
          nBR_SSC(iFLt,6) = 4

          nRate(iFlt) = 0
c         Read slip-rates
          read (10,*) nSR(iFlt)
          nBR_SSC(iFLt,7) = nSR(iFlt)
          if ( nSR(iFlt) .gt. 0 ) then
            read (10,*) (x(k),k=1,nSR(iFlt))
            read (10,*) (wt_sr(iFlt,k),k=1,nSR(iFlt))
          endif
          do i=1,nSR(iFLt)
            nRate(iFlt) = nRate(iFlt) + 1
            RateType(iFLt,nRate(iFlt)) = 1
          enddo

          do i=1,nActRate(iFlt)
            nRate(iFlt) = nRate(iFlt) + 1
            RateType(iFLt,nRate(iFlt)) = 2
          enddo

c         Read recurrence intervals
          read (10,*) nRecInt(iFlt)
          nBR_SSC(iFLt,9) = nRecInt(iFlt)
          if ( nRecInt(iFlt) .gt. 0 ) then
            read (10,*) (x(k),k=1,nRecInt(iFlt))
            read (10,*) (wt_recInt(iFlt,k),k=1,nRecInt(iFlt))
          endif
          do i=1,nRecInt(iFLt)
            nRate(iFLt) = nRate(iFlt) + 1
            RateType(iFLt,nrate(iFLt)) = 3
          enddo

c         Read moment-rates
          read (10,*) nMoRate(iFlt)
          nBR_SSC(iFLt,10) = nMoRate(iFlt)
          if ( nMoRate(iFlt) .gt. 0 ) then
            read (10,*) (x(k),k=1,nMoRate(iFlt))
            read (10,*) (x(k),k=1,nMoRate(iFlt))
            read (10,*) (wt_MoRate(iFlt,k),k=1,nMoRate(iFlt))
          endif
          do i=1,nRecInt(iFLt)
            nRate(iFLt) = nRate(iFlt) + 1
            RateType(iFLt,nrate(iFLt)) = 4
          enddo
          indexRate(iFlt,1) = 0
          indexRate(iFlt,2) = indexRate(iFlt,1) + nSR(iFlt)
          indexrate(iFlt,3) = indexRate(iFlt,2) + nActRate(iFlt)
          indexrate(iFlt,4) = indexRate(iFlt,3) + nRecInt(iFlt)
                                   
c         Read Mag recurrence weights (char and exp)
          read (10,*) nMagRecur(iFlt)
          nBR_SSC(iFLt,4) = nMagRecur(iflt)
          read (10,*) (x(i),i=1,nMagRecur(iFlt))
          read (10,*) (magRecurWt(iFlt,i),i=1,nMagRecur(iFlt))

c         Read past corresponding magnitude parameters. 
          do iRecur=1,nMagRecur(iFlt)
            read (10,*) x(iRecur), x(iRecur), x(iRecur)
          enddo

c         Read seismogenic thickness
          if ( isourceType .ne. 5) then
            read (10,*) nThick1(iFlt)
            nBR_SSC(iFLt,2) = nThick1(iFlt)
            read (10,*) (x(i),i=1,nThick1(iFlt))
            read (10,*) (faultThickWt(iFlt,i),i=1,nThick1(iFlt))
          else
            nThick1(iFlt) = 1
            faultThickWt(iFlt,1) = 1.
          endif
         
c         Read depth pdf
          read (10,*) iDepthModel       

c         Read Mag method (scaling relations or set values)
          read (10,*) iOverRideMag

c         Read reference mags for each fault thickness
          iThick = 1
          nBR_SSC(iFLt,5) = 0 
          do iThick1=1,nThick1(iFlt)
            read (10,*) nRefMag(iFlt,iThick)
            read (10,*) (x2(iThick,i),i=1,nRefMag(iFlt,iThick))
            nBR_SSC(iFLt,5) = nBR_SSC(iFLt,5) + nRefMag(iFlt,iThick)
            read (10,*) (refMagWt(iFlt,iThick,i),i=1,nRefMag(iFlt,iThick))
            iThick = iThick + 1              
          enddo

c         Read Past remaining input for this fault
          read (10,*) minMag
          read (10,*) (x(k),k=1,2), sigArea
          read (10,*) (x(k),k=1,2), sigWidth

c         Read ftype Models
          read (10,*) nFtypeModels(iFlt)
          nBR_SSC(iFLt,3) = 0
          do iFM=1,nFtypeModels(iFlt)
            read (10,*) ftmodelwt(iFlt,iFM)
            read (10,*) nFtype(iFlt,iFM)
            read (10,*) (ftype1(iFM,k),k=1,nFtype(iFlt,iFM))
            read (10,*) (ftype_wt1(iFlt,iFM,k), k=1,nFtype(iFlt,iFM))
            nBR_SSC(iFLt,3) = nBR_SSC(iFLt,3) + nFtype(iFlt,iFM)
          enddo

c       End of Loop over iFlt2 - number of segments    
        enddo

c     End of Loop over iFlt
      enddo
      nFltTotal = iFlt

      return
      end

