c     Set array dimensions

      integer MAX_SITE, MAX_FLT, MAX_SEG, MAX_INTEN, MAX_PROB,
     1        MAX_DIP, MAXPARAM, MAX_MAG, MAX_DIST, MAX_EPS, MAX_N1, MAX_FIles, MAX_N2, MAX_XCosT,
     2        MAX_WIDTH, MAX_SAMPLE, MAX_RISK, MAXRUP, MAX_FIXED_MAG, MAX_MAGDIM, MAX_FTYPE,
     3        MAX_ATTEN, MAX_ATTENTYPE, MAX_BR, MAX_NODE
      integer MAX_Med, MAX_epistemic, MAX_sigma, MAX_SYN

      PARAMETER ( MAX_SITE  = 1, MAX_FLT = 160, MAX_SEG  = 10,
     1            MAX_INTEN = 30, MAX_PROB = 20, MAX_DIP=5,
     2            MAXPARAM = 200, MAX_MAG=15, MAX_DIST=15,
     3            MAX_EPS=15, MAX_N1=150, MAX_Files=3, MAX_N2=6,MAX_Xcost=10,
     4            MAX_WIDTH=12, MAX_SAMPLE=10000, MAX_RISK=4000)
      PARAMETER (MAXRUP=4, MAX_FIXED_MAG=3, MAX_MAGDIM=3, MAX_FTYPE=3, MAX_ATTEN=100,
     1            MAX_ATTENTYPE=4, MAX_BR=20, MAX_NODE=20)

      parameter (MAX_Med=20,MAX_epistemic=15, MAX_sigma=10, MAX_SYN=1)


