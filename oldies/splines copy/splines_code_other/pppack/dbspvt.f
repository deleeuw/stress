      subroutine DBSPVT ( t, k, index, x, left, biatx, dbiatx )
c  Derived from deBoor's DBSPVB by Eric Grosse, 4 May 1987, 23 Nov 1989
c  Calling sequence is unchanged except for addition of dbiatx,
c  which is defined by  dbiatx(i,m) =  d B(left-k+i) / d t(left-k+m)
c  for 1<=i<=k and 1<=m<=2*k.
c  For my own sanity, jhigh was renamed k.  It is still the spline order.
      parameter(jmax = 20)
      parameter(mmax = 40)
      integer index,k,left,i,j,jp1,m,mm
      DOUBLE PRECISION dbiatx(k,1), dsaved(mmax), dterm(mmax), denom
      DOUBLE PRECISION biatx(k),t(1),x,deltal(jmax),deltar(jmax),saved,t
     +erm
      data j/1/
      save j,deltal,deltar
                                        go to (10,20), index
   10 j = 1
      biatx(1) = 1.D0
      mm = 2*k
      do 1001 m=1,mm
1001     dbiatx(1,m) = 0.D0
      if (j .ge. k)                 go to 99
c
   20    jp1 = j + 1
         deltar(j) = t(left+j) - x
         deltal(j) = x - t(left+1-j)
         saved = 0.D0
         do 1002 m=1,mm
1002        dsaved(m) = 0.D0
         do 26 i=1,j
            denom = 1/(deltar(i) + deltal(jp1-i))
            term = biatx(i)*denom
            do 1003 m=1,mm
               dterm(m) = dbiatx(i,m)*denom
               if(i  .eq.m-k) dterm(m) = dterm(m) - term*denom
               if(i-j.eq.m-k) dterm(m) = dterm(m) + term*denom
1003           continue
            biatx(i) = saved + deltar(i)*term
            do 1004 m=1,mm
               dbiatx(i,m) = dsaved(m) + deltar(i)*dterm(m)
               if(i.eq.m-k) dbiatx(i,m) = dbiatx(i,m) + term
1004           continue
            saved = deltal(jp1-i)*term
            do 1005 m=1,mm
               dsaved(m) = deltal(jp1-i)*dterm(m)
               if(i-j.eq.m-k) dsaved(m) = dsaved(m) - term
1005           continue
   26       continue
         biatx(jp1) = saved
         do 1006 m=1,mm
1006        dbiatx(jp1,m) = dsaved(m)
         j = jp1
         if (j .lt. k)              go to 20
c
   99                                   return
      end

