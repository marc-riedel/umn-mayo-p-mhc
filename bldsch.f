c **************************************************************
c
c This file contains the subroutines:  bldsch bldbbh
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************

      subroutine bldsch(nml)

c .................................................
c PURPOSE: calculate coordinates for selected sidechains molecule 'nml'
c          
c OUTPUT:  xat,yat,zat,xbaat,ybaat,zbaat,xtoat,ytoat,
c          ztoat (via 'eyring')
c
c          1st backbone atom of 1st residue of 'nml': 
c
c          - it's position: from 'gbpr(1-3,nml)'
c          - it's axes: from 'setgbl' according to
c            global angles 'gbpr(4-5,nml)'
c
c CALLS: eyring, fnd3ba,setgbl,setsys
c .................................................

      include 'INCL.H'

      dimension xg(3),zg(3)


      call fnd3ba(nml,i1,i2,i3)
c ------------------------------ first 3 bb atoms of 'nml'
      ixrfpt(1,nml)=i1
      ixrfpt(2,nml)=i2
      ixrfpt(3,nml)=i3

c ------------------------------ position of 1st bb atom
      xat(i1) = gbpr(1,nml)
      yat(i1) = gbpr(2,nml)
      zat(i1) = gbpr(3,nml)

      rfpt(1,nml)=xat(i1)
      rfpt(2,nml)=yat(i1)
      rfpt(3,nml)=zat(i1)

      call setgbl(nml,i1,i2,i3, xg, zg)

      xbaat(i1) = zg(1)
      ybaat(i1) = zg(2)
      zbaat(i1) = zg(3)

      xtoat(i1) = xg(1)
      ytoat(i1) = xg(2)
      ztoat(i1) = xg(3)

      ifirs=irsml1(nml)
c      NEED SIDE CHAIN LIST

      do i=ifirs,irsml2(nml)
        jj=iatrs1(i)+7            ! start with cg
        ll=iatrs2(i)-2            ! end before c
        do j=jj,ll       
          jow=iowat(j)
          call eyring(j,jow)
        enddo
      enddo


      return
      end


c **************************************************************
c
c This file contains the subroutines:  bldbbh
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************

      subroutine bldbbh(nml)

c .................................................
c PURPOSE: calculate coordinates for backbone hydrogens of molecule 'nml'
c          
c OUTPUT:  xat,yat,zat,xbaat,ybaat,zbaat,xtoat,ytoat,
c          ztoat (via 'eyring')
c
c          1st backbone atom of 1st residue of 'nml': 
c
c          - it's position: from 'gbpr(1-3,nml)'
c          - it's axes: from 'setgbl' according to
c            global angles 'gbpr(4-5,nml)'
c
c CALLS: eyring, fnd3ba,setgbl,setsys
c .................................................

      include 'INCL.H'

      dimension xg(3),zg(3)

c       ---------------LAST AND FIRST RESIDUES HAVE TO BE FIXED 
      ifirs=irsml1(nml)+1
      print*,ifirs,irsml2(nml)
      do i=ifirs,irsml2(nml)
        ii=iatrs1(i)+1
        if(ityat(ii) .eq. 4)call eyring(ii,iowat(ii))
        do ii=iatrs1(i)+2,iatrs1(i)+6
          if(ityat(ii) .eq. 1)call eyring(ii,iowat(ii))
        enddo
c      print*,ii,ityat(ii)
      enddo
      do i=10,30
      write(6,'(5f8.2,i9)')xat(i),yat(i),zat(i)
      enddo

      return
      end
