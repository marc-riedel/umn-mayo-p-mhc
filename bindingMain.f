!**************************************************************     
! Calculate weights for a multicanonical simulation
!
!
! Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
!                      Shura Hayryan, Chin-Ku Hu
! Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
!                      Jan H. Meinke, Sandipan Mohanty
!
! **************************************************************
! gfortran -c -O2 -g bindingMain.f 

      program bindingMain

      include "INCL.H"
      include "INCP.H"
      
!! Storage space for force-field library paths, seuence-file name, and variable-file name
      character*80 libdir, pdbfile, varfile, newname
      character*80 pdbnms(10)
!! Storage space for end groups
      character grpn*4,grpc*4
      integer iter, nsteps
      double precision eps
      
      pdbnms(10) = "pd10.pdb"
      pdbnms(1) = "pd1.pdb"
      pdbnms(2) = "pd2.pdb"
      pdbnms(3) = "pd3.pdb"
      pdbnms(4) = "pd4.pdb"
      pdbnms(5) = "pd5.pdb"
      pdbnms(6) = "pd6.pdb"
      pdbnms(7) = "pd7.pdb"
      pdbnms(8) = "pd8.pdb"
      pdbnms(9) = "pd9.pdb"


! === Input parts ie EXAMPLES/1bdd.pdb 1bddeq.pdb 10 15000 0.0000001

      write (*,'(a)') ' ENTER PDB file '
      read (*,'(a)',err=1) pdbfile 
  1   write (*,'(/,a,$)') ' file with SEQUENCE:'
      write (*,'(a)') ' ENTER OUTPUT file '
      read (*,'(a)') newname 
      write (*,'(a)') ' ENTER PDB file '
      read (*,*,err=3) iter, nsteps, eps  
  3   write (*,'(/,a,$)') ' BAD PARAMS'
  
! =================================================== Energy setup

!     Directory for SMMP libraries
!     Change the following directory path to where you want to put SMMP
!     libraries of residues. 
      libdir='SMMP/'

!!     Choose energy type with the following switch
!        0  => ECEPP2 or ECEPP3 depending on the value of sh2
!        1  => FLEX 
!        2  => Lund force field
!        3  => ECEPP with Abagyan corrections
!
      ientyp = 0

      sh2=.false.         ! .true. for ECEPP/2; .false. for ECEPP3
      epsd=.false.        ! .true. for  distance-dependent  dielectric
                          !  permittivity

      itysol= 0    !  0: vacuum
                   ! >0: numerical solvent energy
                   ! <0: analytical solvent energy & gradients

      call init_energy(libdir)

! ================================================= Structure setup

      grpn = 'nh2' ! N-terminal group
      grpc = 'cooh'! C-terminal group

      iabin = 0  ! =0: read from PDB-file
                 ! =1: ab Initio from sequence (& variables)

      varfile = ' '
      ntlml = 0
      print*,pdbfile,newname
      print*,vlvr(1:10)
c      call energy()
      print*,"end energy 1"
      call init_mhc(iabin,grpn,grpc,pdbfile,varfile)
c      call outpdb(1, "newname.pdb")
       
c  nmvr(i),nursvr(i), vlvr(i)*crd
!      iter = 10
!     nsteps = 15000
!      eps = 1.0d-7
      nml = 1
      do i=1,20 
        print*,nmvr(i),nursvr(i),vlvr(i)*crd,vlvr(i) ; enddo
      
c      print*,ivrml1(nml),ivrml1(nml)+nvrml(nml)-1
      do i = ivrml1(nml),ivrml1(nml)+nvrml(nml)-1
         if(nmvr(i).ne."phi".and.nmvr(i).ne."psi")
     #            print*,nmvr(i),nursvr(i),vlvr(i) 
      enddo

      do j3=1,10
c      print*,j
      do i = ivrml1(nml),ivrml1(nml)+nvrml(nml)-1
         if(nmvr(i).ne."phi".and.nmvr(i).ne."psi")then
           vlvr(i) = vlvr(i)+0.02
           if(vlvr(i) > pi)vlvr(i)=vlvr(i)-pi2
           iat=iatvr(i)
           ity=ityvr(i)
           if (ity.eq.3) then   ! torsion (assure phase=const.)
             iow=iowat(iat)   ! (cannot be 1st atom of 'nml')
             do j=1,nbdat(iow)
               jat=ibdat(j,iow)
               if (iowat(jat).eq.iow) then ! excl. ring
                 to=vlvr(i)
                 if (abs(to).gt.pi) to=to-sign(pi2,to)
                 toat(jat)=to
                 sntoat(jat)=sin(to)
                 cstoat(jat)=cos(to)
               endif
             enddo
           endif
           endif
      enddo
         call bldsch(1)
      call outpdb(nml, pdbnms(j3))
      enddo

      end program bindingMain

