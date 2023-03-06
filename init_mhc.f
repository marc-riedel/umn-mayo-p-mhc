c **************************************************************
c
c This file contains the subroutines: init_mhc
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku 
c Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
c                      Jan H. Meinke, Sandipan Mohanty
c
c **************************************************************
c FIXME: Data in varfile determines which molecule is changed.

      subroutine init_mhc(iabin,grpn,grpc,seqfile,varfile)

c ----------------------------------------------------------
c PURPOSE: construct starting structure of molecule(s)
c
c          iabin = 1  : ab Initio using sequence & 
c                       variables given in input files
c          iabin != 1 : sequence, variable information
c                       from PDB-file
c
c          grpn:        N-terminal group
c          grpc:        C-terminal group
c
c CALLS:   addend,bldmol,c_alfa,getmol,iendst, mklist, nursvr,
C          pdbread,pdbvars,redseq,redvar,setmvs
C
c ----------------------------------------------------------

      include 'INCL.H'
      include 'INCP.H'

cf2py character*80 optional, intent(in) :: seqfile = ' '
cf2py character*80 optional, intent(in) :: varfile = ' ' 
      
      character grpn*4,grpc*4
      character navr*3, nars*4  
      character seqfile*80, varfile*80
      integer ontlml
      logical readFromStdin

      ontlml = 1
      readFromStdin = .false.

      write (*,*) 'init_molecule: Solvent: ', itysol
         if (iendst(seqfile).le.1) then
 3          write (*,'(/,a)', ADVANCE='NO') ' PDB-file:'
            seqfil=' '
            read (*,'(a)',err=3) seqfil
         else
            seqfil = seqfile
         endif
         write (*,*) 'PDB structure ',seqfil(1:iendst(seqfil))
         print *, 'calling readpdb with ',seqfile
         call pdbread(seqfil,ier)
         print*,ntlml,nchp,nrsp,natp,nrs,nat 
         do k=1,natp 
           kp=ixatp(k)
           xat(kp)=xatp(k) ; yat(kp)=yatp(k) ; zat(kp)=zatp(k) 
         enddo 
c         do k=1,10 ; print*,k,xat(k),yat(k),zat(k) ; enddo 
c         do k=1,10 ; print*,k,nvr,vlvr(k),ityvr(k) ; enddo 
         call outpdb(1, "new0.pdb")
         
         if (ier.ne.0) stop

         call pdbvars()
c in pdbvars  getmol is called c OUTPUT:   molecule  - ivrml1,nvrml
c           residues  - iatrs1,ixatrs,iatrs2,ivrrs1,nvrrs
c           atoms     - nmat,ityat,cgat,blat,baat,csbaat,snbaat,
c                       toat,cstoat,sntoat
c           bonds     - nbdat,iowat,iyowat,ibdat(1-mxbd,),iybdat(1-mxbd,)
c                       ! 1st atom of 'nml': iowat indicates 1st bond
c                          to a FOLLOWING atom (not previous) !
c           variables - ityvr,iclvr,iatvr,nmvr
        
         do k=1,iatrs2(irsml2(1)) 
           kp=ixatp(k)
           if(kp .ne. 0)then
           xat(k)=xatp(kp) ; yat(k)=yatp(kp) ; zat(k)=zatp(kp) 
           endif
         enddo 
      do i=10,30
c      write(6,'(5f8.2,i9)')xat(i),yat(i),zat(i),xatp(i),yatp(i),ityat(i)
      enddo
c         print*,"pnnt24352",ntlml,nrsp,natp,ntlat,iatrs2(irsml2(1))
         call bldbbh(1)
         call bldsch(1)
      do i=10,30
c      write(6,'(5f8.2,i9)')xat(i),yat(i),zat(i),xatp(i),yatp(i),ityat(i)
      enddo
         call outpdb(1, "new2.pdb")
      return
      end
      
      
