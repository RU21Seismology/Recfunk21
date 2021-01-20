c  program rfmigrate_vl (for VERY LIGHT)
c
c 04/23/16 VLL adapted this code to computers not in JJP's office
c 1. plotit removed
c 2. eispack subroutines added at the end of the code so no need for libraries
c ALSO introduced a shift into the pick time so that picks of phases could be made at their _onset_
c and the pre-pick part of the data taper was estimated dynamically
c 11/15/19 VLL 
c	Paths for all ray parameter files (e.g. taup_P.dat) made "local", e.g. "./taup_P.dat"
c                   When running the code - copy or link these files to the work directory.
c
c
c
cc suggested compilation command:
c gfortran rfmigrate.f -o rfmigrate
c compiler version used: "GNU Fortran (GCC) 3.4.6"
c  ------
c  9/20/00 JJP -- adapted from recfunk  --- UPDATED 08/19/07 JJP
c   updated again 07/22/15 JJP
c  migrates MTC receiver functions in the frequency domain.
c  11/15/19 VLL _translation_ - in this context "migrate" means to correct for ray parameter
c
c
c  requires a stacking model in the anirec format
c  such a model may have anisotropy parameters, 
c  but migration code only uses the isotropic velocities. 
c
c  has kluge to cheat the pre-event noise for synthetic records  3/12/00 JJP
c  check to see if the kluge is commented out, because anirec_synth was updated 
c  in 2015 to have a proper-length pre-event time window.
c
c  this version of the RF code reads a file of data filenames 
c  you have two choices: either read the time intervals in the filename file
c  or read them in the sac header
c  the data must be binary SAC format
c  horizontals must be rotated to radial and transverse
c
c  for start times in the file:
c  the file defaults to "in_recfunk" and has lines of the form:
c
c  1997.070.19.33.bh?  <-- code replaces ? with r,t,z
c  57 52 	       <-- start time of analysis window, window duration (sec)
c  1997.076.08.15.bh?
c  62 62 
c  ...
c  ...
c  ...
c  stop                <-- tells code that data is finished, 799 events max
c
c
c  for start times in the SAC header, the default file is "in_recpick"
c  reads seismic record start times from the sac header
c  will search the SAC header for specific markers of P phases
c  T1 - P, Pdiff    ahead(12)
c  T2 - PKP,PKIKP   ahead(13)
c  T3 - PP          ahead(14)
c  T1=T2=T3=0 ==> use original A-marker ahead(9)
c
c  legacy data has timing in the
c   A-marker (ahead(9)) which doesnt indicate which phase was marked, 
c  For the A-marker, the phase is identified as P/Pdiff for DELTA<120, 
c  PKP/PKIKP for DELTA>120
c
c  code does NOT combine data with different sample rates
c  it is possible to spline interpolate the data files for a uniform sample rate
c  see the code in rfmig_mcboot.f, which should be spliceable into this code.
c  data files are limited to 99K pnts. To increase, see common block /datastuff/
c
c   many intermediate quantities are plotted with PLOTIT as the code proceeds.
c  other intermediate quantities can be plotted by uncommenting calls to PLOTIT
c
c  the code writes the BAZ- and EPICEN-dependent RFs to files
c  in a format easily digested by GMT (traces are separated by ">" lines)
c
c  filenames: out[rt]_baz.grid oum[rt]_epi.grid oum[rt]_baz.grid
c
c  out --> no moveout correction;  oum --> moveout correction
c
c  these files are overwritten everytime you run the program
c  so rename them to save them
c
c  xf77 -o /park/backup/bin/rfmigrate rfmigrate.f /park/backup/Plotxy/plotlib.a /park/backup/Ritz/eislib.a /park/backup/Ritz/jlib.a
c  xf77 -o /Users/jjpark/bin/rfmigrate rfmigrate.f /Users/jjpark/Plotxy/plotlib.a /Users/jjpark/Ritz/eislib.a /Users/jjpark/Ritz/jlib.a
c
      implicit real*4(a-h,o-z)
      implicit integer*4(i-n)
      real*8 ar,ai,el(12),aaa
      complex*16 zc,zero
      complex*16 afft,rf,rfs
c      real*8 sss,crf,drfs,drf,sig,sig2,bf
      real*8 drf
      character*120 name,subtext,string,name_stack
      character*28 title(2)
      character*10 cmps(3)
      character*18 xlabel(2)
      character*14 xlabel2(2)
      character*12 xlabel3(2),rfout(2),rfout1(2),rfout2(2)
      character comp(3)
      common/nnaamm/iwflag
      common/npiprol/anpi
      common/moovwin/nst(3)
      common/stap2/tap(16400,16)
      common/taperz/ta(16400,8,3),tai(16400,8,3)
      common/stap0/ar(16400),ai(16400),bb(16400),time(99000)
      common/specs/afft(4100,8,3),rf(4100,2),sss(4100,6),crf(4100,2)
      common/saveit/rft(4100,2999,2),crft(4100,2),rfs(4100,2999,2),
     x   drfs(4100,2999,2),bazs(2999),gcarc(2999),slow(2999),s2n(2999)
      common/datastuff/a(99000,3),pdat(16400),tdat(16400),drf(4100,2)
c NOTE THAT MODEL PARAMETERS ARE REAL*4 IN THIS CODE, NOT REAL*8 AS IN ANIREC
c ONLY 20 LAYERS OVER HALFSPACE ARE ALLOWED
      common/model/z(31),dz(31),rho(31),vp(31),vs(31),vp2(31),
     x               vp4(31),vs2(31),vss(31)
      common/model2/xmu(31),xla(31),xmu2(31),xla2(31),xla4(31)
      common/modstk/tau(31),itau(31),psd(31),gam(31),xik(31)
      common/migstk/arm(4100,2,31),aim(4100,2,31),dam(4100,2,31)
      common/tt_p/d_p(200),d_pkp(200),s_p(200),s_pkp(200)
      common/tt_pp/d_pp(200),d_pkkp(200),s_pp(200),s_pkkp(200)
      common/header/ahead(158)
      common/distrib/xi(72,36),yi(72,36),sumij(72,36)
      common/chisq/chi(160000),enchi(160000),bbaz(160000),
     x chim(200,2),enchim(200,2)
      dimension iah(158),chead(158),mmonth(12),fr(3),tt(3),ttt(3),
     x  xx(2),yy(2),jbaz(180),jjrec(180)
      equivalence (ahead,iah),(ahead,chead)
      data mmonth/0,31,59,90,120,151,181,212,243,273,304,334/
      data comp/'r','t','z'/
c      data cmps/'Radial    ','Transverse','Vertical  '/
c      data title/'Radial Receiver Function    ',
c     x 'Transverse Receiver Function'/
cc      data xlabel/'H\\sub{R}(f) phase','H\\sub{T}(f) phase'/
c      data xlabel2/'|H\\sub{R}(f)|','|H\\sub{T}(f)|'/
c      data xlabel3/'H\\sub{R}(t)','H\\sub{T}(t)'/
      data rfout/'rfout_radial','rfout_transv'/
      data rfout1/'rfoua_radial','rfoua_transv'/
      data rfout2/'rfoup_radial','rfoup_transv'/
      pi=3.14159265358979
      con=180./pi
      pih=pi/2.
      zero=cmplx(0.,0.)
      fmax=5.
      print *,'enter velocity model for stacking', 
     x ' <space-ret>=stack_model'
      read(5,101) name_stack
      if(name_stack(1:1).eq.' ') then
        open(7,file='stack_model',form='formatted')
      else
        open(7,file=name_stack,form='formatted')
      endif
      read(7,101) subtext
      print *,subtext
      read(7,*) nl
      if(nl.gt.20) then
        print *,'only twenty layers in stack model are allowed!'
        stop
      endif
c  read in theta,phi in degrees - polar coords of fast axis
      nlp=nl+1
      nlm=nl-1
      do i=1,nlp
        read(7,*) theta,phi
c  read depth to ith interface, vp (m/sec), pk-to-pk cos(2th) relative P pert
c  pk-to-pk cos(4th) relative P pert, v_s,  pk-to-pk cos(2th) relative S pert
c  density (kg/m**3)
        read(7,*) z(i),vp(i),vp2(i),vp4(i),vs(i),vs2(i),rho(i)
c  recall that we interpret fractional values of b,c,e 
c  as peak-to-peak relative velocity perts.
c  therefore, e=0.02 is 2% pert to mu from slowest to fastest
c  DO NOT NORMALIZE THE MODEL
        xmu(i)=rho(i)*vs(i)**2
        xmu2(i)=vs2(i)*xmu(i)
        xla(i)=rho(i)*vp(i)**2
        xla2(i)=vp2(i)*xla(i)
        xla4(i)=vp4(i)*xla(i)
        vs(i)=vs(i)
        vp(i)=vp(i)
        rho(i)=rho(i)
        z(i)=z(i)
      end do
      close(7)
      do i=2,nl
        dz(i)=z(i)-z(i-1)
      end do
      dz(1)=z(1)
c
c     11/16/19 VLL insert
c
      sdelay=0.
      pdelay=0.
      do i=1,nl
        sdelay=sdelay+(dz(i)/vs(i))
        pdelay=pdelay+(dz(i)/vp(i))
        tau(i)=sdelay-pdelay
	print 112,z(i),vp(i),vs(i),tau(i)
      end do
  112 format(4f11.2)
      print *, 'organ-pipe mode count at 1 Hz in stack of layers: S & P'
      print *,'These are the S & P traveltimes from the basal interface'
      print 104,sdelay,pdelay
      print *,'enter the depth at which you wish to focus the RFs (km)'
      read(5,*) itarget
      targm=itarget*1000.
      do i=nl,1,-1
        resz=targm-z(i)
	vdelay=tau(i)
        if(targm.gt.z(i)) go to 313
      end do
      i=0
      resz=targm
      vdelay=0.
  313 idelay=i+1
      vdelay=vdelay+resz/vs(idelay)-resz/vp(idelay)
      print *,'delay at vertical incidence = ',vdelay,' sec'
c
c end insert 11/16/19
c
      call read_ttimes
      print *,'enter fmax (Hz)'
      read(5,*) fmax
c  some default values
      anpi=2.5
      nwin=3
      npad=16384
      nnf=npad/2
      do i=1,4100
        crft(i,1)=0.
        crft(i,2)=0.
      end do
      irec=0
 1010 print *,'time intervals in file (0) or in SAC-header (1)'
      read(5,*) isegment
      if(isegment.eq.0) then
        print *,
     x'file of seismic records to read? (space-ret: in_recfunk)'
        read(5,101) name
        if(name(1:1).eq.' ') then
          open(10,file='in_recfunk',form='formatted')
        else
          open(10,file=name,form='formatted')
        endif
      elseif(isegment.eq.1) then
        print *,
     x'file of seismic records to read? (space-ret: in_recpick)'
        read(5,101) name
        if(name(1:1).eq.' ') then
          open(10,file='in_recpick',form='formatted')
        else
          open(10,file=name,form='formatted')
        endif
	print *,'duration of data windows for analysis'
	read(5,*) postwin
      else
        go to 1010
      endif
      print *,'minimum number of events in a binsum (eg 1 or 2)'
      read(5,*) ibskip
      print *,'rotate to LQT coordinates to isolate P and SV? (0=no)'
c     read(5,*) lqt
c
c     11/16/19 VLL this is a piece from another code, changes near-surface velocity
      print *,'enter -1 to change reference crust velocity 7.5 km/sec' 
      read(5,*) lqt
      refvel=7.5
      if(lqt.lt.0) then
        print *,'refvel ='
	read(5,*) refvel
      endif
c
c     11/16/19 VLL start reading data 
c
   10 print *,'enter filename e.g. 1998.361.bhz'
      read(10,101) name
  101 format(a)
      if(name(1:4).eq.'stop') go to 111
c  we assume that the last character is either 'r','t',or 'z'
      do i=80,1,-1
        if(name(i:i).ne.' ') then
          namlen=i
          go to 98
        endif
      end do
      print *,name
      print *,'no filename?'
      stop
   98 continue
      do kcmp=1,3
        name(namlen:namlen)=comp(kcmp)
        print *,name
        call sacread(name,a(1,kcmp),ierr)
        dt=ahead(1) 
        tt(1)=0.
        tt(2)=dt
        tt(3)=dt
        nscan=iah(80)
        if(nscan.gt.99000) then
          print *,'CAREFUL! Data length is greater than 99000:', nscan
          call pause('return to return')
        endif
      end do
c  11/21/19. VLL modifications to make code more universal
      irec=irec+1
      if(irec.gt.2999) then
        print *,'max records exceeded',irec
	stop
      endif
c 11/21/19 NOTE that durations of output RF are hard-wired here
      npre=15./dt
      npost=25./dt
      ntot=npre+npost
      ttt(1)=-npre*dt
      ttt(2)=dt
      ttt(3)=dt
      do i=1,ntot
        time(i)=-ttt(1)-(i-1)*ttt(2)
      end do
c
c
      baz=ahead(53)
      bazs(irec)=ahead(53)
      gcarc(irec)=ahead(54)
      print *,bazs(irec),'= back azim ',ahead(54),'= epicentral dist'
      if(ahead(53).lt.-0.1.or.ahead(53).gt.360.1) then
        print *,'HEY! THE BACK AZIMUTH OF THIS RECORD IS OUT OF RANGE!'
        stop
      elseif(ahead(54).lt.-0.1.or.ahead(54).gt.180.1) then
        print *,'HEY! THE GCARC OF THIS RECORD IS OUT OF RANGE!'
        stop
      endif
c  get p-slowness -- imark branches to determine the slowness rule
      igcarc=gcarc(irec)
      fac=gcarc(irec)-igcarc
      print *,'igcarc,i,gcarc(i),fac',igcarc,irec,gcarc(irec),fac
      if(isegment.eq.0) then
        print *,'P window start time, duration (sec)'
        read(10,*) pstart,postwin
	imark=0
      else
c   READING TIME INTERVAL FROM SAC HEADER
c  ahead(12) is SAC's T1 parameter, for arrival time.  Set by pressing T & 1 in plotpk
c   READING TIME INTERVAL FROM SAC HEADER
c  ahead(12) is SAC's T1 parameter, for arrival time.  Set by pressing T & 1 in plotpk
c
c VL 10/2016 a fix to the incorrectly (precisely on time) picked P/PKP phases
c simple - shift picked window back by 1/5 of the analysis interval, postwin/5
c
        if(ahead(12).gt.1) then
          tmark=ahead(12)-postwin/5
          imark=1
          slow(irec)=(1.-fac)*s_p(igcarc)+fac*s_p(igcarc+1)
c  ahead(13) is SAC's T2 parameter, for arrival time.  Set by pressing T & 2 in plotpk
        elseif(ahead(13).gt.1) then
          tmark=ahead(13)-postwin/5
	  imark=2
          slow(irec)=(1.-fac)*s_pkp(igcarc)+fac*s_pkp(igcarc+1)
c  ahead(14) is SAC's T3 parameter, for arrival time.  Set by pressing T & 3 in plotpk
        elseif(ahead(14).gt.1) then
          tmark=ahead(14)-postwin/5
          slow(irec)=(1.-fac)*s_pp(igcarc)+fac*s_pp(igcarc+1)
	  imark=3
c  ahead(9) is SAC's A parameter, for arrival time.  Set by pressing A in plotpk
        elseif(ahead(9).gt.1) then
          tmark=ahead(9)
	  imark=0
        else
          print *,' hey! no P-wave mark for file ',name
          stop
        endif
        pstart=tmark-ahead(6)-3.
      endif
      if(imark.eq.0) then
        if(igcarc.lt.120) then		! P-phase
          slow(irec)=(1.-fac)*s_p(igcarc)+fac*s_p(igcarc+1)
        else  ! PKP-phase
          slow(irec)=(1.-fac)*s_pkp(igcarc)+fac*s_pkp(igcarc+1)
	endif
      endif
	
c  ddf is cycles/day
      ddf=1./(dt*npad)
      nf=fmax/ddf+1
      if(nf.gt.4100) then
        print *,'nf is too big', nf
	stop
      endif
      fr(1)=0.
      fr(2)=ddf
      fr(3)=ddf
c
c   11/16/19 VLL insert code to shift the horizontal time windows so that RFs target a chosen depth
c
c  the stack-model velocities are in SI units, s_p is in sec/km
      slw=(slow(irec)/1000.)       
c  compute the slowness-dependent time delays at the interfaces of stack-model
c  if slowness of the P wave is too large to imply raybottom beneath target,
c  we skip the record
      if(slw*vp(idelay).ge.1.0) then
        irec=irec-1
	go to 10
      endif
      tshft=0.
      if(idelay.ge.1) then
        do k=1,idelay-1
          coss=sqrt(1.-vs(k)*vs(k)*slw*slw)
          cosp=sqrt(1.-vp(k)*vp(k)*slw*slw)
c  trapdoor to avoid incoming slownesses that are evanescent (e.g. head waves)
c  basically, one stacks with the last acceptable interval velocity
          tshft=tshft+dz(k)*(coss/vs(k)-cosp/vp(k))
        end do
      endif
      coss=sqrt(1.-vs(idelay)*vs(idelay)*slw*slw)
      cosp=sqrt(1.-vp(idelay)*vp(idelay)*slw*slw)
      tshft=tshft+resz*(coss/vs(idelay)-cosp/vp(idelay))
c
c.  11/16/19 VLL end insert 
c
c
      npts=postwin/dt
      nst(3)=pstart/dt
c
      nst(1)=(pstart+tshft)/dt
      nst(2)=(pstart+tshft)/dt
c
      call taper(npts,nwin,el)
      print *,'eigenvalues ',(el(k),k=1,nwin)
c  loop over components
c      print *,'npad,nf,nst,npts',npad,nf,nst,npts
      print *,'npad,nf,nst,npts',npad,nf,(nst(i),i=1,3),npts
      print *,'dt,nscan,ddf',dt,nscan,ddf
c                call pause('return to return')
c
      npts0=min0(nst(3)+1,npts)
c
      if(npts0.lt.npts) then
        call taper(npts0,nwin,el)
        print *,'eigenvalues ',(el(k),k=1,nwin)
      endif
c  vertical/radial rotation determined by slowness of the wave
c  slow=sin(\xi)/PVEL, \xi=asin(slow*PVEL), 
c  assume that PVEL=6.0 km/sec at surface
      zrrot=asin(slow(irec)*refvel)
c  here is the option to rotate the Radial and vertical Components to
c  a P and SV coordinate system    
      if(lqt.ne.0) then
        cs=cos(zrrot)
        sn=sin(zrrot)
c  rotate from Z-R coords to L-Q, similar to P-SV
c        do i=1,nst+npts
c          plog=cs*a(i,3)+sn*a(i,1)
c	  a(i,1)=cs*a(i,1)-sn*a(i,3)
c   	  a(i,3)=plog
c        end do
c
c.   11/16/19.  VLL
c
        do i=1,npts
          bb(i)=cs*a(nst(3)+i,3)+sn*a(nst(3)+i,1)
	  bb(i+npts)=cs*a(nst(1)+i,1)-sn*a(nst(1)+i,3)
        end do
        do i=1,npts
	  a(i+nst(3),3)=bb(i)
	  a(i+nst(1),1)=bb(i+npts)
        end do
c                11/16/19 end insert
      endif
      do kcmp=1,3
c  first: compute spectrum of pre-event noise
        do i=1,nf
          sss(i,kcmp)=0.
        end do
        do k=1,nwin
          do i=1,npts0
c OK to reverse order of data here, as we only use mod-sq spectrum
            ar(i)=a(nst(3)+1-i,kcmp)*tap(i,k)
            ai(i)=0.d0
          end do
          do i=npts0+1,npad
            ar(i)=0.d0
            ai(i)=0.d0
          end do   
c  we use the complex FFT routine fft2 -- uses real*8 input arrays
c      print *,'we use the complex FFT routine fft2'
          call fft2(ar,ai,npad)
          do i=1,nf
            sss(i,kcmp)=sss(i,kcmp)+(ar(i)**2+ai(i)**2)
          end do
        end do
c  if pre-event noise duration is truncated, boost the noise spectrum by
c  relative length factor
        if(npts0.lt.npts) then
          do i=1,nf
            sss(i,kcmp)=sss(i,kcmp)*float(npts)/float(npts0)
          end do
        endif
      end do
      if(npts0.lt.npts) then
        call taper(npts,nwin,el)
        print *,'eigenvalues ',(el(k),k=1,nwin)
      endif
      do kcmp=1,3
        do i=1,nf
          sss(i,kcmp+3)=0.
        end do
c  first: compute spectrum of P coda
        do k=1,nwin
          do i=1,npts

c  11/16/19 VLL changes value of widow start point
c
            ar(i)=a(i+nst(kcmp),kcmp)*tap(i,k)
            bb(i)=a(i+nst(kcmp),kcmp)*tap(i,k)
c   end change
c
            ai(i)=0.d0
          end do
c          call plotit(tt,bb,dum,npts,'tapered data',
c     x     'time',' ',1,0,0.0,0,1)
          do i=npts+1,npad
            ar(i)=0.d0
            ai(i)=0.d0
          end do
c  we use the complex FFT routine fft2 -- uses real*8 input arrays
c      print *,'we use the complex FFT routine fft2'
          call fft2(ar,ai,npad)
          do i=1,nf
            afft(i,k,kcmp)=dcmplx(ar(i),ai(i))
            sss(i,kcmp+3)=sss(i,kcmp+3)+(ar(i)**2+ai(i)**2)
          end do
        end do
      end do
c  simple correlation coefficient with vertical (kcmp=3)
      do l=1,2
        do n=1,nf
          zc=zero
          do k=1,nwin
            zc=zc+conjg(afft(n,k,3))*afft(n,k,l)
          end do
          rf(n,l)=zc/(sss(n,6)+sss(n,3))
          crf(n,l)=zabs(zc)**2/(sss(n,6)*sss(n,3+l))
          drf(n,l)=dsqrt((1.d0-crf(n,l))/(crf(n,l)*(nwin-1)))
     x                                         *zabs(rf(n,l))
          crft(n,l)=crft(n,l)+crf(n,l)
c          print *,n,l,zc,aaa,rf(n,l)
          rfs(n,irec,l)=rf(n,l)
          drfs(n,irec,l)=drf(n,l)
        end do
      end do
      sq2=sqrt(2.)
c      if(isss.ne.12345) go to 787
c      print *,isss
c                call pause('return to return')
      do l=1,2
        do n=1,nf
c	  if(n.ge.nf/2) then
c            fac=cos(pi*(n-nf/2)/nf)**2
c	  else
c            fac=1.0
c	  endif
          fac=cos(pi*(n-1)/(2*nf))**2
          ar(n)=real(rf(n,l))*fac
          ar(npad+2-n)=real(rf(n,l))*fac
          ai(n)=-imag(rf(n,l))*fac
          ai(npad+2-n)=imag(rf(n,l))*fac
          pdat(n)=drf(n,l)*fac/sq2
          tdat(n)=fr(1)+(n-1)*fr(2)
        end do
        do n=nf+1,npad/2+2
          ar(n)=0.d0
          ar(npad+2-n)=0.d0
          ai(n)=0.d0
          ai(npad+2-n)=0.d0
        end do
        dtdat=(tdat(nf)-tdat(1))/30.
        tdata=tdat(1)-dtdat
        tdatb=tdat(nf)+dtdat
        bbmx=0.
        do n=1,npad
          bb(n)=dsqrt(ar(n)**2+ai(n)**2)
          bbmx=amax1(bbmx,bb(n))
        end do
c       call plotit_axes(0.,0.,0.,0.)
        do n=1,npad
          if(bb(n).gt.1e-4) then
            pdat(n)=con*(pdat(n)/bb(n)) 
          else
            pdat(n)=360.
          endif
          bb(n)=con*datan2(ai(n),ar(n))
        end do
c       call plotit_axes(0.,0.,0.,0.)
c  invert the Fourier transform
        call fft2(ar,ai,npad)
c  normalization factor:
c  divide by npad for fft2 routine normalization
c  mult by nnf/nf to compensate for losing high freqs
c  mult by 2 if cosine**2 taper is applied to spectral RF
        fac=2.*float(nnf)/(float(npad)*float(nf))
c        fac=(4./3.)*float(nnf)/(float(npad)*float(nf))
c        fac=float(nnf)/(float(npad)*float(nf))
        print *,'ddf,npad,nnf,nf,fac',ddf,npad,nnf,nf,fac
c
c       11/17/19 this looks like the length of output
        nmm=30./dt
c
c     11/17/19. VLL damn it, hard-wired pre-event time....
c     changing 200 to 1000
c     changing to flexible definition
c        nnm=nmm+200
c	   nnm=nmm+1000
c        do i=1,nmm
c          bb(1000+i)=ar(i)*fac
c          time(1000+i)=(i-1)*dt
c        end do
c        do i=1,1000
c          bb(1001-i)=ar(npad+1-i)*fac
c          time(1001-i)=-i*dt
c        end do
c        ttt(1)=-1000*dt
c        ttt(2)=dt
c        ttt(3)=dt
c
         do i=1,npost
          bb(npre+i)=ar(i)*fac
          time(npre+i)=(i-1)*dt
        end do
        do i=1,npre
          bb(npre+1-i)=ar(npad+1-i)*fac
          time(npre+1-i)=-i*dt
        end do
        ttt(1)=-npre*dt
        ttt(2)=dt
        ttt(3)=dt

c        11/17/19 end of changes
c
c        do i=1,nnm 
	   do i=1,ntot
          rft(i,irec,l)=baz+50.*bb(i)
        end do
c        call plotit(ttt,bb,dum,400,'magnified version',
c     x    'time(sec)','H(t)',1,0,0.0,0,59+3*l)
        nmpt=40./dt
c        call plotit(ttt,bb,dum,nnm,title(l),
c     x    'time(sec)',xlabel3(l),1,0,0.0,0,60+3*l)
c        open(12,file=rfout(l),form='formatted')
c        write(12,1077) (time(i),bb(i),i=1,nnm)
c        close(12)
      end do
  787 continue
 1076 format(3f10.4)
 1077 format(2f10.4)
      go to 10
  111 continue
      close(10)
c  plot the inferred angle between radial and vertical against epicentral dist
c      do i=1,irec
c        if(gcarc(i).lt.80.) then
c	  bb(i)=42.+(20.-gcarc(i))*(20./60.)
c	elseif(gcarc(i).lt.118.) then
c	  bb(i)=16.+(6./900.)*(110.-gcarc(i))**2
c	else
c	  bb(i)=9.0
c	endif
c      end do
c      call plotit(gcarc,bb,dum,irec,' ',' ',' ',2,0,0.0,0,0)
c      call plotit(gcarc,bbaz,dum,irec,'P-SV Rotation Check',
c     x 'Epicentral Distance','Inferred R/Z Angle',2,0,0.1,3,1)
c  average coherence as function of frequency
      do l=1,2
        do i=1,nf
          crft(i,l)=crft(i,l)/irec
        end do
      end do
      print *,irec
c      call plotit_axes(tdata,tdatb,0.0,1.0)
      bb(1)=0.
      bb(2)=fmax
      bb(3)=1./float(nwin)
      bb(4)=bb(3)
c      call plotit(bb,bb(3),dum,2,' ',' ',' ',2,0,0.0,0,0)
c      call plotit(fr,crft(1,2),dum,nf,'Average P-Coda Coherence',
c    x    'Freq(Hz)','Avg C\\sup{2}(f)',1,0,0.05,0,0)
c      call plotit(fr,crft(1,1),dum,nf,'Average P-Coda Coherence',
c     x    'Freq(Hz)','Avg C\\sup{2}(f)',1,0,0.0,0,1)
      print *,irec
      do i=1,nmpt
        time(i)=-ttt(1)-(i-1)*ttt(2)
      end do
      nnn=25./dt
      print *,irec
      kaz=0
      open(12,file='outr_baz.grid',form='formatted')
      open(13,file='outt_baz.grid',form='formatted')
      ibazmax=355
      ibazmin=0
      ibazinc=-5
      print *,'back-azimuth range to plot:'
      print *,ibazmax,ibazmin,ibazinc
      print *,'change baz-interval or baz-spacing? (1=yes)'
      read(5,*) ickbaz
      if(ickbaz.eq.1) then
        print *,'enter back-azimuth range:'
        print *,'ibazmax, ibazmin, ibazinc'
        print *,'will make bazinc negative to step back thru baz values'
        read (5,*) ibazmax, ibazmin, ibazinc
        if((ibazmax-ibazmin)*ibazinc.gt.0) ibazinc=-ibazinc
      endif
      abazinc=abs(ibazinc)
      do ibaz=ibazmax,ibazmin,ibazinc
        baz=ibaz
        do l=1,2
          do n=1,nf
            rf(n,l)=zero
            drf(n,l)=0.
          end do    
          jrec=0    
          do i=1,irec
            test=abs(baz-bazs(i))
            if(test.lt.abazinc.or.test.gt.360.-abazinc) then
              jrec=jrec+1
              do n=1,nf
                sig=drfs(n,i,l)
                rf(n,l)=rf(n,l)+rfs(n,i,l)/sig
                drf(n,l)=drf(n,l)+1.0/sig
              end do  
            endif
          end do
c          print *,'# of records in bin:',jrec
          if(jrec.ge.ibskip) then
            kaz=kaz+1
c  do the expectation of variance for the weighted sum yourself
c  if you doubt the formula here -- basically the weighted terms in rf-sum
c  have unit variance, so that variance of the total variance is jrec/drf**2
            do n=1,nf
c	      if(n.ge.nf/2) then
c                fac=cos(pi*(n-nf/2)/nf)**2
c	      else
c                fac=1.0
c	      endif
              fac=cos(pi*(n-1)/(2*nf))**2
              rf(n,l)=rf(n,l)/drf(n,l)
              ar(n)=fac*real(rf(n,l))
              ai(n)=-fac*imag(rf(n,l))
              ar(npad+2-n)=fac*real(rf(n,l))
              ai(npad+2-n)=fac*imag(rf(n,l))
              drf(n,l)=fac*sqrt(float(jrec))/drf(n,l)
            end do
            do n=nf+1,npad-nf+1
              ar(n)=0.d0
              ai(n)=0.d0
            end do
            write(string,109) baz
            bbmx=0.
            do n=1,nf
              bb(n)=dsqrt(ar(n)**2+ai(n)**2)
              bbmx=amax1(bbmx,bb(n))
              pdat(n)=drf(n,l)/sq2
            end do
            do n=1,npad
              if(bb(n).gt.1e-4) then
                pdat(n)=con*(pdat(n)/bb(n)) 
              else
                pdat(n)=360.
              endif
              bb(n)=con*datan2(ai(n),ar(n))
            end do
c            call plotit_axes(tdata,tdatb,-180.,180.)
c            call plotit(tdat,bb,pdat,nf,'Frequency Domain',
c     x   'frequency (Hz)',xlabel(l),3,0,0.05,0,59+l*3)
            do n=1,nf
              bb(n)=ar(n)
            end do
            call fft2(ar,ai,npad)
c  normalization factor:
c  divide by npad for fft2 routine normalization
c  mult by nnf/nf to compensate for losing high freqs
c  mult by 2 if cosine**2 taper is applied to spectral RF
            fac=2.*float(nnf)/(float(npad)*float(nf))
c            fac=(4./3.)*float(nnf)/(float(npad)*float(nf))


c     11/17/19. VLL damn it, hard-wired pre-event time....
c     changing 200 to 1000
c            do i=1,nmm
c              bb(1000+i)=ar(i)*fac
c            end do
c            do i=1,1000
c              bb(1001-i)=ar(npad+1-i)*fac
c            end do
c            ttt(1)=-1000*dt
c            ttt(2)=dt
c            ttt(3)=dt
c  end of change from 200 to 1000
c
         do i=1,npost
          bb(npre+i)=ar(i)*fac
c          time(npre+i)=(i-1)*dt
        end do
        do i=1,npre
          bb(npre+1-i)=ar(npad+1-i)*fac
c          time(npre+1-i)=-i*dt
        end do
        ttt(1)=-npre*dt
        ttt(2)=dt
        ttt(3)=dt
c
c        11/21/19 changes to fix output times



            kz=(kaz-l)/2+1
c           do i=1,nnm 
		 do i =1,ntot
              rft(i,kz,l)=baz+50.*bb(i)*l
              t3=-time(i)
              iunit=11+l
              write(iunit,1022) t3,baz,bb(i)
            end do
            write(iunit,101) '>'
          endif
        end do
 1022 format(f7.3,f6.1,f9.5)
  102 format(a)
      end do  
      close(12)
      close(13)
      kaz=kaz/2
c      print *,kaz,' traces'   
c      do iagain=1,2
c      call plotit_axes(0.,0.,0.,0.)
c      do n=1,kaz
c        call plotit(rft(150,n,2),time(150),dum,nnn,
c     x   'Composite Transverse RF Section',
c     x   'Back-Azimuth (degrees CW from N)',
c     x  'time(sec)',2,0,0.0,0,21*(n/kaz))
c      end do
            
c     call plotit_axes(0.,0.,0.,0.)
c      do n=1,kaz
c        call plotit(rft(150,n,1),time(150),dum,nnn,
c     x   'Composite Radial RF Section',
c     x   'Back-Azimuth (degrees CW from N)',
c     x  'time(sec)',2,0,0.0,0,22*(n/kaz))
c      end do
c      end do
  109 format('Back Azimuth Centered on ',f4.0,'degrees')
c  OK, LETS TRY THE SAME THING, BUT WITH VELOCITY MIGRATION
      sdelay=0.
      pdelay=0.
      do i=1,nl
        sdelay=sdelay+(dz(i)/vs(i))
        pdelay=pdelay+(dz(i)/vp(i))
        tau(i)=sdelay-pdelay
        itau(i)=tau(i)/dt
        print *,itau(i),tau(i)
      end do
      itau(nlp)=4100
      print *, 'organ-pipe mode count at 1 Hz in stack of layers: S & P'
      print *,'These are the S & P traveltimes from the basal interface'
      print 104,sdelay,pdelay
  104 format(2f10.1)
      open(12,file='oumr_baz.grid',form='formatted')
      open(13,file='oumt_baz.grid',form='formatted')
      kaz=0
      do ibaz=ibazmax,ibazmin,ibazinc
        baz=ibaz
        do l=1,2
          do k=1,nlp
            do n=1,nf
              arm(n,l,k)=0.
              aim(n,l,k)=0.
              dam(n,l,k)=0.
            end do    
          end do    
          jrec=0    
          do i=1,irec
            test=abs(baz-bazs(i))
            if(test.lt.abazinc.or.test.gt.360.-abazinc) then
c            if(test.lt.20..or.test.gt.340.) then
              jrec=jrec+1
c  the stack-model velocities are in SI units, s_p is in sec/km
              slw=(slow(i)/1000.)       
c  compute the slowness-dependent time delays at the interfaces of stack-model
              do k=1,nlp
                coss=sqrt(1.-vs(k)*vs(k)*slw*slw)
                cosp=sqrt(1.-vp(k)*vp(k)*slw*slw)
c  trapdoor to avoid incoming slownesses that are evanescent (e.g. head waves)
c  basically, one stacks with the last acceptable interval velocity
                if(vp(k)*slw.lt.1.) then
                  gam(k)=(vp(k)*coss-vs(k)*cosp)/(vp(k)-vs(k))
                else
                  gam(k)=gam(k-1)
                endif
              end do
              psd(1)=gam(1)*tau(1)
              xik(1)=0.
              do k=2,nlp
                psd(k)=psd(k-1)+gam(k)*(tau(k)-tau(k-1))
                xik(k)=psd(k-1)/gam(k)-tau(k-1)
              end do
c  spline-interpolate the RF in freq domain
              do n=1,nf
                sss(n,1)=real(rfs(n,i,l))
                sss(n,2)=imag(rfs(n,i,l))
                sss(n,3)=drfs(n,i,l)
              end do
              call splneq(nf,sss(1,1),sss(1,4))
              call splneq(nf,sss(1,2),sss(1,5))
              call splneq(nf,sss(1,3),sss(1,6))
c  sum into nlp RF stacks, each for a different layer of the stacking model
              do k=1,nlp
                dfm=ddf/gam(k)
                dxm=1.0/gam(k)
                dcs=cos(2.*pi*dfm*xik(k)) 
                dsn=sin(2.*pi*dfm*xik(k)) 
c  we start at f=0., so that the migration phase factor = 1
c  csm and snm are real/imag parts of the migration phase factor
c  and the position in the spline is x=1  (fm is "migration frequency")
                csm=1.
                snm=0.
                xm=1.
                do n=1,nf
                  sig=evaleq(xm,nf,sss(1,3),sss(1,6),0,1.)
                  dre=evaleq(xm,nf,sss(1,1),sss(1,4),0,1.)
                  dim=evaleq(xm,nf,sss(1,2),sss(1,5),0,1.)
                  dre=dre/(sig*gam(k))
                  dim=dim/(sig*gam(k))
                  arm(n,l,k)=arm(n,l,k)+dre*csm+dim*snm
                  aim(n,l,k)=aim(n,l,k)-dre*snm+dim*csm
                  dam(n,l,k)=dam(n,l,k)+1.0/sig
                  csm=csm*dcs-snm*dsn
                  snm=snm*dcs+csm*dsn
                  xm=xm+dxm
                end do
              end do  
            endif
          end do
          if(jrec.ge.ibskip) then
            kaz=kaz+1
	    jjrec((kaz+1)/2)=jrec
	    jbaz((kaz+1)/2)=ibaz
c  do the expectation of variance for the weighted sum yourself
c  if you doubt the formula here -- basically the weighted terms in rf-sum
c  have unit variance, so that variance of the total variance is jrec/drf**2
c    start of k=1,nlp loop ccccccccc
            do k=1,nlp
cccccccccccccccccccccccccccccccccccc
              do n=1,nf
c	        if(n.ge.nf/2) then
c                  fac=cos(pi*(n-nf/2)/nf)**2
c	        else
c                  fac=1.0
c	        endif
                fac=cos(pi*(n-1)/(2*nf))**2
                arm(n,l,k)=arm(n,l,k)/dam(n,l,k)
                aim(n,l,k)=aim(n,l,k)/dam(n,l,k)
                ar(n)=fac*arm(n,l,k)
                ai(n)=-fac*aim(n,l,k)
                ar(npad+2-n)=fac*arm(n,l,k)
                ai(npad+2-n)=fac*aim(n,l,k)
                dam(n,l,k)=fac*sqrt(float(jrec))/dam(n,l,k)
              end do
              do n=nf+1,npad-nf+1
                ar(n)=0.d0
                ai(n)=0.d0
              end do
              write(string,109) baz
              bbmx=0.
              do n=1,nf
                bb(n)=dsqrt(ar(n)**2+ai(n)**2)
                bbmx=amax1(bbmx,bb(n))
                pdat(n)=dam(n,l,k)/sq2
              end do
              do n=1,npad
                if(bb(n).gt.1e-4) then
                  pdat(n)=con*(pdat(n)/bb(n)) 
                else
                  pdat(n)=360.
                endif
                bb(n)=con*datan2(ai(n),ar(n))
              end do
              do n=1,nf
                bb(n)=ar(n)
              end do
c              call plotit_axes(0.,0.,0.,0.)
              call fft2(ar,ai,npad)
c  normalization factor:
c  divide by npad for fft2 routine normalization
c  mult by nnf/nf to compensate for losing high freqs
c  mult by 2 if cosine**2 taper is applied to spectral RF
              fac=2.*float(nnf)/(float(npad)*float(nf))
c              fac=(4./3.)*float(nnf)/(float(npad)*float(nf))


c     11/17/19. VLL damn it, hard-wired pre-event time....
c     changing 200 to 1000
c
c              do i=1,nmm
c                bb(1000+i)=ar(i)*fac
c              end do
c              do i=1,1000
c                bb(1001-i)=ar(npad+1-i)*fac
c              end do
            do i=1,npost
              bb(npre+i)=ar(i)*fac
            end do
            do i=1,npre
              bb(npre+1-i)=ar(npad+1-i)*fac
            end do

c              ttt(1)=-1000*dt
		  ttt(1)=-npre*dt
              ttt(2)=dt
              ttt(3)=dt
c
c            end of change 

              kz=(kaz-l)/2+1
c  we store nlp migrated receiver functions in rft(i,kz+k,l), k=1,nlp
c              do i=1,nnm 
		   do i=1,ntot
                rft(i,kz+k,l)=bb(i)
              end do
 1004 format(9g13.3)
c    end of k=1,nlp loop cccccccc
            end do   
ccccccccccccccccccccccccccccccccc  
c            do i=1,nnm
             do i=1,ntot
              rft(i,kz,l)=0.
            end do         
            i1=1
c   11/17/19 another place where something was hard-wired to 200, changed to 1000
c   11/21/19 changed to flexible values "ntot" and "npre"
            do k=1,nlp
              i2=min(ntot,npre+itau(k))
              print *,i1,i2
              do i=i1,i2 
                rft(i,kz,l)=rft(i,kz,l)+rft(i,kz+k,l)
              end do
              i1=npre+1+itau(k)
            end do
          endif
c
c   end of change
c
          print *,' l = ',l
        end do   ! l-loop over radial & transverse
c  Now thru the k=1,nlp loop, we assemble the segments of the migrated RF
c  the noncausal potion of the RF is copied from the shallow-layer migration 
        if(jrec.ge.ibskip) then
c          call plotit_axes(0.,0.,0.,0.)
          xx(1)=tau(nlp)
          xx(2)=tau(nlp)
          yy(1)=-0.1
          yy(2)=0.1
c          do i=1,nnm
           do i=1,ntot
            bb(i)=rft(i,kz,1)
            rft(i,kz,1)=baz+rft(i,kz,1)*50.
            t3=-time(i)
            write(12,1022) t3,baz,bb(i)
          end do
c  note that we boost the transverse RF by factor of 2 for plotting in 
c the interactive plotting routine plotit, but TRF
c  is written to disk without amplification
c          do i=1,nnm
           do i=1,ntot
            bb(i)=rft(i,kz,2)
            rft(i,kz,2)=baz+rft(i,kz,2)*100.
            t3=-time(i)
            write(13,1022) t3,baz,bb(i)
          end do
          write(12,101) '>'
          write(13,101) '>'
        endif
      end do  
      close(12)
      close(13)

      kaz=kaz/2
      print *,kaz,' traces'   
c      do iagain=1,2
c        call plotit_axes(0.,0.,0.,0.)
c        do n=1,kaz
c          call plotit(rft(150,n,2),time(150),dum,nnn,
c     x   'Composite Transverse RF Section',
c     x   'Back-Azimuth (degrees CW from N)',
c     x  'time(sec)',2,0,0.0,0,21*(n/kaz))
c        end do
            
c        call plotit_axes(0.,0.,0.,0.)
c        do n=1,kaz
c          call plotit(rft(150,n,1),time(150),dum,nnn,
c     x   'Composite Radial RF Section',
c     x   'Back-Azimuth (degrees CW from N)',
c     x  'time(sec)',2,0,0.0,0,22*(n/kaz))
c        end do
c      end do
c  OK, LETS TRY THE SAME THING FOR THE EPICENTRAL SWEEP, 
c  BUT WITH VELOCITY MIGRATION  -- THIS FORMS A CHECK ON THE ALGORITHM
      open(12,file='oumr_epi.grid',form='formatted')
      open(13,file='oumt_epi.grid',form='formatted')
cc  NOW PLOT TRACES THAT ARE BINNED WITH EPICENTRAL DISTANCE & MiGRATED
      print *,'binned events versus baz'
      write(6,1008) (jjrec(i),i=1,kaz)
      write(6,1008) (jbaz(i),i=1,kaz)
 1008 format(20i4)
      print *,'compute RFs binned with epicentral distance'
      print *,'enter back-azimuth limits ib1,ib2 (integers!)'
      print *,' ib1=ib2 -> 0,360 and 360-wraparound if ib1 > ib2 '
      read(5,*) ib1,ib2
      if(ib1.eq.ib2) then
        ick=1
        baz1=0.
        baz2=360.
      elseif(ib1.lt.ib2) then
        ick=1
        baz1=ib1
        baz2=ib2
      else
        ick=2
        baz1=ib1
        baz2=ib2
      endif
      kaz=0
c  code for changing the epicentral bin width
      iep1=0
      iep2=175
      idep=5
      print *,'epicentral bins from',iep1,' to',iep2
      print *,'with halfbin width',idep
      print *,' do you want to change this spacing? (1=yes)'
      read(5,*) ickk
      if(ickk.eq.1) then
 2020   print *,'enter delta range and spacing (iep1,iep2,idep)'
        read(5,*) iep1,iep2,idep
        nep=(iep2-iep1)/idep+1
        if(nep.gt.500.or.nep.lt.1) then
          print *,'HEY! There are ',nep,' traces!'
          go to 2020
	endif
      endif
c  end code for changing the epicentral bin width
      do iepi=iep1,iep2,idep
        epi=iepi
        do l=1,2
          do n=1,nf
            rf(n,l)=zero
            drf(n,l)=0.
          end do    
          do k=1,nlp
            do n=1,nf
              arm(n,l,k)=0.
              aim(n,l,k)=0.
              dam(n,l,k)=0.
            end do    
          end do    
          jrec=0    
          do i=1,irec
            test=abs(epi-gcarc(i))
            if(test.lt.float(idep)) then
c  branch to limit the back-azimuth range for this sum
              ick1=0
              if(ick.eq.1) then
                if(bazs(i).ge.baz1.and.bazs(i).le.baz2) ick1=1
              else
                if(bazs(i).ge.baz1.or.bazs(i).le.baz2) ick1=1
              endif
              if(ick1.eq.1) then               
                jrec=jrec+1
c  the stac  k-model velocities are in SI units, s_p is in sec/km
                slw=(slow(i)/1000.)       
c  compute   the slowness-dependent time delays at the interfaces of stack-model
                do k=1,nlp
                  coss=sqrt(1.-vs(k)*vs(k)*slw*slw)
                  cosp=sqrt(1.-vp(k)*vp(k)*slw*slw)
c  trapdoor to avoid incoming slownesses that are evanescent (e.g. head waves)
c  basically, one stacks with the last acceptable interval velocity
                  if(vp(k)*slw.lt.1.) then
                    gam(k)=(vp(k)*coss-vs(k)*cosp)/(vp(k)-vs(k))
                  else
                    gam(k)=gam(k-1)
                  endif
                end do
                psd(1)=gam(1)*tau(1)
                xik(1)=0.
                do k=2,nlp
                  psd(k)=psd(k-1)+gam(k)*(tau(k)-tau(k-1))
                  xik(k)=psd(k-1)/gam(k)-tau(k-1)
                end do
c  spline-interpolate the RF in freq domain
                do n=1,nf
                  sss(n,1)=real(rfs(n,i,l))
                  sss(n,2)=imag(rfs(n,i,l))
                  sss(n,3)=drfs(n,i,l)
                end do
                call splneq(nf,sss(1,1),sss(1,4))
                call splneq(nf,sss(1,2),sss(1,5))
                call splneq(nf,sss(1,3),sss(1,6))
c  sum into nlp RF stacks, each for a different layer of the stacking model
                do k=1,nlp
                  dfm=ddf/gam(k)
                  dxm=1.0/gam(k)
                  dcs=cos(2.*pi*dfm*xik(k)) 
                  dsn=sin(2.*pi*dfm*xik(k)) 
c  we start at f=0., so that the migration phase factor = 1
c  csm and snm are real/imag parts of the migration phase factor
c  and the position in the spline is x=1  (fm is "migration frequency")
                  csm=1.
                  snm=0.
                  xm=1.
                  do n=1,nf
                    sig=evaleq(xm,nf,sss(1,3),sss(1,6),0,1.)
                    dre=evaleq(xm,nf,sss(1,1),sss(1,4),0,1.)
                    dim=evaleq(xm,nf,sss(1,2),sss(1,5),0,1.)
                    dre=dre/(sig*gam(k))
                    dim=dim/(sig*gam(k))
                    arm(n,l,k)=arm(n,l,k)+dre*csm+dim*snm
                    aim(n,l,k)=aim(n,l,k)-dre*snm+dim*csm
                    dam(n,l,k)=dam(n,l,k)+1.0/sig
                    csm=csm*dcs-snm*dsn
                    snm=snm*dcs+csm*dsn
                    xm=xm+dxm
                  end do
                end do  
              endif
            endif
          end do
          if(jrec.ge.ibskip) then
            kaz=kaz+1
c  do the expectation of variance for the weighted sum yourself
c  if you doubt the formula here -- basically the weighted terms in rf-sum
c  have unit variance, so that variance of the total variance is jrec/drf**2
c    start of k=1,nlp loop ccccccccc
            do k=1,nlp
cccccccccccccccccccccccccccccccccccc
            do n=1,nf
c	      if(n.ge.nf/2) then
c                fac=cos(pi*(n-nf/2)/nf)**2
c	      else
c                fac=1.0
c	      endif
              fac=cos(pi*(n-1)/(2*nf))**2
              arm(n,l,k)=arm(n,l,k)/dam(n,l,k)
              aim(n,l,k)=aim(n,l,k)/dam(n,l,k)
              ar(n)=fac*arm(n,l,k)
              ai(n)=-fac*aim(n,l,k)
              ar(npad+2-n)=fac*arm(n,l,k)
              ai(npad+2-n)=fac*aim(n,l,k)
              dam(n,l,k)=fac*sqrt(float(jrec))/dam(n,l,k)
            end do
            do n=nf+1,npad-nf+1
              ar(n)=0.d0
              ai(n)=0.d0
            end do
            write(string,108) epi
            bbmx=0.
            do n=1,nf
              bb(n)=dsqrt(ar(n)**2+ai(n)**2)
              bbmx=amax1(bbmx,bb(n))
              pdat(n)=dam(n,l,k)/sq2
            end do
            do n=1,npad
              if(bb(n).gt.1e-4) then
                pdat(n)=con*(pdat(n)/bb(n)) 
              else
                pdat(n)=360.
              endif
              bb(n)=con*datan2(ai(n),ar(n))
            end do
            do n=1,nf
              bb(n)=ar(n)
            end do
c            call plotit_axes(0.,0.,0.,0.)
            call fft2(ar,ai,npad)
c  normalization factor:
c  divide by npad for fft2 routine normalization
c  mult by nnf/nf to compensate for losing high freqs
c  mult by 2 if cosine**2 taper is applied to spectral RF
            fac=2.*float(nnf)/(float(npad)*float(nf))
c            fac=(4./3.)*float(nnf)/(float(npad)*float(nf))


c     11/17/19. VLL damn it, hard-wired pre-event time....
c     changing 200 to 1000
c            do i=1,nmm
c              bb(1000+i)=ar(i)*fac
c            end do
c            do i=1,1000
c              bb(1001-i)=ar(npad+1-i)*fac
c            end do
c            ttt(1)=-1000*dt
            do i=1,npost
              bb(npre+i)=ar(i)*fac
            end do
            do i=1,npre
              bb(npre+1-i)=ar(npad+1-i)*fac
            end do
            ttt(1)=-npre*dt
            ttt(2)=dt
            ttt(3)=dt
c
c           end of change
c
            kz=(kaz-l)/2+1
c  we store nlp migrated receiver functions in rft(i,kz+k,l), k=1,nlp
c           do i=1,nnm 
            do i=1,ntot
              rft(i,kz+k,l)=bb(i)
            end do
c    end of k=1,nlp loop cccccccc
            end do   
ccccccccccccccccccccccccccccccccc  
c            do i=1,nnm
		do i=1,ntot
              rft(i,kz,l)=0.
            end do 
 
c     11/17/19. VLL damn it, hard-wired pre-event time....
c     changing 200 to 1000
c     changing to "npre"
       
            i1=1
            do k=1,nlp
              i2=min(ntot,npre+itau(k))
              print *,i1,i2
              do i=i1,i2 
                rft(i,kz,l)=rft(i,kz,l)+rft(i,kz+k,l)
              end do
              i1=npre+1+itau(k)
            end do
c   end change
c
          endif
          print *,' l = ',l
        end do   ! l-loop over radial & transverse
c  Now thru the k=1,nlp loop, we assemble the segments of the migrated RF
c  the noncausal potion of the RF is copied from the shallow-layer migration 
        if(jrec.ge.ibskip) then
c          call plotit_axes(0.,0.,0.,0.)
          xx(1)=tau(nlp)
          xx(2)=tau(nlp)
          yy(1)=-0.1
          yy(2)=0.1
c          do i=1,nnm
		do i=1,ntot
            bb(i)=rft(i,kz,1)
            rft(i,kz,1)=epi+rft(i,kz,1)*50.
            t3=-time(i)
            write(12,1022) t3,epi,bb(i)
          end do
c  note that we boost the transverse RF by factor of 2 for PLOTIT display
c          do i=1,nnm
		do i=1,ntot
            bb(i)=rft(i,kz,2)
            rft(i,kz,2)=epi+rft(i,kz,2)*100.
            t3=-time(i)
            write(13,1022) t3,epi,bb(i)
          end do
          write(12,101) '>'
          write(13,101) '>'
        endif
      end do  
      close(12)
      close(13)
      kaz=kaz/2
      print *,kaz,' traces'   
c      call plotit_axes(0.,0.,0.,0.)
  108 format('Epicenters Centered on ',f4.0,'degrees')
c      do n=1,kaz
c        call plotit(rft(150,n,2),time(150),dum,nnn,
c     x   'Composite Transverse RF Section',
c     x   'Epicentral Distance (degrees)',
c     x  'time(sec)',2,0,0.0,0,21*(n/kaz))
c      end do
            
c      call plotit_axes(0.,0.,0.,0.)
c      do n=1,kaz
c        call plotit(rft(150,n,1),time(150),dum,nnn,
c     x   'Composite Radial RF Section',
c     x   'Epicentral Distance (degrees)',
c     x  'time(sec)',2,0,0.0,0,22*(n/kaz))
c      end do

      stop
      end
c
      subroutine taper(n,nwin,el)
c
c  to generate slepian tapers
c  ta is a real*4 array
c
c         j. park
c
      real*8 el,a,z,pi,ww,cs,ai,an,eps,rlu,rlb
      real*8 dfac,drat,gamma,bh,ell
      common/nnaamm/iwflag
      common/npiprol/anpi
      common/tapsum/tapsum(20),ell(20)
      common/work/ip(16400)        
      common/taperzz/z(262144)  ! we use this common block for scratch space
      common/stap2/ta(16400,8)
      dimension a(16400,8),el(10)
      data iwflag/0/,pi/3.14159265358979d0/
      equivalence (a(1,1),ta(1,1))
      an=dfloat(n)
      ww=dble(anpi)/an
      cs=dcos(2.d0*pi*ww)
c initialize matrix for eispack subroutine
c      print *,'ww,cs,an',ww,cs,an
      do i=0,n-1
        ai=dfloat(i)
        a(i+1,1)=-cs*((an-1.d0)/2.d0-ai)**2
        a(i+1,2)=-ai*(an-ai)/2.d0
        a(i+1,3)=a(i+1,2)**2        ! required by eispack routine
      end do
      eps=1.e-13
      m11=1
      call tridib(n,eps,a(1,1),a(1,2),a(1,3),rlb,rlu,m11,nwin,el,ip,
     x       ierr,a(1,4),a(1,5))
c      print *,ierr,rlb,rlu
      print *,'eigenvalues for the eigentapers'
c      print *,(el(i),i=1,nwin)
      call tinvit(n,n,a(1,1),a(1,2),a(1,3),nwin,el,ip,z,ierr,
     x            a(1,4),a(1,5),a(1,6),a(1,7),a(1,8))      
c  we calculate the eigenvalues of the dirichlet-kernel problem
c  i.e. the bandwidth retention factors
c  from slepian 1978 asymptotic formula, gotten from thomson 1982 eq 2.5
c  supplemented by the asymptotic formula for k near 2n from slepian 1978 eq 61
      dfac=an*pi*ww
      drat=8.d0*dfac
      dfac=4.d0*dsqrt(pi*dfac)*dexp(-2.d0*dfac)
      do k=1,nwin
        el(k)=1.d0-dfac
        dfac=dfac*drat/k  ! is this correct formula? yes,but fails as k -> 2n
      end do
c      print *,'eigenvalues for the eigentapers (small k)'
c      print *,(el(i),i=1,nwin)
      gamma=dlog(8.d0*an*dsin(2.d0*pi*ww))+0.5772156649d0
      do k=1,nwin
        bh=-2.d0*pi*(an*ww-dfloat(k-1)/2.d0-.25d0)/gamma
        ell(k)=1.d0/(1.d0+dexp(pi*bh))
      end do
      do i=1,nwin
        el(i)=dmax1(ell(i),el(i))
      end do     
c   normalize the eigentapers to preserve power for a white process
c   i.e. they have rms value unity
      do k=1,nwin
        kk=(k-1)*n
        tapsum(k)=0.
        tapsq=0.
        do i=1,n
          aa=z(kk+i)
          ta(i,k)=aa
          tapsum(k)=tapsum(k)+aa
          tapsq=tapsq+aa*aa
        end do
        aa=sqrt(tapsq/n)
        tapsum(k)=tapsum(k)/aa
        do i=1,n
          ta(i,k)=ta(i,k)/aa
        end do
      end do
c      print *,'tapsum',(tapsum(i),i=1,nwin)
  101 format(80a)
c  refft preserves amplitudes with zeropadding
c  zum beispiel: a(i)=1.,i=1,100 will transform at zero freq to b(f=0.)=100
c  no matter how much zero padding is done
c  therefore we need not doctor the taper normalization,
c  but wait until the fft to force the desired units
      iwflag=1
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c               routine
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fft2(ar,ai,n)
c  fft routine with 2 real input arrays rather than complex - see ffttwo
c  for comments. 
c  fft2 subroutine mults by exp(i\omega t)
c  OUTPUT: f=0 in ar(1), f=f_N in ar(n/2+1)
c          ar(i)=ar(n+2-i), ai(i)=-ai(n+2-i)
c  fft2 is NOT a unitary transform, mults the series by sqrt(n)
c  the inverse FT can be effected by running fft2 on the conjugate of
c  the FFT-expansion, then taking the the conjugate of the output, 
c  and dividing thru by N. to wit:
c
c   assume Xr, Xi is in freq domain, xr, xi in the time domain
c
c   (Xr,Xi)=fft2(xr,xi,N)
c   (xr,-xi)=fft2(Xr,-Xi,N)/N
c
      implicit real*8 (a-h,o-z)
      implicit integer*4 (i-n)
      dimension ar(1),ai(1)
      mex=dlog(dble(float(n)))/.693147d0
      nv2=n/2
      nm1=n-1
      j=1
      do 7 i=1,nm1
      if(i .ge. j) go to 5
      tr=ar(j)
      ti=ai(j)
      ar(j)=ar(i)
      ai(j)=ai(i)
      ar(i)=tr
      ai(i)=ti
   5  k=nv2
   6  if(k .ge. j) go to 7
      j=j-k
      k=k/2
      go to 6
   7  j=j+k
      pi=3.14159265358979d0
      do 20 l=1,mex
      le=2**l
      le1=le/2
      ur=1.0
      ui=0.
      wr=dcos(pi/le1 )
      wi=dsin (pi/le1)
      do 20 j=1,le1
      do 10 i=j,n,le
      ip=i+le1
      tr=ar(ip)*ur - ai(ip)*ui
      ti=ai(ip)*ur + ar(ip)*ui
      ar(ip)=ar(i)-tr
      ai(ip)=ai(i) - ti
      ar(i)=ar(i)+tr
      ai(i)=ai(i)+ti
  10  continue
      utemp=ur
      ur=ur*wr - ui*wi
      ui=ui*wr + utemp*wi
  20  continue
      return
      end
      subroutine read_ttimes
      implicit real*4(a-h,o-z)
      implicit integer*4(i-n)
      common/tt_p/d_p(200),d_pkp(200),s_p(200),s_pkp(200)
      common/tt_pp/d_pp(200),d_pkkp(200),s_pp(200),s_pkkp(200)
c  read in the P- and PKP-slowness values from taup_time 
c  this output has ray slowness in seconds/degree, so we must convert to
c  phase velocity in km/sec  -- pvel=(6371*\pi/180)/slowness
c  the code interpolates across the upper-mantle triplication for DELTA=19-27
c  The user chooses P or PKP or PP or PKKP based on the T1 parameter in the SAC header
c  The values of Taup_time specify epicentral depth 50km, 
c  so these values should be used for shallow quakes
c  PKKP slownesses NOT implemented at this time
c
c  no need to initialize, all taup files are 180 lines
      open(7,file='./taup_P.dat',
     x  form='formatted')
      do i=1,180
        read(7,*,end=110) d_p(i),s_p(i)
      end do
  110 close(7)
      print *,i,' = counter at end of ttimes file'
  120 n_p=180
c  fill in the head-wave slowness for DELTA.lt.14
c  13.2 sec/degree is 8.423+eps km/sec
      do i=1,13
        s_p(i)=13.2
      end do
c  interpolates the triplication
c  linearly interpolate slowness between Delta=17 and p=11.06 sec/degree
c                                    and Delta=28 and p=8.92 sec/degree
      do i=18,27
        s_p(i)=11.06-(i-17)*(11.06-8.92)/11.
      end do
      open(7,file='./taup_PKP.dat',
     x  form='formatted')
c  this file starts with DELTA=1  - trending from PKiKP to PKIKP to PKP and back to PKIKP
      do i=1,180
        read(7,*,end=130) d_pkp(i),s_pkp(i)
      end do
  130 close(7)
      print *,i,' = counter at end of ttimes file'
  140 n_pkp=180
      open(7,file='./taup_PP.dat',
     x  form='formatted')
      do i=1,180
        read(7,*) d_pp(i),s_pp(i)
      end do
      close(7)
c  fill in the head-wave slowness for DELTA.lt.28
c  13.2 sec/degree is 8.423+eps km/sec will be faster that any Moho headwave
      do i=1,27
        s_pp(i)=13.2
      end do
c  interpolates the triplication
c  linearly interpolate slowness between Delta=34 and p=11.06 sec/degree
c                                    and Delta=56 and p=8.92 sec/degree
      do i=35,55
        s_pp(i)=11.06-(i-34)*(11.06-8.92)/22.
      end do
      n_pp=180
      
      factor=180./(3.14159265358979*6371.)
      do i=1,180
        s_p(i)=s_p(i)*factor
        s_pkp(i)=s_pkp(i)*factor
        s_pp(i)=s_pp(i)*factor
        s_pkkp(i)=s_pkkp(i)*factor
      end do
      print *,'n_p,n_pkp',n_p,n_pkp
c      call plotit_axes(0.,0.,0.,0.)
c      call plotit(d_p,s_p,dum,n_p,' ',' ',' ',2,0,0,0,0)
c      call plotit(d_pkp,s_pkp,dum,n_pkp,' ',
c     x 'Epicentral Distance (degrees)','Slowness (sec/km)',2,0,0,0,21)
c      call plotit(d_pp,s_pp,dum,n_pp,' ',' ',' ',2,0,0,0,0)
c      call plotit(d_pkkp,s_pkkp,dum,n_pkkp,' ',
c     x 'Epicentral Distance (degrees)','Slowness (sec/km)',2,0,0,0,22)
      return
      end
c
c
c
c
c  SUBROUTINES FROM EISPACK AND OTHER NECESSARY COMPONENTS
c
c
c
      SUBROUTINE TINVIT(NM,N,D,E,E2,M,W,IND,Z,
     X                  IERR,RV1,RV2,RV3,RV4,RV6)
C
      INTEGER I,J,M,N,P,Q,R,S,II,IP,JJ,NM,ITS,TAG,IERR,GROUP
      REAL*8 D(N),E(N),E2(N),W(M),Z(NM,M),
     X       RV1(N),RV2(N),RV3(N),RV4(N),RV6(N)
      REAL*8 U,V,UK,XU,X0,X1,EPS2,EPS3,EPS4,NORM,ORDER,MACHEP
      REAL*8 DSQRT,DABS,DFLOAT
      INTEGER IND(M)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE INVERSE ITERATION TECH-
C     NIQUE IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVECTORS OF A TRIDIAGONAL
C     SYMMETRIC MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES,
C     USING INVERSE ITERATION.
C
C     ON INPUT:
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT;
C
C        N IS THE ORDER OF THE MATRIX;
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY;
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E,
C          WITH ZEROS CORRESPONDING TO NEGLIGIBLE ELEMENTS OF E.
C          E(I) IS CONSIDERED NEGLIGIBLE IF IT IS NOT LARGER THAN
C          THE PRODUCT OF THE RELATIVE MACHINE PRECISION AND THE SUM
C          OF THE MAGNITUDES OF D(I) AND D(I-1).  E2(1) MUST CONTAIN
C          0.0D0 IF THE EIGENVALUES ARE IN ASCENDING ORDER, OR 2.0D0
C          IF THE EIGENVALUES ARE IN DESCENDING ORDER.  IF  BISECT,
C          TRIDIB, OR  IMTQLV  HAS BEEN USED TO FIND THE EIGENVALUES,
C          THEIR OUTPUT E2 ARRAY IS EXACTLY WHAT IS EXPECTED HERE;
C
C        M IS THE NUMBER OF SPECIFIED EIGENVALUES;
C
C        W CONTAINS THE M EIGENVALUES IN ASCENDING OR DESCENDING ORDER;
C
C        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
C          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
C          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
C          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.
C
C     ON OUTPUT:
C
C        ALL INPUT ARRAYS ARE UNALTERED;
C
C        Z CONTAINS THE ASSOCIATED SET OF ORTHONORMAL EIGENVECTORS.
C          ANY VECTOR WHICH FAILS TO CONVERGE IS SET TO ZERO;
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          -R         IF THE EIGENVECTOR CORRESPONDING TO THE R-TH
C                     EIGENVALUE FAILS TO CONVERGE IN 5 ITERATIONS;
C
C        RV1, RV2, RV3, RV4, AND RV6 ARE TEMPORARY STORAGE ARRAYS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC
C                ON S360 ::::::::::
C  FOR F_FLOATING DEC FORTRAN
C      DATA MACHEP/1.1D-16/
C  FOR G_FLOATING DEC FORTRAN
       DATA MACHEP/1.25D-15/
C
      IERR = 0
      IF (M .EQ. 0) GO TO 1001
      TAG = 0
      ORDER = 1.0D0 - E2(1)
      Q = 0
C     :::::::::: ESTABLISH AND PROCESS NEXT SUBMATRIX ::::::::::
  100 P = Q + 1
C
      DO 120 Q = P, N
         IF (Q .EQ. N) GO TO 140
         IF (E2(Q+1) .EQ. 0.0D0) GO TO 140
  120 CONTINUE
C     :::::::::: FIND VECTORS BY INVERSE ITERATION ::::::::::
  140 TAG = TAG + 1
      S = 0
C
      DO 920 R = 1, M
         IF (IND(R) .NE. TAG) GO TO 920
         ITS = 1
         X1 = W(R)
         IF (S .NE. 0) GO TO 510
C     :::::::::: CHECK FOR ISOLATED ROOT ::::::::::
         XU = 1.0D0
         IF (P .NE. Q) GO TO 490
         RV6(P) = 1.0D0
         GO TO 870
  490    NORM = DABS(D(P))
         IP = P + 1
C
         DO 500 I = IP, Q
  500    NORM = NORM + DABS(D(I)) + DABS(E(I))
C     :::::::::: EPS2 IS THE CRITERION FOR GROUPING,
C                EPS3 REPLACES ZERO PIVOTS AND EQUAL
C                ROOTS ARE MODIFIED BY EPS3,
C                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW ::::::::::
         EPS2 = 1.0D-3 * NORM
         EPS3 = MACHEP * NORM
         UK = DBLE(Q-P+1)
         EPS4 = UK * EPS3
         UK = EPS4 / DSQRT(UK)
         S = P
  505    GROUP = 0
         GO TO 520
C     :::::::::: LOOK FOR CLOSE OR COINCIDENT ROOTS ::::::::::
  510    IF (DABS(X1-X0) .GE. EPS2) GO TO 505
         GROUP = GROUP + 1
         IF (ORDER * (X1 - X0) .LE. 0.0D0) X1 = X0 + ORDER * EPS3
C     :::::::::: ELIMINATION WITH INTERCHANGES AND
C                INITIALIZATION OF VECTOR ::::::::::
  520    V = 0.0D0
C
         DO 580 I = P, Q
            RV6(I) = UK
            IF (I .EQ. P) GO TO 560
            IF (DABS(E(I)) .LT. DABS(U)) GO TO 540
C     :::::::::: WARNING -- A DIVIDE CHECK MAY OCCUR HERE IF
C                E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY ::::::::::
            XU = U / E(I)
            RV4(I) = XU
            RV1(I-1) = E(I)
            RV2(I-1) = D(I) - X1
            RV3(I-1) = 0.0D0
            IF (I .NE. Q) RV3(I-1) = E(I+1)
            U = V - XU * RV2(I-1)
            V = -XU * RV3(I-1)
            GO TO 580
  540       XU = E(I) / U
            RV4(I) = XU
            RV1(I-1) = U
            RV2(I-1) = V
            RV3(I-1) = 0.0D0
  560       U = D(I) - X1 - XU * V
            IF (I .NE. Q) V = E(I+1)
  580    CONTINUE
C
         IF (U .EQ. 0.0D0) U = EPS3
         RV1(Q) = U
         RV2(Q) = 0.0D0
         RV3(Q) = 0.0D0
C     :::::::::: BACK SUBSTITUTION
C                FOR I=Q STEP -1 UNTIL P DO -- ::::::::::
  600    DO 620 II = P, Q
            I = P + Q - II
            RV6(I) = (RV6(I) - U * RV2(I) - V * RV3(I)) / RV1(I)
            V = U
            U = RV6(I)
  620    CONTINUE
C     :::::::::: ORTHOGONALIZE WITH RESPECT TO PREVIOUS
C                MEMBERS OF GROUP ::::::::::
         IF (GROUP .EQ. 0) GO TO 700
         J = R
C
         DO 680 JJ = 1, GROUP
  630       J = J - 1
            IF (IND(J) .NE. TAG) GO TO 630
            XU = 0.0D0
C
            DO 640 I = P, Q
  640       XU = XU + RV6(I) * Z(I,J)
C
            DO 660 I = P, Q
  660       RV6(I) = RV6(I) - XU * Z(I,J)
C
  680    CONTINUE
C
  700    NORM = 0.0D0
C
         DO 720 I = P, Q
  720    NORM = NORM + DABS(RV6(I))
C
         IF (NORM .GE. 1.0D0) GO TO 840
C     :::::::::: FORWARD SUBSTITUTION ::::::::::
         IF (ITS .EQ. 5) GO TO 830
         IF (NORM .NE. 0.0D0) GO TO 740
         RV6(S) = EPS4
         S = S + 1
         IF (S .GT. Q) S = P
         GO TO 780
  740    XU = EPS4 / NORM
C
         DO 760 I = P, Q
  760    RV6(I) = RV6(I) * XU
C     :::::::::: ELIMINATION OPERATIONS ON NEXT VECTOR
C                ITERATE ::::::::::
  780    DO 820 I = IP, Q
            U = RV6(I)
C     :::::::::: IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE
C                WAS PERFORMED EARLIER IN THE
C                TRIANGULARIZATION PROCESS ::::::::::
            IF (RV1(I-1) .NE. E(I)) GO TO 800
            U = RV6(I-1)
            RV6(I-1) = RV6(I)
  800       RV6(I) = U - RV4(I) * RV6(I-1)
  820    CONTINUE
C
         ITS = ITS + 1
         GO TO 600
C     :::::::::: SET ERROR -- NON-CONVERGED EIGENVECTOR ::::::::::
  830    IERR = -R
         XU = 0.0D0
         GO TO 870
C     :::::::::: NORMALIZE SO THAT SUM OF SQUARES IS
C                1 AND EXPAND TO FULL ORDER ::::::::::
  840    U = 0.0D0
C
         DO 860 I = P, Q
  860    U = U + RV6(I)**2
C
         XU = 1.0D0 / DSQRT(U)
C
  870    DO 880 I = 1, N
  880    Z(I,R) = 0.0D0
C
         DO 900 I = P, Q
  900    Z(I,R) = RV6(I) * XU
C
         X0 = X1
  920 CONTINUE
C
      IF (Q .LT. N) GO TO 100
 1001 RETURN
C     :::::::::: LAST CARD OF TINVIT ::::::::::
      END
      SUBROUTINE TRIDIB(N,EPS1,D,E,E2,LB,UB,M11,M,W,IND,IERR,RV4,RV5)
C
      INTEGER I,J,K,L,M,N,P,Q,R,S,II,M1,M2,M11,M22,TAG,IERR,ISTURM
      REAL*8 D(N),E(N),E2(N),W(M),RV4(N),RV5(N)
      REAL*8 U,V,LB,T1,T2,UB,XU,X0,X1,EPS1,MACHEP
      REAL*8 DABS,DMAX1,DMIN1,DFLOAT
      INTEGER IND(M)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BISECT,
C     NUM. MATH. 9, 386-393(1967) BY BARTH, MARTIN, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 249-256(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVALUES OF A TRIDIAGONAL
C     SYMMETRIC MATRIX BETWEEN SPECIFIED BOUNDARY INDICES,
C     USING BISECTION.
C
C     ON INPUT:
C
C        N IS THE ORDER OF THE MATRIX;
C
C        EPS1 IS AN ABSOLUTE ERROR TOLERANCE FOR THE COMPUTED
C          EIGENVALUES.  IF THE INPUT EPS1 IS NON-POSITIVE,
C          IT IS RESET FOR EACH SUBMATRIX TO A DEFAULT VALUE,
C          NAMELY, MINUS THE PRODUCT OF THE RELATIVE MACHINE
C          PRECISION AND THE 1-NORM OF THE SUBMATRIX;
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY;
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2(1) IS ARBITRARY;
C
C        M11 SPECIFIES THE LOWER BOUNDARY INDEX FOR THE DESIRED
C          EIGENVALUES;
C
C        M SPECIFIES THE NUMBER OF EIGENVALUES DESIRED.  THE UPPER
C          BOUNDARY INDEX M22 IS THEN OBTAINED AS M22=M11+M-1.
C
C     ON OUTPUT:
C
C        EPS1 IS UNALTERED UNLESS IT HAS BEEN RESET TO ITS
C          (LAST) DEFAULT VALUE;
C
C        D AND E ARE UNALTERED;
C
C        ELEMENTS OF E2, CORRESPONDING TO ELEMENTS OF E REGARDED
C          AS NEGLIGIBLE, HAVE BEEN REPLACED BY ZERO CAUSING THE
C          MATRIX TO SPLIT INTO A DIRECT SUM OF SUBMATRICES.
C          E2(1) IS ALSO SET TO ZERO;
C
C        LB AND UB DEFINE AN INTERVAL CONTAINING EXACTLY THE DESIRED
C          EIGENVALUES;
C
C        W CONTAINS, IN ITS FIRST M POSITIONS, THE EIGENVALUES
C          BETWEEN INDICES M11 AND M22 IN ASCENDING ORDER;
C
C        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
C          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
C          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
C          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.;
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          3*N+1      IF MULTIPLE EIGENVALUES AT INDEX M11 MAKE
C                     UNIQUE SELECTION IMPOSSIBLE,
C          3*N+2      IF MULTIPLE EIGENVALUES AT INDEX M22 MAKE
C                     UNIQUE SELECTION IMPOSSIBLE;
C
C        RV4 AND RV5 ARE TEMPORARY STORAGE ARRAYS.
C
C     NOTE THAT SUBROUTINE TQL1, IMTQL1, OR TQLRAT IS GENERALLY FASTER
C     THAN TRIDIB, IF MORE THAN N/4 EIGENVALUES ARE TO BE FOUND.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC
C                ON S360 ::::::::::
C  FOR F_FLOATING DEC FORTRAN
C      DATA MACHEP/1.1D-16/
C  FOR G_FLOATING DEC FORTRAN
       DATA MACHEP/1.25D-15/
C
      IERR = 0
      TAG = 0
      XU = D(1)
      X0 = D(1)
      U = 0.0D0
C     :::::::::: LOOK FOR SMALL SUB-DIAGONAL ENTRIES AND DETERMINE AN
C                INTERVAL CONTAINING ALL THE EIGENVALUES ::::::::::
      DO 40 I = 1, N
         X1 = U
         U = 0.0D0
         IF (I .NE. N) U = DABS(E(I+1))
         XU = DMIN1(D(I)-(X1+U),XU)
         X0 = DMAX1(D(I)+(X1+U),X0)
         IF (I .EQ. 1) GO TO 20
         IF (DABS(E(I)) .GT. MACHEP * (DABS(D(I)) + DABS(D(I-1))))
     X      GO TO 40
   20    E2(I) = 0.0D0
   40 CONTINUE
C
      X1 = DMAX1(DABS(XU),DABS(X0)) * MACHEP * DBLE(N)
      XU = XU - X1
      T1 = XU
      X0 = X0 + X1
      T2 = X0
C     :::::::::: DETERMINE AN INTERVAL CONTAINING EXACTLY
C                THE DESIRED EIGENVALUES ::::::::::
      P = 1
      Q = N
      M1 = M11 - 1
      IF (M1 .EQ. 0) GO TO 75
      ISTURM = 1
   50 V = X1
      X1 = XU + (X0 - XU) * 0.5D0
      IF (X1 .EQ. V) GO TO 980
      GO TO 320
   60 IF (S - M1) 65, 73, 70
   65 XU = X1
      GO TO 50
   70 X0 = X1
      GO TO 50
   73 XU = X1
      T1 = X1
   75 M22 = M1 + M
      IF (M22 .EQ. N) GO TO 90
      X0 = T2
      ISTURM = 2
      GO TO 50
   80 IF (S - M22) 65, 85, 70
   85 T2 = X1
   90 Q = 0
      R = 0
C     :::::::::: ESTABLISH AND PROCESS NEXT SUBMATRIX, REFINING
C                INTERVAL BY THE GERSCHGORIN BOUNDS ::::::::::
  100 IF (R .EQ. M) GO TO 1001
      TAG = TAG + 1
      P = Q + 1
      XU = D(P)
      X0 = D(P)
      U = 0.0D0
C
      DO 120 Q = P, N
         X1 = U
         U = 0.0D0
         V = 0.0D0
         IF (Q .EQ. N) GO TO 110
         U = DABS(E(Q+1))
         V = E2(Q+1)
  110    XU = DMIN1(D(Q)-(X1+U),XU)
         X0 = DMAX1(D(Q)+(X1+U),X0)
         IF (V .EQ. 0.0D0) GO TO 140
  120 CONTINUE
C
  140 X1 = DMAX1(DABS(XU),DABS(X0)) * MACHEP
      IF (EPS1 .LE. 0.0D0) EPS1 = -X1
      IF (P .NE. Q) GO TO 180
C     :::::::::: CHECK FOR ISOLATED ROOT WITHIN INTERVAL ::::::::::
      IF (T1 .GT. D(P) .OR. D(P) .GE. T2) GO TO 940
      M1 = P
      M2 = P
      RV5(P) = D(P)
      GO TO 900
  180 X1 = X1 * DBLE(Q-P+1)
      LB = DMAX1(T1,XU-X1)
      UB = DMIN1(T2,X0+X1)
      X1 = LB
      ISTURM = 3
      GO TO 320
  200 M1 = S + 1
      X1 = UB
      ISTURM = 4
      GO TO 320
  220 M2 = S
      IF (M1 .GT. M2) GO TO 940
C     :::::::::: FIND ROOTS BY BISECTION ::::::::::
      X0 = UB
      ISTURM = 5
C
      DO 240 I = M1, M2
         RV5(I) = UB
         RV4(I) = LB
  240 CONTINUE
C     :::::::::: LOOP FOR K-TH EIGENVALUE
C                FOR K=M2 STEP -1 UNTIL M1 DO --
C                (-DO- NOT USED TO LEGALIZE -COMPUTED GO TO-) ::::::::::
      K = M2
  250    XU = LB
C     :::::::::: FOR I=K STEP -1 UNTIL M1 DO -- ::::::::::
         DO 260 II = M1, K
            I = M1 + K - II
            IF (XU .GE. RV4(I)) GO TO 260
            XU = RV4(I)
            GO TO 280
  260    CONTINUE
C
  280    IF (X0 .GT. RV5(K)) X0 = RV5(K)
C     :::::::::: NEXT BISECTION STEP ::::::::::
  300    X1 = (XU + X0) * 0.5D0
         IF ((X0 - XU) .LE. (2.0D0 * MACHEP *
     X      (DABS(XU) + DABS(X0)) + DABS(EPS1))) GO TO 420
C     :::::::::: IN-LINE PROCEDURE FOR STURM SEQUENCE ::::::::::
  320    S = P - 1
         U = 1.0D0
C
         DO 340 I = P, Q
            IF (U .NE. 0.0D0) GO TO 325
            V = DABS(E(I)) / MACHEP
            IF (E2(I) .EQ. 0.0D0) V = 0.0D0
            GO TO 330
  325       V = E2(I) / U
  330       U = D(I) - X1 - V
            IF (U .LT. 0.0D0) S = S + 1
  340    CONTINUE
C
         GO TO (60,80,200,220,360), ISTURM
C     :::::::::: REFINE INTERVALS ::::::::::
  360    IF (S .GE. K) GO TO 400
         XU = X1
         IF (S .GE. M1) GO TO 380
         RV4(M1) = X1
         GO TO 300
  380    RV4(S+1) = X1
         IF (RV5(S) .GT. X1) RV5(S) = X1
         GO TO 300
  400    X0 = X1
         GO TO 300
C     :::::::::: K-TH EIGENVALUE FOUND ::::::::::
  420    RV5(K) = X1
      K = K - 1
      IF (K .GE. M1) GO TO 250
C     :::::::::: ORDER EIGENVALUES TAGGED WITH THEIR
C                SUBMATRIX ASSOCIATIONS ::::::::::
  900 S = R
      R = R + M2 - M1 + 1
      J = 1
      K = M1
C
      DO 920 L = 1, R
         IF (J .GT. S) GO TO 910
         IF (K .GT. M2) GO TO 940
         IF (RV5(K) .GE. W(L)) GO TO 915
C
         DO 905 II = J, S
            I = L + S - II
            W(I+1) = W(I)
            IND(I+1) = IND(I)
  905    CONTINUE
C
  910    W(L) = RV5(K)
         IND(L) = TAG
         K = K + 1
         GO TO 920
  915    J = J + 1
  920 CONTINUE
C
  940 IF (Q .LT. N) GO TO 100
      GO TO 1001
C     :::::::::: SET ERROR -- INTERVAL CANNOT BE FOUND CONTAINING
C                EXACTLY THE DESIRED EIGENVALUES ::::::::::
  980 IERR = 3 * N + ISTURM
 1001 LB = T1
      UB = T2
      RETURN
C     :::::::::: LAST CARD OF TRIDIB ::::::::::
      END
      subroutine sacread(name,a,ierr)
c  read a SAC-format file
      real*4 ahead
      character*120 name
      character*4 chead(158)
      common/header/ahead(158)
      dimension a(1)
      dimension iah(158)
      equivalence (ahead,iah),(ahead,chead)
      ierr=0
      open(8,file=name,access='direct',recl=512,err=999)
c  read the 158-word header, can use its info 
      read(8,rec=1,err=998)(ahead(i),i=1,128)
      print *,(ahead(i),i=50,55)
      print *,(iah(70+i),i=1,10)
      nscan=iah(80)
      irec=2
      read(8,rec=2,err=998)(ahead(i+128),i=1,30),(a(j),j=1,98)
c      print *,'rec=2'
      if(chead(145).ne.' GRN'.and.chead(145).ne.' grn') then
        nword=nscan+158
        nrec=(nword-1)/128+1
        last=nword-128*(nrec-1)
      else       
        print *,'reading greens functions' 
        nword=6*nscan+158
        nrec=(nword-1)/128+1
        last=nword-128*(nrec-1)
      endif
      if(nrec.gt.3) then
        do irec=3,nrec-1
c          print *,'rec=',irec
          read(8,rec=irec,err=998)(a((irec-3)*128+98+i),i=1,128)
        end do
      endif
c      print *,'rec=',nrec
      read(8,rec=nrec,err=998)(a((nrec-3)*128+98+i),i=1,last)
      close(8)
      return
  999 print *,'open error'
      ierr=1
      return
  998 print *,'read error: reading',irec,' out of',nrec
      ierr=1
      print *,irec
      close(8)
      return
      end

c----------------------extra subroutines that are a part of "Jlib" ---------------------

      subroutine splneq(nn, u, s)
cz$$$$ calls no other routines
c  finds coeffs for a spline interpolation of equally spaced data
c  based on 'spline interpolation ona  digital computer' by r.f.thompson
c  nn  number of data points (may be negative - see d1,d2 below)
c  u  array of function values to be interpolated, assumed to samples at equal i
c     intervals of a twice differentiable function.
c  s  array to hold computed values of 2nd derivative of spline fit at sample
c     points.  these values are required by evaleq  to interpolate
c  if the user wishes to force specific values on the derivatives at the end
c  points, he should put h*du(1)/dx  nad  h*du(n)/dx  in  s(1),s(2), then call
c  splneq  with nn=-number of terms in series. h = sample spacing in x.
c  normally the derivatives are found by fitting a parabola through the
c  1st and last 3 points.
c  if the number of terms is between 1 and 3, straight-line interpolation is don
      implicit integer*4 (i-n)
c     dimension u(1),s(1),a(13)
c   change to make gfortran stop complaining, VLL 04/20/16
       dimension u(3),s(3),a(13)
c
      n=iabs(nn)
      if (n.le.3) go to 5000
      d1=-0.5*u(3)  +2.0*u(2)  -1.5*u(1)
      dn= 0.5*u(n-2)-2.0*u(n-1)+1.5*u(n)
      if (nn.gt.0) go to 1000
      d1=s(1)
      dn=s(2) 
 1000 a(1)=2.0
      a(2)=3.5
      s(1)=u(2)-u(1)-d1
      s(2)=u(1)-2.0*u(2)+u(3)-0.5*s(1)
      n1=n-1
      do 3000 i=3,n1
      if (i.gt.13) go to 3000
      k=i
      a(k)=4.-1.0/a(k-1)
 3000 s(i)=u(i-1)-2.*u(i)+u(i+1)-s(i-1)/a(k-1)
      s(n)=u(n1)-u(n)+dn-s(n1)/a(k)
      s(n)=6.0*s(n)/(2.0-1.0/a(k))
      n2=n-2
c  compute 2nd derivatives by back-substitution
c  the array  a  tends to a constant (2+sqrt(3)) so only 13 elements are needed
      do 4000 j=1,n2
      i=n-j
      k=min0(i,k)
 4000 s(i)=(6.0*s(i)-s(i+1))/a(k)
      s(1)=3.0*s(1)-0.5*s(2)
      return
c  series too short for cubic spline.  fit straight lines.
 5000 do 5500 i=1,n
 5500 s(i)=0.0
      return
      end                                                               

      function evaleq(y, nn, u, s,ick,h)
cz$$$$ calls no other routines
c  performs spline interpolation of equally spaced data.
c  based on 'spline interpolation on a digital computer'( by r.f.thompson.
c  evaluates a spline interpolate in a set of equally spaced samples.
c  the routine  splneq  should be called first, to establish the array  s .
c  y  the  coordinate at which interpolate is required, with y=1 for 1st
c     sample point, y=2 for the 2nd, etc.  if actual spacing is  h  and  x1 is
c     the 1st sample coordinate use  y = 1.0 + (x-x1)/h
c  nn  number of samples of function in original set.
c  u  array containing function samples.
c  s  array of normalized 2nd derivatives, computed by  splneq.  the derivatives
c     have been multiplied by h**2, where h is the sample spacing.
c  if  y  is out of the range (1,nn), the 1st or last sample value is used.
      dimension u(1),s(1)
      data z3/.333333333333/,z6/.16666666666666667/
c
      if (y.le.1.0) go to 1500
      if (y.ge.float(nn)) go to 2000
      k1=y
      k=k1+1
      dk=k-y
      dk1=y-k1
      if(ick.eq.1) go to 1000
      ff1=s(k1)*dk*dk*dk
      ff2=s(k)*dk1*dk1*dk1
      evaleq=(dk*(6.0*u(k1)-s(k1))+ dk1*(u(k)*6.0-s(k)) + ff1 +ff2)/6.0
      return
c  evaluate the first derivative
 1000 a1=dk1*(1.-dk1/2.)-z3
      a2=dk1*dk1/2.-z6
      a3=u(k)-u(k1)+a1*s(k1)+s(k)*a2
      evaleq=a3/h
      return
c  out of range.  supply constant values
 1500 evaleq=u(1)
      return
 2000 evaleq=u(nn)
      return
      end          

c
c
      subroutine pause(string)
      character*(*) string
      print *,string
      read(5,*) 
      return
      end


