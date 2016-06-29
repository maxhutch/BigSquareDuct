c Trbulent DUCT flow, run_time statistics
c To monitor the mid plane flux of the streamwise velocity, the nelx,nely,nelz must be specified
c in the usrcheck. nelx,nely must be even
c A specific navier5.f is needed. This version contains the routines of comp_derivat and avg* 
c For different meshes the midplane must be adjusted in subroutines planar_s and planar_t manualy
C=======================================================================

      subroutine avg_stat_all

      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'   

      parameter (rline = lx1*lelx)
      parameter (sline = ly1*lely)
      common /plane_R/  wavg_pl_r(rline),rw1(rline),rw2(rline)
      common /plane_S/  wavg_pl_s(sline),sw1(sline),sw2(sline)
      integer mid_r, mid_s
      real dwdzmean
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

      integer nstat
      parameter (nstat = 71) ! Number of statistical fields to be saved

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

      integer nv
      parameter (nv = 10)   ! Interval of time record      
      ! Important: 
      ! remiander of time-steps/nv should be = 0 
      ! remainder of IOSTEP_avg/nv should be = 0

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

      integer n2ptc
      parameter (n2ptc = 250000000)   ! Interval to dump files for 2-point-correlations     

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C
      parameter (kx1=lx1,ky1=ly1,kz1=ly1,kx2=lx2,ky2=ly2,kz2=ly2)
      common /avgcmnr/ atime,timel
      common /c_p0/ p0(lx1,ly1,lz1,lelt)
      common /cvflow_r/ flow_rate,base_flow,domain_length,xsec
     $                , scale_vf(3)

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C
C ---------------------------------------------------------------------- C
C -------------- Define the various quantities on a 3D array------------ C
C ---------------------------------------------------------------------- C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

      real stat(kx1*ky1*kz1*lelt,nstat)

      real pm1(kx1*ky1*kz1*lelt)     
      real wk1(kx1*ky1*kz1)
      real wk2(kx1*ky1*kz1)

      real duidxj(lx1*ly1*lz1,lelt,1:3*ldim)
      real ur(lx1*ly1*lz1),us(lx1*ly1*lz1),ut(lx1*ly1*lz1)
      real vr(lx1*ly1*lz1),vs(lx1*ly1*lz1),vt(lx1*ly1*lz1)
      real wr(lx1*ly1*lz1),ws(lx1*ly1*lz1),wt(lx1*ly1*lz1)

      real u_avg, uu_avg
      real xlmin, xlmax, domain_x
      real ylmin, ylmax, domain_y
      real zlmin, zlmax, domain_z

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C
C ---------------------------------------------------------------------- C
C -------------- Define the various quantities on a 2D array ----------- C
C ---------------------------------------------------------------------- C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

      real stat_xy(ky1*lely*kx1*lelx,nstat)
      real xx(kx1*lelx),yy(ky1*lely),zz(kz1*lelz)
      real uzz(kx1*ky1*kz1*lelt)

      real w1(kx1,ky1,lelx,lely),w2(kx1,ky1,lelx,lely)
      real w1x(lx1*lelx), w2x(lx1*lelx)
      real w1y(ly1*lely), w2y(ly1*lely)
      real w1z(lz1*lelz), w2z(lz1*lelz)

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

      logical ifverbose
      integer icalld
      save    icalld
      data    icalld  /0/

      integer icall2
      save    icall2
      data    icall2 /-9/

      real atime,timel,times

      integer indts, ind_1, ind_2, nrec, ss
      save    indts, ind_1, ind_2, nrec, ss
      save    times
      save    domain_x, domain_y, domain_z

      character*80 pippo
      character*80 val1, val2, val3, val4, val5, val6
      character*80 val7, val8, val9, val10
      character*80 inputname1, inputname2, inputname3

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C
C     Note that planar_r and s routines must be adjusted for nelx/2+1 and nely/2+1
C     if and only  if you've changed your grid in the x-y plane!
C     Remember to change the number of nelx and nely in the size file correctly
      nelx = 16       ! Number of elements in x,y, and z directions.
      nely = 16       ! NOTE, this may vary from one mesh to the next.
      nelz = 62      ! First adjust this based on your grid then do the same for plan_r and s
      ntot = nx1*ny1*nz1*nelv 
      mid_r = nx1*nelx / 2 + 1
      mid_s = ny1*nely / 2 + 1

      nto2    = nx2*ny2*nz2*nelv
      ntot_2d = kx1*nx1*ky1*ny1
c     checking if the element number in the x and y directions are even
      if (istep.eq.0) then
       if(mod(nely,2).ne.0.or.mod(nelx,2).ne.0) then
        if(nid.eq.0) write(*,*) 'ABORT: nelx and nely have to be even!'
        call exitt
       endif
       if(nelx.gt.lelx .or. nely.gt.lely .or. nelz.gt.lelz) then
        if(nid.eq.0) write(*,*) 'ABORT: nel_xyz > lel_xyz!'
        call exitt
       endif
      endif


      if (istep.ne.icall2) then
        call invers2(jacmi,jacm1,ntot)
        icall2=istep
      endif

      if (kx1.eq.1) then
        write(6,*) nid,
     $    'Error. Uncomment kx1 param. stmt. in avg_all, navier5.f'
        return
      endif

      if (icalld.eq.0) then
        icalld = icalld + 1
 
        atime  = 0.
        timel  = time

        xlmin = glmin(xm1,ntot)
        xlmax = glmax(xm1,ntot)
        domain_x = xlmax - xlmin

        ylmin = glmin(ym1,ntot)          
        ylmax = glmax(ym1,ntot)
        domain_y = ylmax - ylmin

        zlmin = glmin(zm1,ntot)          
        zlmax = glmax(zm1,ntot)
        domain_z = zlmax - zlmin

        call rzero(stat,ntot*nstat)
        call rzero(stat_xy,ntot_2d*nstat)
      endif

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

      dtime = time  - timel

      iastep = param(68) 
      if  (iastep.eq.0) iastep=param(15)   ! same as iostep
      if  (iastep.eq.0) iastep=500

      ifverbose=.false.
      if (istep.le.10) ifverbose=.false.
      if  (mod(istep,iastep).eq.0) ifverbose=.false. ! get rid of output for now

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

      if (istep.eq.1) then
        nrec  = 0
        ind_1 = nv
        ind_2 = n2ptc 
        ss    = 1
      endif

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

      if (mod(istep,ind_2).eq.0.and.istep.ge.1) then

        ind_2 = ind_2 + n2ptc

        param(66)=6.
        call outpost(vx,vy,vz,p0,p0,'2pt')
        param(66)=6.

      endif

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C  
     
      if (mod(istep,ind_1).eq.0.and.istep.ge.1) then
      
        if (ss.eq.1) then
            times = time
            ss    = 2
        endif

        nrec  = nrec + 1
        ind_1 = ind_1 + nv
        atime = atime + dtime
        
        call mappr(pm1,pr,wk1,wk2) ! map pressure to mesh 1 (vel. mesh) 
        call copy(p0,pm1,ntot)
        
        pmean = -surf_mean(pm1,1,'W  ',ierr)
        call cadd(p0,pmean,ntot)
        
        call comp_derivat(duidxj,vx,vy,vz,ur,us,ut,vr,vs,vt,wr,ws,wt)

        if (atime.ne.0..and.dtime.ne.0.) then
          beta  = dtime/atime
          alpha = 1.-beta

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C
C ---------------------------------------------------------------------- C
C ------  Write time history ------------------------------------------- C
C ---------------------------------------------------------------------- C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

          call planar_average_r(wavg_pl_r,vz,rw1,rw2)
          call planar_average_s(wavg_pl_s,vz,sw1,sw2)
c      dwdzmean = surf_mean(duidxj(:,1,3),1,'W  ',ierr)

          if (nid.eq.0) then	 	
           open(unit=37,access='append',file='history')
           if (ss.eq.2) then
            write(37,'(A)') 
     +       '  Time   <dp/dz>  <Re_bulk_mid_X>   <Re_bulk_mid_X>'
            ss = 3
           endif
          write(37,3) time,scale_vf(3),wavg_pl_r(mid_r),wavg_pl_s(mid_s)
          endif
          if (nid.eq.0.and.istep.eq.NTSTEPS) then
            close(37)
          endif

    3     format(1p15e17.9)

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C 

      call avg1(stat(1,1),vx,alpha,beta,ntot,'velx',ifverbose)     ! <u>
      call avg1(stat(1,2),vy,alpha,beta,ntot,'vely',ifverbose)     ! <v>
      call avg1(stat(1,3),vz,alpha,beta,ntot,'velz',ifverbose)     ! <w>
      call avg1(stat(1,4),p0,alpha,beta,ntot,'pres',ifverbose)     ! <p>
      
      call avg2(stat(1,5),vx,alpha,beta,ntot,'urms',ifverbose)     ! <uu> (u: instantaneous) 
      call avg2(stat(1,6),vy,alpha,beta,ntot,'vrms',ifverbose)     ! <vv> (v: instantaneous)
      call avg2(stat(1,7),vz,alpha,beta,ntot,'wrms',ifverbose)     ! <ww> (w: instantaneous)
      call avg2(stat(1,8),p0,alpha,beta,ntot,'prms',ifverbose)     ! <pp> (p: instantaneous)

      call avg3(stat(1,9),vx,vy,alpha,beta,ntot,'uvrm',ifverbose)  ! <uv> (u, v: instantaneous)
      call avg3(stat(1,10),vy,vz,alpha,beta,ntot,'vwrm',ifverbose) ! <vw> (v, w: instantaneous)
      call avg3(stat(1,11),vz,vx,alpha,beta,ntot,'wurm',ifverbose) ! <uw> (u, w: instantaneous)

      call avg4(stat(1,12),vx,p0,alpha,beta,ntot,'pu',ifverbose)   ! <pu> (p, u: instantaneous)
      call avg4(stat(1,13),vy,p0,alpha,beta,ntot,'pv',ifverbose)   ! <pv> (p, v: instantaneous)           
      call avg4(stat(1,14),vz,p0,alpha,beta,ntot,'pw',ifverbose)   ! <pw> (p, w: instantaneous)

      call avg5(stat(1,15),p0,duidxj,alpha,beta,ntot,              ! <pdudx> (p, dudx: instantaneous) 
     & 'pux',ifverbose)
      call avg5(stat(1,16),p0,duidxj,alpha,beta,ntot,           ! <pdudy> (p, dudx: instantaneous) 
     & 'puy',ifverbose)
      call avg5(stat(1,17),p0,duidxj,alpha,beta,ntot,           ! <pdudz> (p, dudx: instantaneous)
     & 'puz',ifverbose)

      call avg5(stat(1,18),p0,duidxj,alpha,beta,ntot,      ! <pdvdx> (p, dvdx: instantaneous) 
     & 'pvx',ifverbose)        
      call avg5(stat(1,19),p0,duidxj,alpha,beta,ntot,      ! <pdvdy> (p, dvdy: instantaneous)
     &           'pvy',ifverbose)        
      call avg5(stat(1,20),p0,duidxj,alpha,beta,ntot,              ! <pdvdz> (p, dudz: instantaneous)   
     & 'pvz',ifverbose)        

      call avg5(stat(1,21),p0,duidxj,alpha,beta,ntot,              ! <pdwdx> (p, dwdx: instantaneous) 
     & 'pwx',ifverbose)        
      call avg5(stat(1,22),p0,duidxj,alpha,beta,ntot,      ! <pdwdy> (p, dwdy: instantaneous)
     & 'pwy',ifverbose)        
      call avg5(stat(1,23),p0,duidxj,alpha,beta,ntot,     ! <pdwdz> (p, dwdz: instantaneous)
     & 'pwz',ifverbose)        

      call avg6(stat(1,24),vx,vx,vx,alpha,beta,ntot,'u3',ifverbose)    ! <uuu> (u: instantaneous) 
      call avg6(stat(1,25),vy,vy,vy,alpha,beta,ntot,'v3',ifverbose)    ! <vvv> (v: instantaneous) 
      call avg6(stat(1,26),vz,vz,vz,alpha,beta,ntot,'w3',ifverbose)    ! <www> (w: instantaneous) 
      call avg6(stat(1,27),p0,p0,p0,alpha,beta,ntot,'p3',ifverbose)    ! <ppp> (p: instantaneous) 

      call avg6(stat(1,28),vx,vx,vy,alpha,beta,ntot,'u2v',ifverbose)   ! <uuv> (u, v: instantaneous) 
      call avg6(stat(1,29),vx,vx,vz,alpha,beta,ntot,'u2w',ifverbose)   ! <uuw> (u, w: instantaneous)  
      call avg6(stat(1,30),vy,vy,vx,alpha,beta,ntot,'v2v',ifverbose)   ! <vvu> (v, u: instantaneous)	
      call avg6(stat(1,31),vy,vy,vz,alpha,beta,ntot,'v2w',ifverbose)   ! <vvw> (v, w: instantaneous) 

      call avg6(stat(1,32),vz,vz,vx,alpha,beta,ntot,'w2u',ifverbose)   ! <wwu> (w, u: instantaneous)
      call avg6(stat(1,33),vz,vz,vy,alpha,beta,ntot,'w2v',ifverbose)   ! <wwv> (w, v: instantaneous)
      call avg6(stat(1,34),vx,vy,vz,alpha,beta,ntot,'uvw',ifverbose)   ! <uvw> (u, v, w: instantaneous) 	

      call avg8(stat(1,35),vx,vx,vx,vx,alpha,beta,ntot,'u4',      ! <uuuu> (u: instantaneous)     
     & ifverbose)     
      call avg8(stat(1,36),vy,vy,vy,vy,alpha,beta,ntot,'v4',      ! <vvvv> (v: instantaneous)
     & ifverbose)     
      call avg8(stat(1,37),vz,vz,vz,vz,alpha,beta,ntot,'w4',      ! <wwww> (w: instantaneous)	 
     & ifverbose)     
      call avg8 (stat(1,38),p0,p0,p0,p0,alpha,beta,ntot,	      ! <pppp> (p: instantaneous) 
     & 'p4',ifverbose)

      call avg8(stat(1,39),vx,vx,vx,vy,alpha,beta,ntot,'u4',      ! <uuuv> (u: instantaneous)     
     & ifverbose)     
      call avg8(stat(1,40),vx,vx,vy,vy,alpha,beta,ntot,'v4',      ! <uuvv> (v: instantaneous)
     & ifverbose)     
      call avg8(stat(1,41),vx,vy,vy,vy,alpha,beta,ntot,'w4',      ! <uvvv> (w: instantaneous)	 
     & ifverbose)  

      call avg7(stat(1,42),duidxj,alpha,beta,ntot,'e11',ifverbose)   ! e11: <d2u2/dx2 + d2u2/dy2 + d2u2/dz2> (u: instantaneous)
      call avg7(stat(1,43),duidxj,alpha,beta,ntot,'e22',ifverbose)   ! e22: <d2v2/dx2 + d2v2/dy2 + d2v2/dz2> (v: instantaneous) 
      call avg7(stat(1,44),duidxj,alpha,beta,ntot,'e33',ifverbose)   ! e33: <d2w2/dx2 + d2w2/dy2 + d2w2/dz2> (w: instantaneous)

      call avg7(stat(1,45),duidxj,alpha,beta,ntot,'e12',ifverbose)   ! e12: <du/dx.dv/dx + du/dy.dv/dy + du/dz.dv/dz> (u, v: instantaneous)       
      call avg7(stat(1,46),duidxj,alpha,beta,ntot,'e13',ifverbose)   ! e13: <du/dx.dw/dx + du/dy.dw/dy + du/dz.dw/dz> (u, w: instantaneous) 
      call avg7(stat(1,47),duidxj,alpha,beta,ntot,'e23',ifverbose)   ! e23: <dv/dx.dw/dx + dv/dy.dw/dy + dv/dz.dw/dz> (v, w: instantaneous) 

      call avg5(stat(1,48),p0,duidxj,alpha,beta,ntot,'omz',ifverbose) ! <omz>
      call avg5(stat(1,49),p0,duidxj,alpha,beta,ntot,'ozz',ifverbose) ! <omz*omz>

      call avg5(stat(1,50),p0,duidxj,alpha,beta,ntot,'aaa',ifverbose) ! <dw/dx*dw/dx>
      call avg5(stat(1,51),p0,duidxj,alpha,beta,ntot,'bbb',ifverbose) ! <dw/dy*dw/dy>
      call avg5(stat(1,52),p0,duidxj,alpha,beta,ntot,'ccc',ifverbose) ! <dw/dx*dw/dy>

      call avg5(stat(1,53),p0,duidxj,alpha,beta,ntot,'ddd',ifverbose) ! <du/dx*du/dx>
      call avg5(stat(1,54),p0,duidxj,alpha,beta,ntot,'eee',ifverbose) ! <du/dy*du/dy>
      call avg5(stat(1,55),p0,duidxj,alpha,beta,ntot,'fff',ifverbose) ! <du/dx*du/dy>
     
      call avg5(stat(1,56),p0,duidxj,alpha,beta,ntot,'ggg',ifverbose) ! <dv/dx*dv/dx>
      call avg5(stat(1,57),p0,duidxj,alpha,beta,ntot,'hhh',ifverbose) ! <dv/dy*dv/dy>
      call avg5(stat(1,58),p0,duidxj,alpha,beta,ntot,'iii',ifverbose) ! <dv/dx*dv/dy>
     
      call avg5(stat(1,59),p0,duidxj,alpha,beta,ntot,'jjj',ifverbose) ! <du/dx*dv/dx>
      call avg5(stat(1,60),p0,duidxj,alpha,beta,ntot,'kkk',ifverbose) ! <du/dy*dv/dy>
      call avg5(stat(1,61),p0,duidxj,alpha,beta,ntot,'lll',ifverbose) ! <du/dx*dv/dy>
      call avg5(stat(1,62),p0,duidxj,alpha,beta,ntot,'mmm',ifverbose) ! <du/dy*dv/dx>

      call avg5(stat(1,63),p0,duidxj,alpha,beta,ntot,'nnn',ifverbose) ! <du/dx>
      call avg5(stat(1,64),p0,duidxj,alpha,beta,ntot,'ooo',ifverbose) ! <du/dy>
      call avg5(stat(1,65),p0,duidxj,alpha,beta,ntot,'ppp',ifverbose) ! <du/dz>

      call avg5(stat(1,66),p0,duidxj,alpha,beta,ntot,'qqq',ifverbose) ! <dv/dx>
      call avg5(stat(1,67),p0,duidxj,alpha,beta,ntot,'rrr',ifverbose) ! <dv/dy>
      call avg5(stat(1,68),p0,duidxj,alpha,beta,ntot,'sss',ifverbose) ! <dv/dz>

      call avg5(stat(1,69),p0,duidxj,alpha,beta,ntot,'ttt',ifverbose) ! <dw/dx>
      call avg5(stat(1,70),p0,duidxj,alpha,beta,ntot,'uuu',ifverbose) ! <dw/dy>
      call avg5(stat(1,71),p0,duidxj,alpha,beta,ntot,'vvv',ifverbose) ! <dw/dz>

          endif
      endif

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C 
C ------------------------------------------------------------------------------- C
C ------- average the statistical quantities in the homogeneous directions ------ C
C   planar_average_z; average in the homogeneous streamwise (axial) z-direction   C
C ------------------------------------------------------------------------------- C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C 

      if (mod(istep,iastep).eq.0.and.istep.gt.1) then

        call z_average(stat_xy(1,1),stat(1,1),w1,w2)    ! <u>
        call z_average(stat_xy(1,2),stat(1,2),w1,w2)    ! <v>
        call z_average(stat_xy(1,3),stat(1,3),w1,w2)    ! <w>
        call z_average(stat_xy(1,4),stat(1,4),w1,w2)    ! <p>

        call z_average(stat_xy(1,5),stat(1,5),w1,w2)    ! <uu>  (u: instantaneous)
        call z_average(stat_xy(1,6),stat(1,6),w1,w2)    ! <vv>  (v: instantaneous)
        call z_average(stat_xy(1,7),stat(1,7),w1,w2)    ! <ww>  (w: instantaneous)
        call z_average(stat_xy(1,8),stat(1,8),w1,w2)    ! <pp>  (p: instantaneous)
        
        call z_average(stat_xy(1,9),stat(1,9),w1,w2)    ! <uv>  (u, v: instantaneous)
        call z_average(stat_xy(1,10),stat(1,10),w1,w2)  ! <vw>  (v, w: instantaneous)
        call z_average(stat_xy(1,11),stat(1,11),w1,w2)  ! <uw>  (u, w: instantaneous)

        call z_average(stat_xy(1,12),stat(1,12),w1,w2)  ! <pu>  (p, u: instantaneous)
        call z_average(stat_xy(1,13),stat(1,13),w1,w2)  ! <pv>  (p, v: instantaneous)
        call z_average(stat_xy(1,14),stat(1,14),w1,w2)  ! <pw>  (p, w: instantaneous)
       
        call z_average(stat_xy(1,15),stat(1,15),w1,w2)  ! <pdudx>  (p, dudx: instantaneous)
        call z_average(stat_xy(1,16),stat(1,16),w1,w2)  ! <pdudy>  (p, dudy: instantaneous)
        call z_average(stat_xy(1,17),stat(1,17),w1,w2)  ! <pdudz>  (p, dudz: instantaneous)
       
        call z_average(stat_xy(1,18),stat(1,18),w1,w2)  ! <pdvdx>  (p, dvdx: instantaneous)
        call z_average(stat_xy(1,19),stat(1,19),w1,w2)  ! <pdvdy>  (p, dvdy: instantaneous)
        call z_average(stat_xy(1,20),stat(1,20),w1,w2)  ! <pdvdz>  (p, dvdz: instantaneous)
       
        call z_average(stat_xy(1,21),stat(1,21),w1,w2)  ! <pdwdx>  (p, dwdx: instantaneous)
        call z_average(stat_xy(1,22),stat(1,22),w1,w2)  ! <pdwdy>  (p, dwdy: instantaneous)
        call z_average(stat_xy(1,23),stat(1,23),w1,w2)  ! <pdwdz>  (p, dwdz: instantaneous)
       
        call z_average(stat_xy(1,24),stat(1,24),w1,w2)  ! <uuu>  (u: instantaneous)   
        call z_average(stat_xy(1,25),stat(1,25),w1,w2)  ! <vvv>  (v: instantaneous) 
        call z_average(stat_xy(1,26),stat(1,26),w1,w2)  ! <www>  (w: instantaneous) 
        call z_average(stat_xy(1,27),stat(1,27),w1,w2)  ! <ppp>  (p: instantaneous) 
       
        call z_average(stat_xy(1,28),stat(1,28),w1,w2)  ! <uuv>  (u, v: instantaneous)   
        call z_average(stat_xy(1,29),stat(1,29),w1,w2)  ! <uuw>  (u, w: instantaneous) 
        call z_average(stat_xy(1,30),stat(1,30),w1,w2)  ! <vvu>  (v, u: instantaneous)   
        call z_average(stat_xy(1,31),stat(1,31),w1,w2)  ! <vvw>  (v, w: instantaneous)
        call z_average(stat_xy(1,32),stat(1,32),w1,w2)  ! <wwu>  (w, u: instantaneous)   
        call z_average(stat_xy(1,33),stat(1,33),w1,w2)  ! <wwv>  (w, v: instantaneous)
        call z_average(stat_xy(1,34),stat(1,34),w1,w2)  ! <uvw>  (u, v, w: instantaneous) 
       
        call z_average(stat_xy(1,35),stat(1,35),w1,w2)  ! <uuuu>  (u: instantaneous)   
        call z_average(stat_xy(1,36),stat(1,36),w1,w2)  ! <vvvv>  (v: instantaneous) 
        call z_average(stat_xy(1,37),stat(1,37),w1,w2)  ! <wwww>  (w: instantaneous)   
        call z_average(stat_xy(1,38),stat(1,38),w1,w2)  ! <pppp>  (p: instantaneous)
       
        call z_average(stat_xy(1,39),stat(1,39),w1,w2)  ! <uuuv>  (u: instantaneous)   
        call z_average(stat_xy(1,40),stat(1,40),w1,w2)  ! <uuvv>  (v: instantaneous) 
        call z_average(stat_xy(1,41),stat(1,41),w1,w2)  ! <uvvv>  (w: instantaneous)
       
        call z_average(stat_xy(1,42),stat(1,42),w1,w2)  ! <e11>: <d2u2/dx2 + d2u2/dy2 + d2u2/dz2>  (u: instantaneous)
        call z_average(stat_xy(1,43),stat(1,43),w1,w2)  ! <e22>: <d2v2/dx2 + d2v2/dy2 + d2v2/dz2>  (v: instantaneous)
        call z_average(stat_xy(1,44),stat(1,44),w1,w2)  ! <e33>: <d2w2/dx2 + d2w2/dy2 + d2w2/dz2>  (w: instantaneous)
        call z_average(stat_xy(1,45),stat(1,45),w1,w2)  ! <e12>: <du/dx.dv/dx + du/dy.dv/dy + du/dz.dv/dz>  (u, v: instantaneous)
        call z_average(stat_xy(1,46),stat(1,46),w1,w2)  ! <e13>: <du/dx.dw/dx + du/dy.dw/dy + du/dz.dw/dz>  (u, w: instantaneous)
        call z_average(stat_xy(1,47),stat(1,47),w1,w2)  ! <e23>: <dv/dx.dw/dx + dv/dy.dw/dy + dv/dz.dw/dz>  (v, w: instantaneous)
       
        call z_average(stat_xy(1,48),stat(1,48),w1,w2)  ! <omz>
        call z_average(stat_xy(1,49),stat(1,49),w1,w2)  ! <omz*omz> 
       
        call z_average(stat_xy(1,50),stat(1,50),w1,w2)  ! <dw/dx*dw/dx>
        call z_average(stat_xy(1,51),stat(1,51),w1,w2)  ! <dw/dy*dw/dy>
        call z_average(stat_xy(1,52),stat(1,52),w1,w2)  ! <dw/dx*dw/dy>
       
        call z_average(stat_xy(1,53),stat(1,53),w1,w2)  ! <du/dx*du/dx>
        call z_average(stat_xy(1,54),stat(1,54),w1,w2)  ! <du/dy*du/dy>
        call z_average(stat_xy(1,55),stat(1,55),w1,w2)  ! <du/dx*du/dy>
       
        call z_average(stat_xy(1,56),stat(1,56),w1,w2)  ! <dv/dx*dv/dx>
        call z_average(stat_xy(1,57),stat(1,57),w1,w2)  ! <dv/dy*dv/dy>
        call z_average(stat_xy(1,58),stat(1,58),w1,w2)  ! <dv/dx*dv/dy>
       
        call z_average(stat_xy(1,59),stat(1,59),w1,w2)  ! <du/dx*dv/dx>
        call z_average(stat_xy(1,60),stat(1,60),w1,w2)  ! <du/dy*dv/dy>
        call z_average(stat_xy(1,61),stat(1,61),w1,w2)  ! <du/dx*dv/dy>
        call z_average(stat_xy(1,62),stat(1,62),w1,w2)  ! <du/dy*dv/dx>
       
        call z_average(stat_xy(1,63),stat(1,63),w1,w2)  ! <du/dx>
        call z_average(stat_xy(1,64),stat(1,64),w1,w2)  ! <du/dy>
        call z_average(stat_xy(1,65),stat(1,65),w1,w2)  ! <du/dz>
       
        call z_average(stat_xy(1,66),stat(1,66),w1,w2)  ! <dv/dx>
        call z_average(stat_xy(1,67),stat(1,67),w1,w2)  ! <dv/dy>
        call z_average(stat_xy(1,68),stat(1,68),w1,w2)  ! <dv/dz>
       
        call z_average(stat_xy(1,69),stat(1,69),w1,w2)  ! <dw/dx>
        call z_average(stat_xy(1,70),stat(1,70),w1,w2)  ! <dw/dy>
        call z_average(stat_xy(1,71),stat(1,71),w1,w2)  ! <dw/dz>

        atime = 0.
      endif

      timel = time

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C 
C ------------------------------------------------------------------------------- C
C ------------ Write statistics to file ----------------------------------------- C
C ------------------------------------------------------------------------------- C
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C 

      if(nid.eq.0.and.istep.eq.1) indts = 0
      iostep_avg = param(68)
      if(nid.eq.0 .and. istep.gt.0 .and. 
     &   mod(istep,iostep_avg).eq.0) then

        indts = indts + 1
        write(pippo,'(i4.4)') indts
CX          inputname1 = 'ZSTAT/stat'//trim(pippo)
CX	    inputname1 = 'stat'
        inputname1 = 'stat'//trim(pippo)


C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C
C ------------------------------------------------------------------------------- C
C ----- Inputname1 -------------------------------------------------------------- C
C ------------------------------------------------------------------------------- C	    
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

        open(unit=33,form='unformatted',file=inputname1)

        my=ny1*nely
        mx=nx1*nelx
        m=my*mx

        write(val1,'(1p15e17.9)') 1/param(2)  		     ! Reynolds number	  
        write(val2,'(1p15e17.9)') domain_x,domain_y,domain_z ! domain size
        write(val3,'(9i9)') nelx,nely,nelz                   ! number of elements 
        write(val4,'(9i9)') nx1-1,ny1-1,nz1-1                ! polynomial order
        write(val5,'(9i9)')       nstat                      ! number of saved statistics 
        write(val6,'(1p15e17.9)') times                      ! start time
        write(val7,'(1p15e17.9)') time + nv*DT               ! end time
        write(val8,'(1p15e17.9)') nrec*DT                    ! average time
        write(val9,'(1p15e17.9)') DT                         ! time step
        write(val10,'(9i9)')      nrec                       ! number of time records


        write(33) '(Re ='//trim(val1)
     &   //') (Lx, Ly, Lz ='//trim(val2)
     &   //') (nelx, nely, nelz ='//trim(val3)
     &   //') (Polynomial order ='//trim(val4)
     &   //') (Nstat ='//trim(val5)
     &   //') (start time ='//trim(val6)
     &   //') (end time ='//trim(val7)
     &   //') (average time ='//trim(val8)
     &   //') (time step ='//trim(val9)
     &   //') (nrec ='//trim(val10)
     &   //')'

            write(33) 1/param(2),
     &      domain_x, domain_y, domain_z,
     &      nelx    , nely    , nelz,
     &      nx1-1   , ny1-1   , nz1-1,
     &      nstat,
     &      times,
     &      time + nv*DT,
     &      nrec*DT,
     &      DT,
     &      nrec


        do i=1,nstat
          write(33) (stat_xy(j,i),j=1,m)  
        enddo

        close(33)
     
        nrec = 0
        times = time + nv*DT

      endif

      return
      end

C=======================================================================
      subroutine planar_average_r(ua,u,w1,w2)
c
c     Compute s-t planar average of quantity u()
c
      include 'SIZE'
      include 'GEOM'
      include 'PARALLEL'
      include 'WZ'
      include 'ZPER'
c
      real ua(nx1,nelx),u(nx1,ny1,nx1,nelv),w1(nx1,nelx),w2(nx1,nelx)
      integer e,eg,ex,ey,ez
c
      nx = nx1*nelx
      call rzero(ua,nx)
      call rzero(w1,nx)
c
      do e=1,nelt
c
         eg = lglel(e)
         call get_exyz(ex,ey,ez,eg,nelx,nely,nelz)
         if (ex.eq.9) then
         do i=1,1
         do k=1,nz1
         do j=1,ny1
            zz = (1.-zgm1(i,1))/2.  ! = 1 for i=1, = 0 for k=nx1
            aa = zz*area(j,k,4,e) + (1-zz)*area(j,k,2,e)  ! wgtd jacobian
            w1(i,ex) = w1(i,ex) + aa
            ua(i,ex) = ua(i,ex) + aa*u(i,j,k,e)
         enddo
         enddo
         enddo
         endif
      enddo
c
      call gop(ua,w2,'+  ',nx)
      call gop(w1,w2,'+  ',nx)
c
      do i=1,nx
         ua(i,1) = ua(i,1) / w1(i,1)   ! Normalize
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine planar_average_s(ua,u,w1,w2)
c
c     Compute r-t planar average of quantity u()
c
      include 'SIZE'
      include 'GEOM'
      include 'PARALLEL'
      include 'WZ'
      include 'ZPER'
c
      real ua(ny1,nely),u(nx1,ny1,nx1,nelv)
     $     ,w1(ny1,nely),w2(ny1,nely)
      integer e,eg,ex,ey,ez
c
      ny = ny1*nely
      call rzero(ua,ny)
      call rzero(w1,ny)
c
      do e=1,nelt
         eg = lglel(e)
         call get_exyz(ex,ey,ez,eg,nelx,nely,nelz)
         
         if (ey.eq.9) then
         do k=1,nz1
         do j=1,1
         do i=1,nx1
            zz = (1.-zgm1(j,2))/2.  ! = 1 for i=1, = 0 for k=nx1
            aa = zz*area(i,k,1,e) + (1-zz)*area(i,k,3,e)  ! wgtd jacobian
            w1(j,ey) = w1(j,ey) + aa
            ua(j,ey) = ua(j,ey) + aa*u(i,j,k,e)
         enddo
         enddo
         enddo
         endif
      enddo
c
      call gop(ua,w2,'+  ',ny)
      call gop(w1,w2,'+  ',ny)
c
      do i=1,ny
         ua(i,1) = ua(i,1) / w1(i,1)   ! Normalize
      enddo

      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      function my_surf_mean(u_prime,ifld,bc_in,ierr)

      include 'SIZE'
      include 'TOTAL'
      
      real u_prime(lx1*ly1*lz1,lelt,1:3*ldim)
      real u(1)

      integer e,f
      character*3 bc_in
      do kk=1,n
         u(kk) = u_prime(kk,1,3)
      end do

      usum = 0
      asum = 0

      nface = 2*ndim
      do e=1,nelv
      do f=1,nface
         if (cbc(f,e,ifld).eq.bc_in) then
            call fcsum2(usum_f,asum_f,u,e,f)
            usum = usum + usum_f
            asum = asum + asum_f
         endif
      enddo
      enddo

      usum = glsum(usum,1)  ! sum across processors
      asum = glsum(asum,1)

      surf_mean = usum
      ierr      = 1

      if (asum.gt.0) surf_mean = usum/asum
      if (asum.gt.0) ierr      = 0

      return
      end

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C
C										  C
C     This routine to modify element vertices					  C	
C     Adam Peplinski; 	     	     						  C	
C     This is the first user subroutine called; I use it to			  C
C     set my own parameters not to mess with .rea file.	 			  C
C     I open .rea.urs file and read all the parameters				  C
C										  C	
C     Parameters used in this code						  C
C     UPARAM(1) - /= 0 if restart						  C
C     UPARAM(2) - checkpiont dump frequency (like IOSTEP)			  C
C     UPARAM(3) - 0 create restart file names using SESSION.restart		  C
C                 1 use standard file names rs8restart?....			  C
C     UPARAM(4) - history pionts frequency  					  C
C     UPARAM(5) - jet/flow velocity ratio					  C
C     UPARAM(6) - pipe radius	    						  C
C     UPARAM(7) - 1 to use vel. field from restart as basefield u,v; 		  C
C                   otherwise blasius 	   	      				  C
C     UPARAM(8) - position of delta*=1						  C
C										  C	
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C


      SUBROUTINE usrdat
      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'    ! For nelx,nely,nelz - needed for z_average

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

      integer i, len, ierr

      integer unparam
      parameter (unparam=20)
      real UPARAM(unparam)
      common /userparam/ UPARAM

      character*132 fname 
      character*1 fnam1(132)
      equivalence (fnam1,fname)

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

      call rzero(UPARAM,unparam)

c     Open parameter file and read contents
      ierr=0
      if (NID.eq.0) then
         call blank(fname,132)
         len = ltrunc(REAFLE,132)
         call chcopy(fnam1(1),REAFLE,len)
         call chcopy(fnam1(len+1),'.usr',4)
         write(6,*) 'Openning uresr parameter file: ',fname
         open (unit=59,file=fname,err=30, status='old')
         read(59,*,err=30)      ! skip header
         read(59,*,err=30) len  ! number of lines to read
         goto 31
 30      ierr=1
 31   endif
      call err_chk(ierr,'Error reading .rea.usr file.$')

c     send number of parameters
      call bcast(len ,ISIZE)

c     compare with array length
      if (len.gt.unparam) then
         if(NID.eq.0) write(6,*) 'ERROR: too many parameters in ',fname
         call exitt
      endif

c     read and distribute user parameters
      ierr=0
      if (NID.eq.0) then
         do i=1, len
            read(59,*,err=40) UPARAM(i)
            write(6,45) i,UPARAM(i)
         enddo
         close(59)
         goto 41
 40      ierr=1
 41   endif
      call err_chk(ierr,'Error reading .rea.usr file.$')

 45   FORMAT('UPARAM(',I2,') = ',G13.5)

      call bcast(UPARAM ,unparam*WDSIZE)

      nelx = 16
      nely = 16
      nelz = 62

      RETURN
      END
C======================================================================
C=======================================================================
c-----------------------------------------------------------------------
      subroutine avg1(avg,f,alpha,beta,n,name,ifverbose)
      include 'SIZE'
      include 'TSTEP'
c
      real avg(n),f(n)
      character*4 name
      logical ifverbose
c
      do k=1,n
         avg(k) = alpha*avg(k) + beta*f(k)
      enddo
c
      if (ifverbose) then
         avgmax = glmax(avg,n)
         avgmin = glmin(avg,n)
         if (nid.eq.0) write(6,1) istep,time,avgmin,avgmax
     $                           ,alpha,beta,name
    1    format(i9,1p5e13.5,1x,a4,' av1mnx')
      endif
c
      return
      end
c-----------------------------------------------------------------------

      subroutine avg2(avg,f,alpha,beta,n,name,ifverbose)
      include 'SIZE'
      include 'TSTEP'
c
      real avg(n),f(n)
      character*4 name
      logical ifverbose
c
      do k=1,n
         avg(k) = alpha*avg(k) + beta*f(k)*f(k)
      enddo
c
      if (ifverbose) then
         avgmax = glmax(avg,n)
         avgmin = glmin(avg,n)
         if (nid.eq.0) write(6,1) istep,time,avgmin,avgmax
     $                           ,alpha,beta,name
    1    format(i9,1p5e13.5,1x,a4,' av2mnx')
      endif
c
      return
      end
c-----------------------------------------------------------------------

      subroutine avg3(avg,f,g,alpha,beta,n,name,ifverbose)
      include 'SIZE'
      include 'TSTEP'
c
      real avg(n),f(n),g(n)
      character*4 name
      logical ifverbose
c
      do k=1,n
         avg(k) = alpha*avg(k) + beta*f(k)*g(k)
      enddo
c
      if (ifverbose) then
         avgmax = glmax(avg,n)
         avgmin = glmin(avg,n)
         if (nid.eq.0) write(6,1) istep,time,avgmin,avgmax
     $                           ,alpha,beta,name
    1    format(i9,1p5e13.5,1x,a4,' av3mnx')
      endif
c
      return
      end
c-----------------------------------------------------------------------

      subroutine avg4(avg,f,g,alpha,beta,n,name,ifverbose)
      include 'SIZE'
      include 'TSTEP'
      include 'SOLN'

C      real pm1(lx1*ly1*lz1*lelv)
C      real wk1(lx1*ly1*lz1)
C      real wk2(lx1*ly1*lz1)

      real avg(n),f(n),g(n)
      character*4 name
      logical ifverbose

C      call mappr(pm1,pr,wk1,wk2) ! map pressure to mesh 1 (vel. mesh) 
c
      do k=1,n
         avg(k) = alpha*avg(k) + beta*f(k)*g(k)
      enddo
c
      if (ifverbose) then
         avgmax = glmax(avg,n)
         avgmin = glmin(avg,n)
         if (nid.eq.0) write(6,1) istep,time,avgmin,avgmax
     $                           ,alpha,beta,name
    1    format(i9,1p5e13.5,1x,a4,' av3mnx')
      endif
c
      return
      end
c-----------------------------------------------------------------------

      subroutine avg5(avg,f,duidxj,alpha,beta,n,name,ifverbose)
c
      include 'SIZE'
      include 'TOTAL'
c
      real duidxj(lx1*ly1*lz1,lelt,1:3*ldim)
CX      real ur(lx1*ly1*lz1),us(lx1*ly1*lz1),ut(lx1*ly1*lz1)
CX      real vr(lx1*ly1*lz1),vs(lx1*ly1*lz1),vt(lx1*ly1*lz1)
CX      real wr(lx1*ly1*lz1),ws(lx1*ly1*lz1),wt(lx1*ly1*lz1)
c     
      integer e
      real avg(n),f(n) 
      character*3 name
      logical ifverbose
c      
CX      call comp_derivat(duidxj,vx,vy,vz,ur,us,ut,vr,vs,vt,wr,ws,wt)
c
      if (name .eq. 'pux') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*f(k)*duidxj(k,1,1)
         enddo
      elseif (name .eq. 'pvy') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*f(k)*duidxj(k,1,2)
         enddo
      elseif (name .eq. 'pwz') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*f(k)*duidxj(k,1,3)
         enddo
      elseif (name .eq. 'puy') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*f(k)*duidxj(k,1,4)
         enddo
      elseif (name .eq. 'pvz') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*f(k)*duidxj(k,1,5)
         enddo
      elseif (name .eq. 'pwx') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*f(k)*duidxj(k,1,6)
         enddo
      elseif (name .eq. 'puz') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*f(k)*duidxj(k,1,7)
         enddo
      elseif (name .eq. 'pvx') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*f(k)*duidxj(k,1,8)
         enddo
      elseif (name .eq. 'pwy') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*f(k)*duidxj(k,1,9)
         enddo
      elseif (name .eq. 'omz') then  ! <omz>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*(duidxj(k,1,8)-duidxj(k,1,4))
         enddo   
      elseif (name .eq. 'ozz') then  ! <omz*omz>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*
     $     (duidxj(k,1,8)-duidxj(k,1,4))**2
         enddo 
      elseif (name .eq. 'aaa') then  ! <dw/dx*dw/dx>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,6)*duidxj(k,1,6)
         enddo  
      elseif (name .eq. 'bbb') then  ! <dw/dy*dw/dy>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,9)*duidxj(k,1,9)
         enddo 
      elseif (name .eq. 'ccc') then  ! <dw/dx*dw/dy>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,6)*duidxj(k,1,9)
         enddo 
      elseif (name .eq. 'ddd') then  ! <du/dx*du/dx>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,1)*duidxj(k,1,1)
         enddo 
      elseif (name .eq. 'eee') then  ! <du/dy*du/dy>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,4)*duidxj(k,1,4)
         enddo    
      elseif (name .eq. 'fff') then  ! <du/dx*du/dy>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,1)*duidxj(k,1,4)
         enddo  
      elseif (name .eq. 'ggg') then  ! <dv/dx*dv/dx>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,8)*duidxj(k,1,8)
         enddo  
      elseif (name .eq. 'hhh') then  ! <dv/dy*dv/dy>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,2)*duidxj(k,1,2)
         enddo 
      elseif (name .eq. 'iii') then  ! <dv/dx*dv/dy>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,8)*duidxj(k,1,2)
         enddo 
      elseif (name .eq. 'jjj') then  ! <du/dx*dv/dx>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,1)*duidxj(k,1,8)
         enddo 
      elseif (name .eq. 'kkk') then  ! <du/dy*dv/dy>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,4)*duidxj(k,1,2)
         enddo 
      elseif (name .eq. 'lll') then  ! <du/dx*dv/dy>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,1)*duidxj(k,1,2)
         enddo 
      elseif (name .eq. 'mmm') then  ! <du/dy*dv/dx>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,4)*duidxj(k,1,8)
         enddo    
      elseif (name .eq. 'nnn') then  ! <du/dx>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,1)
         enddo  
      elseif (name .eq. 'ooo') then  ! <du/dy>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,4)
         enddo  
      elseif (name .eq. 'ppp') then  ! <du/dz>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,7)
         enddo  
      elseif (name .eq. 'qqq') then  ! <dv/dx>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,8)
         enddo  
      elseif (name .eq. 'rrr') then  ! <dv/dy>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,2)
         enddo  
      elseif (name .eq. 'sss') then  ! <dv/dz>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,5)
         enddo  
      elseif (name .eq. 'ttt') then  ! <dw/dx>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,6)
         enddo 
      elseif (name .eq. 'uuu') then  ! <dw/dy>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,9)
         enddo 
      elseif (name .eq. 'vvv') then  ! <dw/dz>
         do k=1,n
            avg(k) = alpha*avg(k) + beta*duidxj(k,1,3)
         enddo 
c
      endif
c
      if (ifverbose) then
         avgmax = glmax(avg,n)
         avgmin = glmin(avg,n)
         if (nid.eq.0) write(6,1) istep,time,avgmin,avgmax
     $                           ,alpha,beta,name
    1    format(i9,1p5e13.5,1x,a4,' av3mnx')
      endif
c
      return
      end
c----------------------------------------------------------------------

      subroutine avg6(avg,f,g,h,alpha,beta,n,name,ifverbose)
      include 'SIZE'
      include 'TSTEP'
c
      real avg(n),f(n),g(n),h(n)
      character*4 name
      logical ifverbose
c
      do k=1,n
         avg(k) = alpha*avg(k) + beta*f(k)*g(k)*h(k)
      enddo
c
      if (ifverbose) then
         avgmax = glmax(avg,n)
         avgmin = glmin(avg,n)
         if (nid.eq.0) write(6,1) istep,time,avgmin,avgmax
     $                           ,alpha,beta,name
    1    format(i9,1p5e13.5,1x,a4,' av3mnx')
      endif
c
      return
      end
c-----------------------------------------------------------------------

      subroutine avg7(avg,duidxj,alpha,beta,n,name,ifverbose)
c
      include 'SIZE'
      include 'TOTAL'
c
      real duidxj(lx1*ly1*lz1,lelt,1:3*ldim)
CX      real ur(lx1*ly1*lz1),us(lx1*ly1*lz1),ut(lx1*ly1*lz1)
CX      real vr(lx1*ly1*lz1),vs(lx1*ly1*lz1),vt(lx1*ly1*lz1)
CX      real wr(lx1*ly1*lz1),ws(lx1*ly1*lz1),wt(lx1*ly1*lz1)
      real avg(n)
      character*3 name
      logical ifverbose
c      

CX      call comp_derivat(duidxj,vx,vy,vz,ur,us,ut,vr,vs,vt,wr,ws,wt)
c
      if (name .eq. 'e11') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*(duidxj(k,1,1)*duidxj(k,1,1) +
     $           duidxj(k,1,4)*duidxj(k,1,4) + 
     $           duidxj(k,1,7)*duidxj(k,1,7))
         enddo
      elseif (name .eq. 'e22') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*(duidxj(k,1,8)*duidxj(k,1,8) +
     $           duidxj(k,1,2)*duidxj(k,1,2) + 
     $           duidxj(k,1,5)*duidxj(k,1,5))
         enddo
      elseif (name .eq. 'e33') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*(duidxj(k,1,6)*duidxj(k,1,6) +
     $           duidxj(k,1,9)*duidxj(k,1,9) + 
     $           duidxj(k,1,3)*duidxj(k,1,3))
         enddo
      elseif (name .eq. 'e12') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*(duidxj(k,1,1)*duidxj(k,1,8) +
     $           duidxj(k,1,4)*duidxj(k,1,2) + 
     $           duidxj(k,1,7)*duidxj(k,1,5))
         enddo
      elseif (name .eq. 'e13') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*(duidxj(k,1,1)*duidxj(k,1,6) +
     $           duidxj(k,1,4)*duidxj(k,1,9) + 
     $           duidxj(k,1,7)*duidxj(k,1,3))
         enddo
      elseif (name .eq. 'e23') then
         do k=1,n
            avg(k) = alpha*avg(k) + beta*(duidxj(k,1,8)*duidxj(k,1,6) +
     $           duidxj(k,1,2)*duidxj(k,1,9) + 
     $           duidxj(k,1,5)*duidxj(k,1,3))
         enddo
c
      endif
c
      if (ifverbose) then
         avgmax = glmax(avg,n)
         avgmin = glmin(avg,n)
         if (nid.eq.0) write(6,1) istep,time,avgmin,avgmax
     $                           ,alpha,beta,name
    1    format(i9,1p5e13.5,1x,a4,' av3mnx')
      endif
c
      return
      end
c-----------------------------------------------------------------------

      subroutine avg8(avg,f,g,h,s,alpha,beta,n,name,ifverbose)
      include 'SIZE'
      include 'TSTEP'
c
      real avg(n),f(n),g(n),h(n),s(n)
      character*4 name
      logical ifverbose
c
      do k=1,n
         avg(k) = alpha*avg(k) + beta*f(k)*g(k)*h(k)*s(k)
      enddo
c
      if (ifverbose) then
         avgmax = glmax(avg,n)
         avgmin = glmin(avg,n)
         if (nid.eq.0) write(6,1) istep,time,avgmin,avgmax
     $                           ,alpha,beta,name
    1    format(i9,1p5e13.5,1x,a4,' av3mnx')
      endif
c
      return
      end
c-----------------------------------------------------------------------

      subroutine comp_derivat(duidxj,u,v,w,ur,us,ut,vr,vs,vt,wr,ws,wt)
c
      include 'SIZE'
      include 'TOTAL'
c
      integer e
c
      real duidxj(lx1*ly1*lz1,lelt,1:3*ldim)    ! 9 terms 
      real u  (lx1*ly1*lz1,lelt)
      real v  (lx1*ly1*lz1,lelt)
      real w  (lx1*ly1*lz1,lelt)
      real ur (1) , us (1) , ut (1)
      real vr (1) , vs (1) , vt (1)
      real wr (1) , ws (1) , wt (1)
c
c      common /dudxyj/ jacmi(lx1*ly1*lz1,lelt)
c      real jacmi
c
      n    = nx1-1                          ! Polynomial degree
      nxyz = nx1*ny1*nz1
c
      do e=1,nelv
         call local_grad3(ur,us,ut,u,N,e,dxm1,dxtm1)
         call local_grad3(vr,vs,vt,v,N,e,dxm1,dxtm1)
         call local_grad3(wr,ws,wt,w,N,e,dxm1,dxtm1)
c

      do k=1,nxyz
         duidxj(k,e,1) = jacmi(k,e)*(ur(k)*rxm1(k,1,1,e)+
     $        us(k)*sxm1(k,1,1,e)+
     $        ut(k)*txm1(k,1,1,e))
         duidxj(k,e,2) = jacmi(k,e)*(vr(k)*rym1(k,1,1,e)+
     $        vs(k)*sym1(k,1,1,e)+
     $        vt(k)*tym1(k,1,1,e))
         duidxj(k,e,3) = jacmi(k,e)*(wr(k)*rzm1(k,1,1,e)+
     $        ws(k)*szm1(k,1,1,e)+
     $        wt(k)*tzm1(k,1,1,e))
         duidxj(k,e,4) = jacmi(k,e)*(ur(k)*rym1(k,1,1,e)+
     $        us(k)*sym1(k,1,1,e)+
     $        ut(k)*tym1(k,1,1,e))
         duidxj(k,e,5) = jacmi(k,e)*(vr(k)*rzm1(k,1,1,e)+
     $        vs(k)*szm1(k,1,1,e)+
     $        vt(k)*tzm1(k,1,1,e))
         duidxj(k,e,6) = jacmi(k,e)*(wr(k)*rxm1(k,1,1,e)+
     $        ws(k)*sxm1(k,1,1,e)+
     $        wt(k)*txm1(k,1,1,e))
         duidxj(k,e,7) = jacmi(k,e)*(ur(k)*rzm1(k,1,1,e)+
     $        us(k)*szm1(k,1,1,e)+
     $        ut(k)*tzm1(k,1,1,e))
         duidxj(k,e,8) = jacmi(k,e)*(vr(k)*rxm1(k,1,1,e)+
     $        vs(k)*sxm1(k,1,1,e)+
     $        vt(k)*txm1(k,1,1,e))
         duidxj(k,e,9) = jacmi(k,e)*(wr(k)*rym1(k,1,1,e)+
     $        ws(k)*sym1(k,1,1,e)+
     $        wt(k)*tym1(k,1,1,e))
      enddo
      enddo
c
c
      return
      end
