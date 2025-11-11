      ! Implementation of Continuum Species Transport (CST) for Nek5000
      !
      ! This module can be configured via userParam04 in your .par
      ! (parameter) file.
      !
      ! PARAMETER 04: Continuum Species Transport equation version.
      !   0 = off, passive scalar advection-diffusion only
      !   1 = Marschall (2012)
      !   2 = Li and Su (2025)
      ! Affects 2 terms, species_diff() called by uservp and
      ! species_src() called by userq.
      !
      ! PARAMETER 05: sink term mode for bubble interior.
      !   0 = off, no sink term
      !   1 = soft sink, add dc/dt = -sink_str*c explicit sink inside
      !       psi<0.1 with linear psi activation
      !   2 = hard sink, set c=0 inside psi<0.5
      !
      ! PARAMETER 06: source term mode for liquid bulk.
      !   0 = off, no source term
      !   1 = soft sink, add dc/dt = source_str*(1-c) implicit source
      !       inside psi>0.9 with linear psi activation
c-----------------------------------------------------------------------
      real function species_diff(ix,iy,iz,el)
      ! Calculate the coefficient for the diffusion term in the CST
      ! transport equation.
      implicit none
      integer ix, iy, iz, el
      include 'SIZE'
      include 'TOTAL'
      include 'CASE'

      integer species_equation_version
      real psi

      psi = t(ix,iy,iz,el,ifld_cls-1)
      psi = max(0.0,psi)
      psi = min(1.0,psi)
      species_equation_version = int(uparam(iprm_cst_ver))
      if (species_equation_version .eq. 0 .or.
     $    species_equation_version .eq. 1) then
        ! CST off or Marschall et al (2012) version
        species_diff = ((1.0-psi)*diffratio + psi)/Pe
      elseif (species_equation_version .eq. 2) then
        ! Li and Su (2025) version
        species_diff = (psi/(psi+(1.0-psi)*solubilityratio) +
     $      diffratio*(1.0-psi)/((1.0-psi)+psi/solubilityratio))/Pe
      else
        print *, "Bad species_equation_version",
     $        species_equation_version
        call abort
      endif

      endfunction
c-----------------------------------------------------------------------
      real function species_src(ix,iy,iz,el)
      ! Calculate the source terms in the CST transport equation.
      implicit none
      integer ix, iy, iz, el
      include 'SIZE'
      include 'TOTAL'
      include 'CASE'

      ! Save our computation result so we can compute once (per MPI
      ! rank per timestep) and look up the result afterwards for each
      ! GLL point (ix, iy, iz, el).
      common /speciestransport/ spdiv(lx1,ly1,lz1,lelv)
      real spdiv

      ! Parse user parameter.
      integer species_equation_version
      species_equation_version = int(uparam(iprm_cst_ver))
      if (species_equation_version .eq. 0) then
        species_src = 0.0
        return
      endif

      ! If ix = iy = iz = el = 1 (e.g. only on the first GLL point
      ! on the first local element), do all the work.
      if (ix*iy*iz*el .eq. 1) then
        block
          integer i, ntot
          real stmp(lx1,ly1,lz1,lelv), spx(lx1,ly1,lz1,lelv),
     $        spy(lx1,ly1,lz1,lelv), spz(lx1,ly1,lz1,lelv)

          ! For each local GLL point on each local element:
          ! (This seems to be just shorthand for iterating over all 4
          ! dimensions in 1 loop (due to array ordering).)
          ntot = lx1*ly1*lz1*nelt
          do i=1,ntot
            block
              real psi, term1, term2
              ! Scalar field value t(ix,iy,iz,el,ifield-1) for field
              ! number ifield.  So this is reading a length ntot array
              ! starting at ix=i, iy=1, iz=1, el=1 from CLS field.
              ! psi(:) = t(:, ifield=ifld_cls)
              psi = t(i,1,1,1,ifld_cls-1)
              ! Clip to [0, 1].
              psi = max(0.0,psi)
              psi = min(1.0,psi)
              ! Accumulation coefficient.
              if (species_equation_version .eq. 1) then
                ! Marschall (2012) version
                term1 = 1.0 - diffratio
              elseif (species_equation_version .eq. 2) then
                ! Li and Su (2025) version
                term1 = ((1.0 - diffratio)/solubilityratio)
     $                  / ((psi/solubilityratio+(1.0-psi))**2)
              else
                print *, "Bad species_equation_version",
     $                species_equation_version
                call abort
              endif
              ! CST coefficient.
              term2 = (1.0/solubilityratio - diffratio)
     $                / (psi/solubilityratio + (1.0-psi))
              stmp(i,1,1,1) = (term1 - term2)/Pe
            endblock
          enddo
          ! Multiply stmp by the current c field values.
          call col2(stmp, t(1,1,1,1,ifld_c-1),ntot)

          ! Multiply stmp by our version of "grad(alpha)".
          block
            real delta(lx1,ly1,lz1,lelv)
            ! Calculate surface area density from CLS, as delta.
            ! delta = 1/(4*epsilon) (cosh(clip(atanh(2*psi - 1))))^(-2)
            call deltals(t(1,1,1,1,ifld_cls-1),delta)
            ! Calculate interface unit normals from TLS, as sp{x,y,z}.
            call cls_normals(spx,spy,spz,ifld_tls)
            ! Multiply interface normal * surface area * (accumulation
            ! and CST coefficients), as sp{x,y,z}.
            call col4(spx,spx,delta,stmp,ntot)
            call col4(spy,spy,delta,stmp,ntot)
            call col4(spz,spz,delta,stmp,ntot)
          endblock

          ! Take the divergence of (spx, spy, spz), as spdiv.
          call opdiv(spdiv,spx,spy,spz)
          ! Multiply spdiv by QQ^T (scatter/gather operators) to
          ! convert from global to local coordinates.
          call dssum(spdiv,lx1,ly1,lz1)
          ! Multiply spdiv by binvm1, the "inverse mass matrix for
          ! velocity mesh", as spdiv.
          call col2(spdiv,binvm1,ntot)
        endblock
      endif

      ! For the rest of the elements/GLL points, just look up based on
      ! work already done.
      species_src = spdiv(ix,iy,iz,el)
      endfunction
c-----------------------------------------------------------------------
      subroutine add_artificial_srcsnk(ix,iy,iz,el,qvol,avol)
      ! Artificial sink term for the bubble interior and source term for
      ! the liquid bulk to drive quasi-steady-state MTC & Sherwood number
      ! equilibrium.
      implicit none
      integer, intent(in) :: ix, iy, iz, el
      real, intent(inout) :: qvol, avol
      include 'SIZE'
      include 'TOTAL'
      include 'CASE'

      real, external :: glsum
      real c, psi, volratio
      real sink, sink_dt
      common /speciestransport_sink/ sink, sink_dt
      integer sinkmode, sourcemode

      c = t(ix,iy,iz,el,ifld_c-1)
      psi = t(ix,iy,iz,el,ifld_cls-1)
      if (ix*iy*iz*el.eq.1) then
        sink_dt = sink_dt + dt
      endif

      ! Clip to [0, 1].
      psi = max(0.0, psi)
      psi = min(1.0, psi)

      ! Parse user parameters.
      sinkmode = int(uparam(iprm_sinkmode))
      sourcemode = int(uparam(iprm_sourcemode))

      ! Total term is of the form qvol - avol*c (so that avol can be
      ! specified implicitly to improve stability).
      !
      ! We need to be careful that neither term flips signs (for c<0 or
      ! c>1), so we gate the terms with if-statements on c. The gas/
      ! liquid are meant as sink/source and should never flip signs and
      ! become source/sink.
      !
      ! We also perform some thresholding of psi -> max(0.0, psi-0.9)/
      ! 0.1 and (1-psi) -> max(0.0, 0.1-psi)/0.1 to try to keep the
      ! source and sink terms out of the two-phase boundary region.
      ! Otherwise these terms will try to drive c=0.5 at psi=0.5,
      ! interfering with the mass transfer under study.

      if (sinkmode .eq. 2) then
        ! Hard sink -- set everything to zero within psi<0.5.
        if (psi .le. 0.5) then
          sink = sink + c*bm1(ix,iy,iz,el)
          c = 0.0
          t(ix,iy,iz,el,ifld_c-1) = 0.0
        endif
      elseif (sinkmode .eq. 1) then
        ! Add -c*(1-psi) term to drive bubble (psi=0) towards c=0. Do
        ! this as an explicit term; implicit sinks reduce stability.
        if (c .ge. 0) then
          sink = sink + (max(0.0, 0.1-psi)/0.1)*sink_str*c*
      $          bm1(ix,iy,iz,el)*dt
          qvol = qvol - (max(0.0, 0.1-psi)/0.1)*sink_str*c
        else
          ! Clip c field to avoid flipping sign of CST term. Count
          ! this as negative extraction to preserve conservation.
          sink = sink - (max(0.0, 0.1-psi)/0.1)*sink_str*c*
      $          bm1(ix,iy,iz,el)*dt
          c = 0.0
          t(ix,iy,iz,el,ifld_c-1) = 0.0
        endif
      endif

      if (sourcemode .eq. 1) then
        ! Add (1-c)*psi term to drive liquid bulk (psi=1) towards c=1.
        if (c .le. 1) then
          avol = avol + (max(0.0, psi-0.9)/0.1)*source_str
          qvol = qvol + (max(0.0, psi-0.9)/0.1)*source_str
        endif
      endif

      return
      end
c-----------------------------------------------------------------------
      real function get_sink()
      ! Return the rate of specie removed. This function has side
      ! effects (resets counters) and should not be called repeatedly.
      !
      ! This is called from userchk, gated on ifoutfld ("is this step
      ! a checkpoint?"), before add_artificial_srcsink is called from
      ! useric. So the value it returns is non-inclusive of the extract
      ! from the current (checkpoint) timestep. The extract from the
      ! current timestep will be included in the next checkpoint's
      ! get_sink() call.
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CASE'

      real, external :: glsum
      real sink, sink_dt
      common /speciestransport_sink/ sink, sink_dt
      ! Calculate sink rate and reset sink/sink_dt.
      get_sink = glsum(sink, 1)/sink_dt
      sink = 0.0
      sink_dt = 0.0
      endfunction
c-----------------------------------------------------------------------
      real function integrate_gradc()
      ! Calculate the integral of grad(c).dA over the interface defined
      ! by the level set field, and normalize it by total surface area.
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CASE'

      real, external :: glsum
      real total_area, surfaceintegral_gradc

      block
        real vec_x(lx1,ly1,lz1,lelv), vec_y(lx1,ly1,lz1,lelv),
     $      vec_z(lx1,ly1,lz1,lelv)
        integer ntot
        ntot = lx1*ly1*lz1*nelt

        ! Calculate interface unit normals from TLS.
        !   vec = nhat
        call cls_normals(vec_x, vec_y, vec_z, ifld_tls)

        ! Scale by surface area.
        !   vec = dA
        block
          real delta(lx1,ly1,lz1,lelv)
          ! Calculate surface area density from CLS, as delta.
          ! delta = 1/(4*epsilon) (cosh(clip(atanh(2*psi - 1))))^(-2)
          call deltals(t(1,1,1,1,ifld_cls-1), delta)
          total_area = glsum(delta, ntot)
          ! Multiply vector by surface area density.
          call col2(vec_x, delta, ntot)
          call col2(vec_y, delta, ntot)
          call col2(vec_z, delta, ntot)
        endblock

        ! Dot with grad(c).
        !   vec = grad(c) . dA
        block
          real gradc_x(lx1,ly1,lz1,lelv), gradc_y(lx1,ly1,lz1,lelv),
     $        gradc_z(lx1,ly1,lz1,lelv)
          ! Take grad(c) as gradc_{x,y,z}.
          call gradm1(gradc_x, gradc_y, gradc_z, t(1,1,1,1,ifld_c-1))
          ! Take inner product of vec by grad(c).
          call col2(vec_x, gradc_x, ntot)
          call col2(vec_y, gradc_y, ntot)
          call col2(vec_z, gradc_z, ntot)
        endblock

        ! Sum across elements.
        !  surfaceintegral_gradc = \int grad(c) . dA
        surfaceintegral_gradc = glsum(vec_x, ntot) + glsum(vec_y, ntot)
     $      + glsum(vec_z, ntot)
      endblock

      ! Return normalized grad(c) result.
      integrate_gradc = surfaceintegral_gradc/total_area
      endfunction
c-----------------------------------------------------------------------
      real function integrate_c_bulk()
      ! Calculate the integral of c dV over the liquid and normalize it
      ! by total volume.
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CASE'
      real, external :: glsc2, glsc3
      integer ntot
      real psi(lx1,ly1,lz1,nelt), total_volume, c_bulk_total
      ntot = lx1*ly1*lz1*nelt
      psi = t(1,1,1,1,ifld_cls-1)
      ! Clip to [0, 1].
      psi = max(0.0,psi)
      psi = min(1.0,psi)
      ! Apply same thresholding we apply in the source term.
      psi = max(0.0, psi-0.9)/0.1
      ! Perform volume integrals.
      c_bulk_total = glsc3(psi, t(1,1,1,1,ifld_c-1), bm1, ntot)
      total_volume = glsc2(psi, bm1, ntot)
      ! Compute quotient.
      integrate_c_bulk = c_bulk_total/total_volume
      endfunction
