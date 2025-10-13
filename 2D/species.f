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
      ! PARAMETER 05: dc/dy source term coefficient for jump-periodic
      ! boundary condition.
      !   The total concentration added across the domain height will
      !   be dc/dy*(width*height*depth) for a simple box domain. This
      !   should be a negative value if the bubble is a source for c
      !   (so the jump-periodic BC is the "sink"), or a positive value
      !   if the bubble is sink (so the BC is "source"). Set to 0 to
      !   disable.
      !   In a non-dimensional case where c is in [0, 1] and the bubble
      !   is expected to absorb enough species to deplete the liquid
      !   around it faster than it rises, a good starting guess is
      !   1.0/domain_height so that a total jump of 1.0 is added.
      ! Note that using the jump periodic BC will cause the c field
      ! to no longer be directly interpretable as the concentration
      ! of a vertical pipe section, due to the linear (dc/dy)*y extra
      ! concentration we add as a source term. Thus we add a ctrue
      ! field (calculated only, not transported) to view the true
      ! concentration.
      
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

      real function jump_periodic_src(ix,iy,iz,el)
      ! Calculate the dc/dy source term for our jump-periodic boundary
      ! condition.
      implicit none
      integer ix, iy, iz, el
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'CASE'

      real dcdy
      dcdy = uparam(iprm_dcdy)
      jump_periodic_src = dcdy*y
      endfunction

      subroutine update_ctrue()
      ! Update ctrue (which is c, corrected for dcdy).
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CASE'

      real i, j, k, l, dcdy
      dcdy = uparam(iprm_dcdy)
      do i=1,lx1
        do j=1,ly1
          do k=1,lz1
            do l=1,nelt
              t(i,j,k,l,ifld_ctrue-1) = t(i,j,k,l,ifld_c-1) +
     $            dcdy*ym1(i,j,k,l)
            enddo
          enddo
        enddo
      enddo
      end

      subroutine species_sink()
      ! Set ctrue = 0 in the interior of the bubble.
      ! This translates to setting c = dcdy*y.
      ! Call this before update_ctrue() in usrchk()!
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'CASE'

      real i, j, k, l, ctrue, removed
      do i=1,lx1
        do j=1,ly1
          do k=1,lz1
            do l=1,nelt
              ! Could also blend between existing value and 0.0 by psi,
              ! rather than apply a sharp threshold at psi = 0.5.
              if (t(i,j,k,l,ifld_cls-1).lt.0.5) then
                ctrue = t(i,j,k,l,ifld_c-1) - jump_periodic_src(i,j,k,l)
                removed = removed + ctrue*binvm1(i,j,k,l)
                t(i,j,k,l,ifld_c-1) = jump_periodic_src(i,j,k,l)
              endif
            enddo
          enddo
        enddo
      enddo

      block
        real total_removed
        total_removed = glsum(removed, 1)
        if (nio.eq.0) then
          write(*,*) "Species sink:", total_removed
        endif
      endblock
      end