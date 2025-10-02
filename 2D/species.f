      ! Implementation of Continuum Species Transport (CST) for Nek5000
      !
      ! This module can be configured via userParam04 in your .par
      ! (parameter) file.
      !
      ! PARAMETER 04: Continuum Species Transport equation version.
      !   0 = off, passive scalar advection-diffusion only
      !   1 = Marschall (2012)
      !   2 = Li and Su (2025)
      ! Affects 2 terms, speciesDiff() called by uservp and
      ! speciesSrc() called by userq.
      
      real function speciesDiff(ix,iy,iz,el)
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
      species_equation_version = int(uparam(iprm_cst_ver))
      if (species_equation_version .eq. 0 .or.
     $    species_equation_version .eq. 1) then
        ! CST off or Marschall et al (2012) version
        speciesDiff = ((1.0-psi)*diffratio + psi)/Pe
      elseif (species_equation_version .eq. 2) then
        ! Li and Su (2025) version
        speciesDiff = (psi/(psi+(1.0-psi)*solubilityratio) +
     $      diffratio*(1.0-psi)/((1.0-psi)+psi/solubilityratio))/Pe
      endif

      endfunction

      real function speciesSrc(ix,iy,iz,el)
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
        speciesSrc = 0.0
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
              endif
              ! CST coefficient.
              term2 = (1.0 - diffratio/solubilityratio)
     $                / (psi/solubilityratio + (1.0-psi))
              stmp(i,1,1,1) = term1 - term2
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
      speciesSrc = spdiv(ix,iy,iz,el)/Pe
      endfunction
