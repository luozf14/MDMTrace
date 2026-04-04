subroutine mdmfm_init() bind(C, name="mdmfm_init")
  use iso_c_binding
  implicit none
  integer(c_int) :: flag
  real(c_double) :: thtspec, trgt1, am(4), qvalue, eexc, thetacal(10), ekine
  common /kineblck/ thtspec, trgt1, am, qvalue, eexc, thetacal, ekine
  external :: raytrace

  flag = 0
  call raytrace(flag)
  trgt1 = 0.0d0
end subroutine mdmfm_init

subroutine mdmfm_set_probes(dipole_probe, multipole_probe) bind(C, name="mdmfm_set_probes")
  use iso_c_binding
  implicit none
  real(c_double), value :: dipole_probe
  real(c_double), value :: multipole_probe
  real(c_double) :: data(75, 200), ititle(200)
  real(c_double) :: jeff_params(6)
  real(c_double) :: bqr, bhr, bor, bdr, bddr
  common /BLCK0/ data, ititle

  jeff_params = (/ -0.51927d0, 0.038638d0, 0.028404d0, -0.022797d0, -0.019275d0, 0.755583d0 /)

  bqr = -multipole_probe * 1.0d-4 * jeff_params(6)
  bhr = bqr * jeff_params(2) / jeff_params(1)
  bor = bqr * jeff_params(3) / jeff_params(1)
  bdr = bqr * jeff_params(4) / jeff_params(1)
  bddr = bqr * jeff_params(5) / jeff_params(1)

  data(14, 4) = bqr
  data(15, 4) = bhr
  data(16, 4) = bor
  data(17, 4) = bdr
  data(18, 4) = bddr
  data(15, 5) = dipole_probe * 1.034d0 * 1.0d-4
end subroutine mdmfm_set_probes

subroutine mdmfm_eval_entry_multipole_field(x_local, y_local, z_local, bx_out, by_out, bz_out) &
    bind(C, name="mdmfm_eval_entry_multipole_field")
  use iso_c_binding
  implicit none
  real(c_double), value :: x_local, y_local, z_local
  real(c_double), intent(out) :: bx_out, by_out, bz_out
  real(c_double) :: data(75, 200), ititle(200)
  real(c_double) :: bx, by, bz, k, tc(6), dtc(6)
  real(c_double) :: d, s, bt, grad1, grad2, grad3, grad4, grad5
  real(c_double) :: c0, c1, c2, c3, c4, c5
  real(c_double) :: dh, do_val, dd, ddd, dsh, dso, dsd, dsdd
  real(c_double) :: l, rad, bqd, bhx, boc, bdc, bdd
  real(c_double) :: z11, z12, z21, z22
  real(c_double) :: frh, fro, frd, frdd
  real(c_double) :: z_enter_max, z_exit_min
  integer(c_int) :: in_region
  common /BLCK0/ data, ititle
  common /BLCK10/ bx, by, bz, k, tc, dtc
  common /BLCK90/ d, s, bt, grad1, grad2, grad3, grad4, grad5
  common /BLCK91/ c0, c1, c2, c3, c4, c5
  common /BLCK92/ in_region
  common /BLCK93/ dh, do_val, dd, ddd, dsh, dso, dsd, dsdd
  external :: bpoles

  l = data(12, 4)
  rad = data(13, 4)
  bqd = data(14, 4)
  bhx = data(15, 4)
  boc = data(16, 4)
  bdc = data(17, 4)
  bdd = data(18, 4)
  z11 = data(19, 4)
  z12 = data(20, 4)
  z21 = data(21, 4)
  z22 = data(22, 4)
  frh = data(35, 4)
  fro = data(36, 4)
  frd = data(37, 4)
  frdd = data(38, 4)
  dsh = data(39, 4)
  dso = data(40, 4)
  dsd = data(41, 4)
  dsdd = data(42, 4)

  if (frh == 0.0d0) frh = 1.0d0
  if (fro == 0.0d0) fro = 1.0d0
  if (frd == 0.0d0) frd = 1.0d0
  if (frdd == 0.0d0) frdd = 1.0d0

  z_enter_max = -0.5d0 * l - z12
  z_exit_min = 0.5d0 * l + z21

  if ((x_local ** 2 + y_local ** 2) > rad ** 2) then
    bx_out = 0.0d0
    by_out = 0.0d0
    bz_out = 0.0d0
    return
  end if

  if (z_local < (-0.5d0 * l - z11) .or. z_local > (0.5d0 * l + z22)) then
    bx_out = 0.0d0
    by_out = 0.0d0
    bz_out = 0.0d0
    return
  end if

  d = 2.0d0 * rad
  dh = frh * d
  do_val = fro * d
  dd = frd * d
  ddd = frdd * d
  grad1 = -bqd / rad
  grad2 = bhx / (rad ** 2)
  grad3 = -boc / (rad ** 3)
  grad4 = bdc / (rad ** 4)
  grad5 = -bdd / (rad ** 5)
  k = 0.0d0
  dtc = 0.0d0

  if (z_local <= z_enter_max) then
    in_region = 1
    c0 = data(23, 4)
    c1 = data(24, 4)
    c2 = data(25, 4)
    c3 = data(26, 4)
    c4 = data(27, 4)
    c5 = data(28, 4)
    tc(1) = -x_local
    tc(2) = y_local
    tc(3) = -(z_local + 0.5d0 * l)
  else if (z_local < z_exit_min) then
    in_region = 2
    grad1 = -grad1
    grad3 = -grad3
    grad5 = -grad5
    tc(1) = x_local
    tc(2) = y_local
    tc(3) = z_local - 0.5d0 * l
  else
    in_region = 3
    grad1 = -grad1
    grad3 = -grad3
    grad5 = -grad5
    c0 = data(29, 4)
    c1 = data(30, 4)
    c2 = data(31, 4)
    c3 = data(32, 4)
    c4 = data(33, 4)
    c5 = data(34, 4)
    tc(1) = x_local
    tc(2) = y_local
    tc(3) = z_local - 0.5d0 * l
  end if

  call bpoles
  if (in_region == 1) then
    bx_out = -bx
    by_out = by
    bz_out = -bz
  else
    bx_out = bx
    by_out = by
    bz_out = bz
  end if
end subroutine mdmfm_eval_entry_multipole_field

subroutine mdmfm_eval_dipole_field(x_local, y_local, z_local, bx_out, by_out, bz_out) &
    bind(C, name="mdmfm_eval_dipole_field")
  use iso_c_binding
  implicit none
  real(c_double), value :: x_local, y_local, z_local
  real(c_double), intent(out) :: bx_out, by_out, bz_out
  real(c_double) :: data(75, 200), ititle(200)
  real(c_double) :: bx, by, bz, k, tc(6), dtc(6)
  real(c_double) :: ndx, bet1, gama, delt, csc
  real(c_double) :: rca, dels, br, s2, s3, s4, s5, s6, s7, s8, scor
  real(c_double) :: d, dg, s, bf, bt
  real(c_double) :: c0, c1, c2, c3, c4, c5
  real(c_double) :: rb, xc, zc
  integer(c_int) :: in_region, mtyp
  real(c_double) :: phi, alpha, beta
  real(c_double) :: z11, z12, z21, z22
  real(c_double) :: cosa, sina, copab, sipab, cospb, sinpb, sip2
  real(c_double) :: xb, zb, xc_region, zc_region
  real(c_double) :: radius, dr, theta, phi_rad
  real(c_double) :: cos_local_to_c, sin_local_to_c
  real(c_double), parameter :: deg_to_rad = 3.14159265358979323846d0 / 180.0d0
  real(c_double), parameter :: strip_half_width = 30.0d0
  common /BLCK0/ data, ititle
  common /BLCK10/ bx, by, bz, k, tc, dtc
  common /BLCK20/ ndx, bet1, gama, delt, csc
  common /BLCK21/ rca, dels, br, s2, s3, s4, s5, s6, s7, s8, scor
  common /BLCK22/ d, dg, s, bf, bt
  common /BLCK23/ c0, c1, c2, c3, c4, c5
  common /BLCK24/ rb, xc, zc
  common /BLCK25/ in_region, mtyp
  external :: bdip

  d = data(13, 5)
  dg = data(4, 5)
  mtyp = int(data(5, 5))
  if (mtyp == 0) mtyp = 1
  rb = data(14, 5)
  bf = data(15, 5)
  phi = data(16, 5)
  alpha = data(17, 5)
  beta = data(18, 5)
  ndx = data(19, 5)
  bet1 = data(20, 5)
  gama = data(21, 5)
  delt = data(22, 5)
  z11 = data(25, 5)
  z12 = data(26, 5)
  z21 = data(27, 5)
  z22 = data(28, 5)

  if (abs(y_local) > 0.5d0 * d) then
    bx_out = 0.0d0
    by_out = 0.0d0
    bz_out = 0.0d0
    return
  end if

  radius = sqrt((x_local + rb) ** 2 + z_local ** 2)
  dr = radius - rb
  theta = atan2(z_local, x_local + rb)
  phi_rad = phi * deg_to_rad

  if (dr < -strip_half_width .or. dr > strip_half_width) then
    bx_out = 0.0d0
    by_out = 0.0d0
    bz_out = 0.0d0
    return
  end if

  cosa = cos(alpha * deg_to_rad)
  sina = sin(alpha * deg_to_rad)
  xb = -x_local * cosa - z_local * sina
  zb = x_local * sina - z_local * cosa

  copab = cos((phi - alpha - beta) * deg_to_rad)
  sipab = sin((phi - alpha - beta) * deg_to_rad)
  cospb = cos((0.5d0 * phi - beta) * deg_to_rad)
  sinpb = sin((0.5d0 * phi - beta) * deg_to_rad)
  sip2 = sin(0.5d0 * phi * deg_to_rad)
  cos_local_to_c = cos((phi - beta) * deg_to_rad)
  sin_local_to_c = sin((phi - beta) * deg_to_rad)

  zc_region = -zb * copab + xb * sipab - 2.0d0 * rb * sip2 * cospb
  xc_region = -zb * sipab - xb * copab - 2.0d0 * rb * sip2 * sinpb

  k = 0.0d0
  dtc = 0.0d0
  bt = 0.0d0
  s = 0.0d0

  if (zb >= z12 .and. zb <= z11) then
    in_region = 1
    br = data(41, 5)
    c0 = data(29, 5)
    c1 = data(30, 5)
    c2 = data(31, 5)
    c3 = data(32, 5)
    c4 = data(33, 5)
    c5 = data(34, 5)
    dels = data(45, 5)
    rca = data(47, 5)
    csc = cos(alpha * deg_to_rad)
    scor = data(49, 5)
    s2 = data(51, 5) / rb + rca / 2.0d0
    s3 = data(52, 5) / (rb ** 2)
    s4 = data(53, 5) / (rb ** 3) + (rca ** 3) / 8.0d0
    s5 = data(54, 5) / (rb ** 4)
    s6 = data(55, 5) / (rb ** 5) + (rca ** 5) / 16.0d0
    s7 = data(56, 5) / (rb ** 6)
    s8 = data(57, 5) / (rb ** 7) + (rca ** 7) / 25.6d0
    xc = rb * cos(alpha * deg_to_rad)
    zc = -rb * sin(alpha * deg_to_rad)
    tc(1) = xb
    tc(2) = y_local
    tc(3) = zb
    call bdip
  else if (zc_region >= z21 .and. zc_region <= z22) then
    in_region = 3
    br = data(42, 5)
    c0 = data(35, 5)
    c1 = data(36, 5)
    c2 = data(37, 5)
    c3 = data(38, 5)
    c4 = data(39, 5)
    c5 = data(40, 5)
    dels = data(46, 5)
    rca = data(48, 5)
    csc = cos(beta * deg_to_rad)
    scor = data(50, 5)
    s2 = data(58, 5) / rb + rca / 2.0d0
    s3 = data(59, 5) / (rb ** 2)
    s4 = data(60, 5) / (rb ** 3) + (rca ** 3) / 8.0d0
    s5 = data(61, 5) / (rb ** 4)
    s6 = data(62, 5) / (rb ** 5) + (rca ** 5) / 16.0d0
    s7 = data(63, 5) / (rb ** 6)
    s8 = data(64, 5) / (rb ** 7) + (rca ** 7) / 25.6d0
    xc = -rb * cos(beta * deg_to_rad)
    zc = -rb * sin(beta * deg_to_rad)
    tc(1) = xc_region
    tc(2) = y_local
    tc(3) = zc_region
    call bdip
  else if (theta >= 0.0d0 .and. theta <= phi_rad) then
    in_region = 2
    xc = -rb
    zc = 0.0d0
    tc(1) = x_local
    tc(2) = y_local
    tc(3) = z_local
    call bdip
  else
    bx_out = 0.0d0
    by_out = 0.0d0
    bz_out = 0.0d0
    return
  end if

  if (in_region == 1) then
    bx_out = -cosa * bx + sina * bz
    by_out = by
    bz_out = -sina * bx - cosa * bz
  else
    bx_out = cos_local_to_c * bx - sin_local_to_c * bz
    by_out = by
    bz_out = sin_local_to_c * bx + cos_local_to_c * bz
  end if
end subroutine mdmfm_eval_dipole_field
