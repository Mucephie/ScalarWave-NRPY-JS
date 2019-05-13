{
/*
 * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
 */
 /*
  *  Original SymPy expressions:
  *  "[const double uu_dDD00 = invdx0**2*(-5269*uu/1800 + 5*uu_i0m1_i1_i2/3 - 5*uu_i0m2_i1_i2/21 + 5*uu_i0m3_i1_i2/126 - 5*uu_i0m4_i1_i2/1008 + uu_i0m5_i1_i2/3150 + 5*uu_i0p1_i1_i2/3 - 5*uu_i0p2_i1_i2/21 + 5*uu_i0p3_i1_i2/126 - 5*uu_i0p4_i1_i2/1008 + uu_i0p5_i1_i2/3150),
  *    const double uu_dDD11 = invdx1**2*(-5269*uu/1800 + 5*uu_i0_i1m1_i2/3 - 5*uu_i0_i1m2_i2/21 + 5*uu_i0_i1m3_i2/126 - 5*uu_i0_i1m4_i2/1008 + uu_i0_i1m5_i2/3150 + 5*uu_i0_i1p1_i2/3 - 5*uu_i0_i1p2_i2/21 + 5*uu_i0_i1p3_i2/126 - 5*uu_i0_i1p4_i2/1008 + uu_i0_i1p5_i2/3150),
  *    const double uu_dDD22 = invdx2**2*(-5269*uu/1800 + 5*uu_i0_i1_i2m1/3 - 5*uu_i0_i1_i2m2/21 + 5*uu_i0_i1_i2m3/126 - 5*uu_i0_i1_i2m4/1008 + uu_i0_i1_i2m5/3150 + 5*uu_i0_i1_i2p1/3 - 5*uu_i0_i1_i2p2/21 + 5*uu_i0_i1_i2p3/126 - 5*uu_i0_i1_i2p4/1008 + uu_i0_i1_i2p5/3150)]"
  */
  const double uu_i0_i1_i2m5 = in_gfs[IDX4(UUGF, i0, i1, i2 - 5)];
  const double uu_i0_i1_i2m4 = in_gfs[IDX4(UUGF, i0, i1, i2 - 4)];
  const double uu_i0_i1_i2m3 = in_gfs[IDX4(UUGF, i0, i1, i2 - 3)];
  const double uu_i0_i1_i2m2 = in_gfs[IDX4(UUGF, i0, i1, i2 - 2)];
  const double uu_i0_i1_i2m1 = in_gfs[IDX4(UUGF, i0, i1, i2 - 1)];
  const double uu_i0_i1m5_i2 = in_gfs[IDX4(UUGF, i0, i1 - 5, i2)];
  const double uu_i0_i1m4_i2 = in_gfs[IDX4(UUGF, i0, i1 - 4, i2)];
  const double uu_i0_i1m3_i2 = in_gfs[IDX4(UUGF, i0, i1 - 3, i2)];
  const double uu_i0_i1m2_i2 = in_gfs[IDX4(UUGF, i0, i1 - 2, i2)];
  const double uu_i0_i1m1_i2 = in_gfs[IDX4(UUGF, i0, i1 - 1, i2)];
  const double uu_i0m5_i1_i2 = in_gfs[IDX4(UUGF, i0 - 5, i1, i2)];
  const double uu_i0m4_i1_i2 = in_gfs[IDX4(UUGF, i0 - 4, i1, i2)];
  const double uu_i0m3_i1_i2 = in_gfs[IDX4(UUGF, i0 - 3, i1, i2)];
  const double uu_i0m2_i1_i2 = in_gfs[IDX4(UUGF, i0 - 2, i1, i2)];
  const double uu_i0m1_i1_i2 = in_gfs[IDX4(UUGF, i0 - 1, i1, i2)];
  const double uu = in_gfs[IDX4(UUGF, i0, i1, i2)];
  const double uu_i0p1_i1_i2 = in_gfs[IDX4(UUGF, i0 + 1, i1, i2)];
  const double uu_i0p2_i1_i2 = in_gfs[IDX4(UUGF, i0 + 2, i1, i2)];
  const double uu_i0p3_i1_i2 = in_gfs[IDX4(UUGF, i0 + 3, i1, i2)];
  const double uu_i0p4_i1_i2 = in_gfs[IDX4(UUGF, i0 + 4, i1, i2)];
  const double uu_i0p5_i1_i2 = in_gfs[IDX4(UUGF, i0 + 5, i1, i2)];
  const double uu_i0_i1p1_i2 = in_gfs[IDX4(UUGF, i0, i1 + 1, i2)];
  const double uu_i0_i1p2_i2 = in_gfs[IDX4(UUGF, i0, i1 + 2, i2)];
  const double uu_i0_i1p3_i2 = in_gfs[IDX4(UUGF, i0, i1 + 3, i2)];
  const double uu_i0_i1p4_i2 = in_gfs[IDX4(UUGF, i0, i1 + 4, i2)];
  const double uu_i0_i1p5_i2 = in_gfs[IDX4(UUGF, i0, i1 + 5, i2)];
  const double uu_i0_i1_i2p1 = in_gfs[IDX4(UUGF, i0, i1, i2 + 1)];
  const double uu_i0_i1_i2p2 = in_gfs[IDX4(UUGF, i0, i1, i2 + 2)];
  const double uu_i0_i1_i2p3 = in_gfs[IDX4(UUGF, i0, i1, i2 + 3)];
  const double uu_i0_i1_i2p4 = in_gfs[IDX4(UUGF, i0, i1, i2 + 4)];
  const double uu_i0_i1_i2p5 = in_gfs[IDX4(UUGF, i0, i1, i2 + 5)];
  const double vv = in_gfs[IDX4(VVGF, i0, i1, i2)];
  const double tmpFD0 = -(5269.0 / 1800.0)*uu;
  const double uu_dDD00 = pow(invdx0, 2)*(tmpFD0 + ((5.0 / 3.0))*uu_i0m1_i1_i2 - (5.0 / 21.0)*uu_i0m2_i1_i2 + ((5.0 / 126.0))*uu_i0m3_i1_i2 - (5.0 / 1008.0)*uu_i0m4_i1_i2 + ((1.0 / 3150.0))*uu_i0m5_i1_i2 + ((5.0 / 3.0))*uu_i0p1_i1_i2 - (5.0 / 21.0)*uu_i0p2_i1_i2 + ((5.0 / 126.0))*uu_i0p3_i1_i2 - (5.0 / 1008.0)*uu_i0p4_i1_i2 + ((1.0 / 3150.0))*uu_i0p5_i1_i2);
  const double uu_dDD11 = pow(invdx1, 2)*(tmpFD0 + ((5.0 / 3.0))*uu_i0_i1m1_i2 - (5.0 / 21.0)*uu_i0_i1m2_i2 + ((5.0 / 126.0))*uu_i0_i1m3_i2 - (5.0 / 1008.0)*uu_i0_i1m4_i2 + ((1.0 / 3150.0))*uu_i0_i1m5_i2 + ((5.0 / 3.0))*uu_i0_i1p1_i2 - (5.0 / 21.0)*uu_i0_i1p2_i2 + ((5.0 / 126.0))*uu_i0_i1p3_i2 - (5.0 / 1008.0)*uu_i0_i1p4_i2 + ((1.0 / 3150.0))*uu_i0_i1p5_i2);
  const double uu_dDD22 = pow(invdx2, 2)*(tmpFD0 + ((5.0 / 3.0))*uu_i0_i1_i2m1 - (5.0 / 21.0)*uu_i0_i1_i2m2 + ((5.0 / 126.0))*uu_i0_i1_i2m3 - (5.0 / 1008.0)*uu_i0_i1_i2m4 + ((1.0 / 3150.0))*uu_i0_i1_i2m5 + ((5.0 / 3.0))*uu_i0_i1_i2p1 - (5.0 / 21.0)*uu_i0_i1_i2p2 + ((5.0 / 126.0))*uu_i0_i1_i2p3 - (5.0 / 1008.0)*uu_i0_i1_i2p4 + ((1.0 / 3150.0))*uu_i0_i1_i2p5);
  /*
   * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
   */
   /*
	*  Original SymPy expressions:
	*  "[rhs_gfs[IDX4(UUGF, i0, i1, i2)] = vv,
	*    rhs_gfs[IDX4(VVGF, i0, i1, i2)] = wavespeed**2*(uu_dDD00 + uu_dDD11 + uu_dDD22)]"
	*/
  rhs_gfs[IDX4(UUGF, i0, i1, i2)] = vv;
  rhs_gfs[IDX4(VVGF, i0, i1, i2)] = pow(wavespeed, 2)*(uu_dDD00 + uu_dDD11 + uu_dDD22);
}