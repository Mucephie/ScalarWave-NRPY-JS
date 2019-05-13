const REAL invdx0 = 1.0/dxx[0];
const REAL invdx1 = 1.0/dxx[1];
const REAL invdx2 = 1.0/dxx[2];
#pragma omp parallel for
for(int i2=NGHOSTS; i2<NGHOSTS+Nxx[2]; i2++) {
    for(int i1=NGHOSTS; i1<NGHOSTS+Nxx[1]; i1++) {
        for(int i0=NGHOSTS; i0<NGHOSTS+Nxx[0]; i0++) {
            {
               /* 
                * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
                */
               /*
                *  Original SymPy expressions:
                *  "[const double uu_dDD00 = invdx0**2*(-5*uu/2 + 4*uu_i0m1_i1_i2/3 - uu_i0m2_i1_i2/12 + 4*uu_i0p1_i1_i2/3 - uu_i0p2_i1_i2/12),
                *    const double uu_dDD11 = invdx1**2*(-5*uu/2 + 4*uu_i0_i1m1_i2/3 - uu_i0_i1m2_i2/12 + 4*uu_i0_i1p1_i2/3 - uu_i0_i1p2_i2/12),
                *    const double uu_dDD22 = invdx2**2*(-5*uu/2 + 4*uu_i0_i1_i2m1/3 - uu_i0_i1_i2m2/12 + 4*uu_i0_i1_i2p1/3 - uu_i0_i1_i2p2/12)]"
                */
               const double uu_i0_i1_i2m2 = in_gfs[IDX4(UUGF, i0,i1,i2-2)];
               const double uu_i0_i1_i2m1 = in_gfs[IDX4(UUGF, i0,i1,i2-1)];
               const double uu_i0_i1m2_i2 = in_gfs[IDX4(UUGF, i0,i1-2,i2)];
               const double uu_i0_i1m1_i2 = in_gfs[IDX4(UUGF, i0,i1-1,i2)];
               const double uu_i0m2_i1_i2 = in_gfs[IDX4(UUGF, i0-2,i1,i2)];
               const double uu_i0m1_i1_i2 = in_gfs[IDX4(UUGF, i0-1,i1,i2)];
               const double uu = in_gfs[IDX4(UUGF, i0,i1,i2)];
               const double uu_i0p1_i1_i2 = in_gfs[IDX4(UUGF, i0+1,i1,i2)];
               const double uu_i0p2_i1_i2 = in_gfs[IDX4(UUGF, i0+2,i1,i2)];
               const double uu_i0_i1p1_i2 = in_gfs[IDX4(UUGF, i0,i1+1,i2)];
               const double uu_i0_i1p2_i2 = in_gfs[IDX4(UUGF, i0,i1+2,i2)];
               const double uu_i0_i1_i2p1 = in_gfs[IDX4(UUGF, i0,i1,i2+1)];
               const double uu_i0_i1_i2p2 = in_gfs[IDX4(UUGF, i0,i1,i2+2)];
               const double vv = in_gfs[IDX4(VVGF, i0,i1,i2)];
               const double tmpFD0 = -(5.0 / 2.0)*uu;
               const double uu_dDD00 = pow(invdx0, 2)*(tmpFD0 + ((4.0 / 3.0))*uu_i0m1_i1_i2 - (1.0 / 12.0)*uu_i0m2_i1_i2 + ((4.0 / 3.0))*uu_i0p1_i1_i2 - (1.0 / 12.0)*uu_i0p2_i1_i2);
               const double uu_dDD11 = pow(invdx1, 2)*(tmpFD0 + ((4.0 / 3.0))*uu_i0_i1m1_i2 - (1.0 / 12.0)*uu_i0_i1m2_i2 + ((4.0 / 3.0))*uu_i0_i1p1_i2 - (1.0 / 12.0)*uu_i0_i1p2_i2);
               const double uu_dDD22 = pow(invdx2, 2)*(tmpFD0 + ((4.0 / 3.0))*uu_i0_i1_i2m1 - (1.0 / 12.0)*uu_i0_i1_i2m2 + ((4.0 / 3.0))*uu_i0_i1_i2p1 - (1.0 / 12.0)*uu_i0_i1_i2p2);
               /* 
                * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
                */
               /*
                *  Original SymPy expressions:
                *  "[rhs_gfs[IDX4(UUGF, i0, i1, i2)] = vv,
                *    rhs_gfs[IDX4(VVGF, i0, i1, i2)] = uu_dDD00*wavespeed**2 + uu_dDD11*wavespeed**2 + uu_dDD22*wavespeed**2]"
                */
               const double tmp0 = pow(wavespeed, 2);
               rhs_gfs[IDX4(UUGF, i0, i1, i2)] = vv;
               rhs_gfs[IDX4(VVGF, i0, i1, i2)] = tmp0*uu_dDD00 + tmp0*uu_dDD11 + tmp0*uu_dDD22;
            }
            
            
        } // END LOOP: for(int i0=NGHOSTS; i0<NGHOSTS+Nxx[0]; i0++)
    } // END LOOP: for(int i1=NGHOSTS; i1<NGHOSTS+Nxx[1]; i1++)
} // END LOOP: for(int i2=NGHOSTS; i2<NGHOSTS+Nxx[2]; i2++)
