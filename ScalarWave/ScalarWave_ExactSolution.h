#pragma omp parallel for
for(int i2=0; i2<Nxx_plus_2NGHOSTS[2]; i2++) {
    const REAL xx2=xx[2][i2];
    for(int i1=0; i1<Nxx_plus_2NGHOSTS[1]; i1++) {
        const REAL xx1=xx[1][i1];
        for(int i0=0; i0<Nxx_plus_2NGHOSTS[0]; i0++) {
            const REAL xx0=xx[0][i0];
            {
               /* 
                * NRPy+ Finite Difference Code Generation, Step 1 of 1: Evaluate SymPy expressions and write to main memory:
                */
               /*
                *  Original SymPy expressions:
                *  "[in_gfs[IDX4(UUGF, i0, i1, i2)] = -sin(time*wavespeed - (kk0*xx0 + kk1*xx1 + kk2*xx2)/sqrt(kk0**2 + kk1**2 + kk2**2)) + 2,
                *    in_gfs[IDX4(VVGF, i0, i1, i2)] = -wavespeed*cos(time*wavespeed - (kk0*xx0 + kk1*xx1 + kk2*xx2)/sqrt(kk0**2 + kk1**2 + kk2**2))]"
                */
               const double tmp0 = time*wavespeed - (kk0*xx0 + kk1*xx1 + kk2*xx2)/sqrt(pow(kk0, 2) + pow(kk1, 2) + pow(kk2, 2));
               in_gfs[IDX4(UUGF, i0, i1, i2)] = -sin(tmp0) + 2;
               in_gfs[IDX4(VVGF, i0, i1, i2)] = -wavespeed*cos(tmp0);
            }
            
            
        } // END LOOP: for(int i0=0; i0<Nxx_plus_2NGHOSTS[0]; i0++)
    } // END LOOP: for(int i1=0; i1<Nxx_plus_2NGHOSTS[1]; i1++)
} // END LOOP: for(int i2=0; i2<Nxx_plus_2NGHOSTS[2]; i2++)
