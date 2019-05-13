/*
 *  Original SymPy expressions:
 *  "[parity[0] = 1,
 *    parity[1] = sin(xx1)*sin(xx1_inbounds)*sin(xx2)*sin(xx2_inbounds) + sin(xx1)*sin(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds) + cos(xx1)*cos(xx1_inbounds),
 *    parity[2] = sin(xx1)*sin(xx1_inbounds) + sin(xx2)*sin(xx2_inbounds)*cos(xx1)*cos(xx1_inbounds) + cos(xx1)*cos(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds),
 *    parity[3] = sin(xx2)*sin(xx2_inbounds) + cos(xx2)*cos(xx2_inbounds),
 *    parity[4] = (sin(xx1)*sin(xx1_inbounds)*sin(xx2)*sin(xx2_inbounds) + sin(xx1)*sin(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds) + cos(xx1)*cos(xx1_inbounds))**2,
 *    parity[5] = (sin(xx1)*sin(xx1_inbounds) + sin(xx2)*sin(xx2_inbounds)*cos(xx1)*cos(xx1_inbounds) + cos(xx1)*cos(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds))*(sin(xx1)*sin(xx1_inbounds)*sin(xx2)*sin(xx2_inbounds) + sin(xx1)*sin(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds) + cos(xx1)*cos(xx1_inbounds)),
 *    parity[6] = (sin(xx2)*sin(xx2_inbounds) + cos(xx2)*cos(xx2_inbounds))*(sin(xx1)*sin(xx1_inbounds)*sin(xx2)*sin(xx2_inbounds) + sin(xx1)*sin(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds) + cos(xx1)*cos(xx1_inbounds)),
 *    parity[7] = (sin(xx1)*sin(xx1_inbounds) + sin(xx2)*sin(xx2_inbounds)*cos(xx1)*cos(xx1_inbounds) + cos(xx1)*cos(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds))**2,
 *    parity[8] = (sin(xx2)*sin(xx2_inbounds) + cos(xx2)*cos(xx2_inbounds))*(sin(xx1)*sin(xx1_inbounds) + sin(xx2)*sin(xx2_inbounds)*cos(xx1)*cos(xx1_inbounds) + cos(xx1)*cos(xx1_inbounds)*cos(xx2)*cos(xx2_inbounds)),
 *    parity[9] = (sin(xx2)*sin(xx2_inbounds) + cos(xx2)*cos(xx2_inbounds))**2]"
 */
{
   const double tmp0 = cos(xx1)*cos(xx1_inbounds);
   const double tmp1 = sin(xx1)*sin(xx1_inbounds);
   const double tmp2 = sin(xx2)*sin(xx2_inbounds);
   const double tmp3 = cos(xx2)*cos(xx2_inbounds);
   const double tmp4 = tmp0 + tmp1*tmp2 + tmp1*tmp3;
   const double tmp5 = tmp0*tmp2 + tmp0*tmp3 + tmp1;
   const double tmp6 = tmp2 + tmp3;
   parity[0] = 1;
   parity[1] = tmp4;
   parity[2] = tmp5;
   parity[3] = tmp6;
   parity[4] = pow(tmp4, 2);
   parity[5] = tmp4*tmp5;
   parity[6] = tmp4*tmp6;
   parity[7] = pow(tmp5, 2);
   parity[8] = tmp5*tmp6;
   parity[9] = pow(tmp6, 2);
}
