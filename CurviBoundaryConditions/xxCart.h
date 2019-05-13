/*
 *  Original SymPy expressions:
 *  "[xCart[0] = xx0*sin(xx1)*cos(xx2),
 *    xCart[1] = xx0*sin(xx1)*sin(xx2),
 *    xCart[2] = xx0*cos(xx1)]"
 */
{
   const double tmp0 = xx0*sin(xx1);
   xCart[0] = tmp0*cos(xx2);
   xCart[1] = tmp0*sin(xx2);
   xCart[2] = xx0*cos(xx1);
}
