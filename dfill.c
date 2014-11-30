
/*
 *  function dfill : intialises double precision array to scalar value
 */
  void
  dfill(const int n, const double val, double a[], const int ia){
    int i;

    for (i=0; i<(n-1)*ia+1; i+=ia)
      a[i] = val;
  }
