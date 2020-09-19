/*
    Y(i,j) := A(is,js)

   Method Rindow\OpenBLAS\Math::
    public function slice(
        bool $reverse,
        bool $addMode,
        int $m,
        int $n,
        int $k,
        Buffer $A, int $offsetA, int $incA,
        Buffer $Y, int $offsetY, int $incY,
        int $startAxis0,
        int $sizeAxis0,
        int $startAxis1,
        int $sizeAxis1
        ) : void
 {{{ */
static PHP_METHOD(Math, slice)
{
    php_rindow_openblas_buffer_t* bufferA;
    php_rindow_openblas_buffer_t* bufferY;
    zend_bool reverse;
    zend_bool addMode;
    zend_long m;
    zend_long n;
    zend_long k;
    zval* obja=NULL;
    zend_long offsetA;
    zend_long incA;
    zval* objy=NULL;
    zend_long offsetY;
    zend_long incY;
    zend_long startAxis0;
    zend_long sizeAxis0;
    zend_long startAxis1;
    zend_long sizeAxis1;
    zend_long i;
    zend_long j;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 15, 15)
        Z_PARAM_BOOL(reverse)
        Z_PARAM_BOOL(addMode)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_LONG(k)

        Z_PARAM_OBJECT_OF_CLASS(obja,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(incA)
        Z_PARAM_OBJECT_OF_CLASS(objy,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetY)

        Z_PARAM_LONG(incY)
        Z_PARAM_LONG(startAxis0)
        Z_PARAM_LONG(sizeAxis0)
        Z_PARAM_LONG(startAxis1)
        Z_PARAM_LONG(sizeAxis1)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_M, m)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_K, k)) {
        return;
    }
    if(startAxis0<0){
        zend_throw_exception(spl_ce_InvalidArgumentException, "Argument startAxis0 must be greater than or equal 0.", 0);
        return;
    }
    if(sizeAxis0<=0){
        zend_throw_exception(spl_ce_InvalidArgumentException, "Argument sizeAxis0 must be greater than 0.", 0);
        return;
    }
    if(startAxis1<0){
        zend_throw_exception(spl_ce_InvalidArgumentException, "Argument startAxis1 must be greater than or equal 0.", 0);
        return;
    }
    if(sizeAxis1<=0){
        zend_throw_exception(spl_ce_InvalidArgumentException, "Argument sizeAxis1 must be greater than 0.", 0);
        return;
    }
    // Check Buffer A
    bufferA = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(obja);
    if(m*n*k*incA+offsetA > bufferA->size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch BufferA size and m,n,k", 0);
        return;
    }
    // Check Buffer Y
    bufferY = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(objy);
    if(sizeAxis0*sizeAxis1*k*incY+offsetY > bufferY->size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "BufferY size is too small", 0);
        return;
    }

    if(startAxis0>=m||
        sizeAxis0+startAxis0>m){
        zend_throw_exception(spl_ce_InvalidArgumentException, "Axis0 range is too large for source array.",0);
        return;
    }
    if(startAxis1>=n||
        sizeAxis1+startAxis1>n){
        zend_throw_exception(spl_ce_InvalidArgumentException, "Axis1 range is too large for source array.",0);
        return;
    }
    if(bufferA->dtype != bufferY->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type", 0);
        return;
    }
    for(i=0;i<sizeAxis0;i++) {
        for(j=0;j<sizeAxis1;j++){
            zend_long pa;
            zend_long py;
            pa = (i+startAxis0)*n*k*incA+(j+startAxis1)*k*incA+offsetA;
            py = i*sizeAxis1*k+j*k*incY+offsetY;
            if(bufferA->dtype==php_rindow_openblas_dtype_float32){
                float *a = &(((float *)bufferA->data)[pa]);
                float *y = &(((float *)bufferY->data)[py]);
                if(!reverse) {
                    if(!addMode) {
                        cblas_scopy((blasint)k,
                            a, (blasint)incA,
                            y, (blasint)incY);
                    } else {
                        cblas_saxpy((blasint)k,
                            1.0,
                            a, (blasint)incA,
                            y, (blasint)incY);
                    }
                } else {
                    if(!addMode) {
                        cblas_scopy((blasint)k,
                            y, (blasint)incY,
                            a, (blasint)incA);
                    } else {
                        cblas_saxpy((blasint)k,
                            1.0,
                            y, (blasint)incY,
                            a, (blasint)incA);
                    }
                }
            } else if(bufferA->dtype==php_rindow_openblas_dtype_float64){
                double *a = &(((double *)bufferA->data)[pa]);
                double *y = &(((double *)bufferY->data)[py]);
                if(!reverse) {
                    cblas_dcopy((blasint)k,
                        a, (blasint)incA,
                        y, (blasint)incY);
                } else {
                    cblas_dcopy((blasint)k,
                        y, (blasint)incY,
                        a, (blasint)incA);
                }
            } else {
                zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
                return;
            }
        }
    }
}
/* }}} */
