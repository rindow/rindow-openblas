/*
    Y(i,j) := A(is,js)

   Method Rindow\OpenBLAS\Math::
    public function slice(
        bool $reverse,
        bool $addMode,
        int $m,
        int $n,
        int $k,
        int $size,
        Buffer $A, int $offsetA, int $incA,
        Buffer $Y, int $offsetY, int $incY,
        int $startAxis0,
        int $sizeAxis0,
        int $startAxis1,
        int $sizeAxis1,
        int $startAxis2,
        int $sizeAxis2
        ) : void
 {{{ */
static PHP_METHOD(Math, slice)
{
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    php_interop_polite_math_matrix_linear_buffer_t* bufferY;
    zend_bool reverse;
    zend_bool addMode;
    zend_long m;
    zend_long n;
    zend_long k;
    zend_long size;
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
    zend_long startAxis2;
    zend_long sizeAxis2;
    zend_long i0,i1,i2;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 18, 18)
        Z_PARAM_BOOL(reverse)
        Z_PARAM_BOOL(addMode)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_LONG(k)

        Z_PARAM_LONG(size)
        Z_PARAM_OBJECT(obja) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(incA)
        Z_PARAM_OBJECT(objy) // Interop\Polite\Math\Matrix\LinearBuffer

        Z_PARAM_LONG(offsetY)
        Z_PARAM_LONG(incY)
        Z_PARAM_LONG(startAxis0)
        Z_PARAM_LONG(sizeAxis0)
        Z_PARAM_LONG(startAxis1)

        Z_PARAM_LONG(sizeAxis1)
        Z_PARAM_LONG(startAxis2)
        Z_PARAM_LONG(sizeAxis2)
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
    if(size<=0){
        zend_throw_exception(spl_ce_InvalidArgumentException, "Argument size must be greater than or equal 0.", 0);
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
    if(startAxis2<0){
        zend_throw_exception(spl_ce_InvalidArgumentException, "Argument startAxis2 must be greater than or equal 0.", 0);
        return;
    }
    if(sizeAxis2<=0){
        zend_throw_exception(spl_ce_InvalidArgumentException, "Argument sizeAxis2 must be greater than 0.", 0);
        return;
    }
    // Check Buffer A
    bufferA = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(obja);
    if(php_rindow_openblas_assert_buffer_type(bufferA,"a")) {
        return;
    }
    if(m*n*k*size*incA+offsetA > bufferA->size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch BufferA size and m,n,k,size", 0);
        return;
    }
    // Check Buffer Y
    bufferY = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(objy);
    if(php_rindow_openblas_assert_buffer_type(bufferY,"y")) {
        return;
    }
    if(sizeAxis0*sizeAxis1*size*incY+offsetY > bufferY->size) {
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
    if(startAxis2>=k||
        sizeAxis2+startAxis2>k){
        zend_throw_exception(spl_ce_InvalidArgumentException, "Axis2 range is too large for source array.",0);
        return;
    }
    if(bufferA->dtype != bufferY->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type", 0);
        return;
    }
    for(i0=0;i0<sizeAxis0;i0++) {
        for(i1=0;i1<sizeAxis1;i1++){
            for(i2=0;i2<sizeAxis2;i2++){
                zend_long pa;
                zend_long py;
                pa = (i0+startAxis0)*n*k*size*incA+
                     (i1+startAxis1)*k*size*incA+
                     (i2+startAxis2)*size*incA+
                     offsetA;
                py = i0*sizeAxis1*sizeAxis2*size+
                     i1*sizeAxis2*size*incY+
                     i2*size*incY+
                     offsetY;
                if(bufferA->dtype==php_interop_polite_math_matrix_dtype_float32){
                    float *a = &(((float *)bufferA->data)[pa]);
                    float *y = &(((float *)bufferY->data)[py]);
                    if(!reverse) {
                        if(!addMode) {
                            cblas_scopy((blasint)size,
                                a, (blasint)incA,
                                y, (blasint)incY);
                        } else {
                            cblas_saxpy((blasint)size,
                                1.0,
                                a, (blasint)incA,
                                y, (blasint)incY);
                        }
                    } else {
                        if(!addMode) {
                            cblas_scopy((blasint)size,
                                y, (blasint)incY,
                                a, (blasint)incA);
                        } else {
                            cblas_saxpy((blasint)size,
                                1.0,
                                y, (blasint)incY,
                                a, (blasint)incA);
                        }
                    }
                } else if(bufferA->dtype==php_interop_polite_math_matrix_dtype_float64){
                    double *a = &(((double *)bufferA->data)[pa]);
                    double *y = &(((double *)bufferY->data)[py]);
                    if(!reverse) {
                        if(!addMode) {
                            cblas_dcopy((blasint)size,
                                a, (blasint)incA,
                                y, (blasint)incY);
                        } else {
                            cblas_daxpy((blasint)size,
                                1.0,
                                a, (blasint)incA,
                                y, (blasint)incY);
                        }
                    } else {
                        if(!addMode) {
                            cblas_dcopy((blasint)size,
                                y, (blasint)incY,
                                a, (blasint)incA);
                        } else {
                            cblas_daxpy((blasint)size,
                                1.0,
                                y, (blasint)incY,
                                a, (blasint)incA);
                        }
                    }
                } else {
                    uint8_t *a, *y;
                    int rc;
                    int valueSize = php_rindow_openblas_common_dtype_to_valuesize(bufferA->dtype);
                    a = php_rindow_openblas_get_address(bufferA,pa,valueSize);
                    y = php_rindow_openblas_get_address(bufferY,py,valueSize);
                    if(!reverse) {
                        if(!addMode) {
                            rc = php_rindow_openblas_math_copy(
                                size, bufferA->dtype, a, incA, y, incY);
                        } else {
                            rc = php_rindow_openblas_math_add(
                                size, bufferA->dtype, a, incA, y, incY);
                        }
                    } else {
                        if(!addMode) {
                            rc = php_rindow_openblas_math_copy(
                                size, bufferA->dtype, y, incY, a, incA);
                        } else {
                            rc = php_rindow_openblas_math_add(
                                size, bufferA->dtype, y, incY, a, incA);
                        }
                    }
                    if(rc)
                        return;
                }
            }
        }
    }
}
/* }}} */
