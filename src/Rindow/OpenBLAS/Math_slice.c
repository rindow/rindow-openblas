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
        "m", m)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        "k", k)) {
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

    switch (bufferA->dtype) {
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataY,bufferY,offsetY)
            rindow_matlib_s_slice(reverse,addMode,m,n,k,size,pDataA,incA,pDataY,incY,startAxis0,sizeAxis0,startAxis1,sizeAxis1,startAxis2,sizeAxis2);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataY,bufferY,offsetY)
            rindow_matlib_d_slice(reverse,addMode,m,n,k,size,pDataA,incA,pDataY,incY,startAxis0,sizeAxis0,startAxis1,sizeAxis1,startAxis2,sizeAxis2);
            break;
        }
        default:{
            if(!php_rindow_openblas_common_dtype_is_int(bufferA->dtype)&&
                !php_rindow_openblas_common_dtype_is_bool(bufferA->dtype)) {
                zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
                return;
            }
            void *pDataA = rindow_matlib_common_get_address(bufferA->dtype, bufferA->data,offsetA);
            void *pDataY = rindow_matlib_common_get_address(bufferY->dtype, bufferY->data,offsetY);
            rindow_matlib_i_slice(reverse,addMode,m,n,k,size,bufferA->dtype,pDataA,incA,pDataY,incY,startAxis0,sizeAxis0,startAxis1,sizeAxis1,startAxis2,sizeAxis2);
            break;
        }
    }
}
/* }}} */
