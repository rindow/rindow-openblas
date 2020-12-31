/*
   Y(i) := A(X[i])

   Method Rindow\OpenBLAS\Math::
    public function selectAxis0(
        int $m,
        int $n,
        int $k,
        Buffer $A, int $offsetA, int $ldA,
        Buffer $X, int $offsetX, int $incX,
        Buffer $Y, int $offsetY, int $ldY ) : void
 {{{ */
static PHP_METHOD(Math, selectAxis0)
{
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    php_interop_polite_math_matrix_linear_buffer_t* bufferX;
    php_interop_polite_math_matrix_linear_buffer_t* bufferY;
    zend_long m;
    zend_long n;
    zend_long k;
    zval* a=NULL;
    zend_long offsetA;
    zend_long ldA;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zval* y=NULL;
    zend_long offsetY;
    zend_long ldY;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 12, 12)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_LONG(k)
        Z_PARAM_OBJECT(a) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_OBJECT(y) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetY)
        Z_PARAM_LONG(ldY)
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
    // Check Buffer A
    bufferA = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_buffer_type(bufferA,"a")) {
        return;
    }
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_A, bufferA,m,n,offsetA,ldA)) {
        return;
    }

    // Check Buffer X
    bufferX = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(bufferX,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, bufferX,k,offsetX,incX)) {
        return;
    }

    // Check Buffer Y
    bufferY = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(y);
    if(php_rindow_openblas_assert_buffer_type(bufferY,"y")) {
        return;
    }
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_Y, bufferY,k,n,offsetY,ldY)) {
        return;
    }

    // Check Buffer A and Y
    if(bufferA->dtype!=bufferY->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for A and Y", 0);
        return;
    }
    if(bufferX->dtype==php_interop_polite_math_matrix_dtype_bool) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Data type of BufferX must not be bool", 0);
        return;
    }

    switch (bufferA->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *a = &(((float *)bufferA->data)[offsetA]);
                float *y = &(((float *)bufferY->data)[offsetY]);
                zend_long i,selector;
                for(i=0; i<k; i++,y+=ldY) {
                    if(rindow_openblas_math_get_integer(
                                bufferX->dtype, bufferX->data, offsetX,incX,
                                i, &selector)) {
                        zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type of label number.", 0);
                        return;
                    }
                    if(selector<0||selector>=m) {
                        zend_throw_exception(spl_ce_RuntimeException, "Label number is out of bounds.", 0);
                        return;
                    }
                    if(n==1) {
                        *y = a[selector*ldA];
                    } else {
                        cblas_scopy((blasint)n,
                            &(a[selector*ldA]), (blasint)1,
                            y, (blasint)1);
                    }
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *a = &(((double *)bufferA->data)[offsetA]);
                double *y = &(((double *)bufferY->data)[offsetX]);
                zend_long i,selector;
                for(i=0; i<k; i++,y+=ldY) {
                    if(rindow_openblas_math_get_integer(
                                bufferX->dtype, bufferX->data, offsetX,incX,
                                i, &selector)) {
                        zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type of label number.", 0);
                        return;
                    }
                    if(selector<0||selector>=m) {
                        zend_throw_exception(spl_ce_RuntimeException, "Label number is out of bounds.", 0);
                        return;
                    }
                    if(n==1) {
                        *y = a[selector*ldA];
                    } else {
                        cblas_dcopy((blasint)n,
                            &(a[selector*ldA]), (blasint)1,
                            y, (blasint)1);
                    }
                }
            }
            break;
        default:
            if(!php_rindow_openblas_common_dtype_is_int(bufferA->dtype)&&
                !php_rindow_openblas_common_dtype_is_bool(bufferA->dtype)) {
                zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
                return;
            }
            {
                int valueSize;
                uint8_t *a, *y;
                valueSize = php_rindow_openblas_common_dtype_to_valuesize(bufferA->dtype);

                a = php_rindow_openblas_get_address(bufferA,offsetA,valueSize);
                y = php_rindow_openblas_get_address(bufferY,offsetY,valueSize);
                zend_long i,selector,incrementY,selectSizeA,copySize;
                incrementY = ldY*valueSize;
                selectSizeA = ldA*valueSize;
                copySize = n*valueSize;
                for(i=0; i<k; i++,y+=incrementY) {
                    if(rindow_openblas_math_get_integer(
                                bufferX->dtype, bufferX->data, offsetX,incX,
                                i, &selector)) {
                        zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type of label number.", 0);
                        return;
                    }
                    if(selector<0||selector>=m) {
                        zend_throw_exception(spl_ce_RuntimeException, "Label number is out of bounds.", 0);
                        return;
                    }
                    memcpy(y, &(a[selector*selectSizeA]), copySize);
                }
            }
            break;
    }
}

/*
   Y(i) := A(X[i])

   Method Rindow\OpenBLAS\Math::
    public function selectAxis1(
        int $m,
        int $n,
        Buffer $A, int $offsetA, int $ldA,
        Buffer $X, int $offsetX, int $incX,
        Buffer $Y, int $offsetY, int $incY ) : void
 {{{ */
static PHP_METHOD(Math, selectAxis1)
{
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    php_interop_polite_math_matrix_linear_buffer_t* bufferX;
    php_interop_polite_math_matrix_linear_buffer_t* bufferY;
    zend_long m;
    zend_long n;
    zval* a=NULL;
    zend_long offsetA;
    zend_long ldA;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zval* y=NULL;
    zend_long offsetY;
    zend_long incY;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 11, 11)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(a) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_OBJECT(y) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetY)
        Z_PARAM_LONG(incY)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_M, m)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    // Check Buffer A
    bufferA = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_buffer_type(bufferA,"a")) {
        return;
    }
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_A, bufferA,m,n,offsetA,ldA)) {
        return;
    }

    // Check Buffer X
    bufferX = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(bufferX,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, bufferX,m,offsetX,incX)) {
        return;
    }

    // Check Buffer Y
    bufferY = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(y);
    if(php_rindow_openblas_assert_buffer_type(bufferY,"y")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_Y, bufferY,m,offsetY,incY)) {
        return;
    }

    // Check Buffer A and Y
    if(bufferA->dtype!=bufferY->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for A and Y", 0);
        return;
    }
    if(bufferX->dtype==php_interop_polite_math_matrix_dtype_bool) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Data type of BufferX must not be bool", 0);
        return;
    }

    switch (bufferA->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *a = &(((float *)bufferA->data)[offsetA]);
                float *y = &(((float *)bufferY->data)[offsetY]);
                zend_long i,selector;
                for(i=0; i<m; i++) {
                    if(rindow_openblas_math_get_integer(
                                bufferX->dtype, bufferX->data, offsetX,incX,
                                i, &selector)) {
                        zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type of label number.", 0);
                        return;
                    }
                    if(selector<0||selector>=n) {
                        zend_throw_exception(spl_ce_RuntimeException, "Label number is out of bounds.", 0);
                        return;
                    }
                    y[i*incY] = a[i*ldA+selector];
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *a = &(((double *)bufferA->data)[offsetA]);
                double *y = &(((double *)bufferY->data)[offsetX]);
                zend_long i,selector;
                for(i=0; i<m; i++) {
                    if(rindow_openblas_math_get_integer(
                                bufferX->dtype, bufferX->data, offsetX,incX,
                                i, &selector)) {
                        zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type of label number.", 0);
                        return;
                    }
                    if(selector<0||selector>=n) {
                        zend_throw_exception(spl_ce_RuntimeException, "Label number is out of bounds.", 0);
                        return;
                    }
                    y[i*incY] = a[i*ldA+selector];
                }
            }
            break;
        default:
            zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
            return;
    }
}
/* }}} */
