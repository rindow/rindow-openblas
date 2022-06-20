/*
   B(n,k) := A(X(n),k)

   Method Rindow\OpenBLAS\Math::
   public function gather(
       bool $reverse,
       bool $addMode,
       int $n,
       int $k,
       int $numClass,
       Buffer $X, int $offsetX,
       Buffer $A, int $offsetA,
       Buffer $B, int $offsetB
       ) : void
 {{{ */
static PHP_METHOD(Math, gather)
{
    zend_bool reverse;
    zend_bool addMode;
    zend_long n;
    zend_long k;
    zend_long numClass;
    zval* x;
    zend_long offsetX;
    zval* a;
    zend_long offsetA;
    zval* b;
    zend_long offsetB;
    php_interop_polite_math_matrix_linear_buffer_t* bufferX;
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    php_interop_polite_math_matrix_linear_buffer_t* bufferB;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 11, 11)
        Z_PARAM_BOOL(reverse)
        Z_PARAM_BOOL(addMode)
        Z_PARAM_LONG(n)
        Z_PARAM_LONG(k)
        Z_PARAM_LONG(numClass)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_OBJECT(a) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetA)
        Z_PARAM_OBJECT(b) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetB)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        "k", k)) {
        return;
    }
    if(numClass<=0) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Argument numClass must be greater than or equal 0.", 0);
        return;
    }

    // Check Buffer X
    bufferX = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(bufferX,"x")) {
        return;
    }
    if(offsetX<0) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Argument offsetX must be greater than or equal 0.", 0);
        return;
    }
    if(offsetX+n > bufferX->size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Matrix X specification too large for buffer.", 0);
        return;
    }

    // Check Buffer A
    bufferA = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_buffer_type(bufferA,"a")) {
        return;
    }
    if(offsetA<0) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Argument offsetA must be greater than or equal 0.", 0);
        return;
    }
    if(offsetA+numClass*k > bufferA->size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Matrix A specification too large for buffer.", 0);
        return;
    }

    // Check Buffer B
    bufferB = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(b);
    if(php_rindow_openblas_assert_buffer_type(bufferB,"b")) {
        return;
    }
    if(offsetB<0) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Argument offsetB must be greater than or equal 0.", 0);
        return;
    }
    if(offsetB+n*k > bufferB->size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Matrix B specification too large for buffer.", 0);
        return;
    }

    // Check Buffer A and Y
    if(bufferA->dtype!=bufferB->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for A and B", 0);
        return;
    }
    if(bufferX->dtype==php_interop_polite_math_matrix_dtype_bool) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Data type of BufferX must not be bool", 0);
        return;
    }

    switch (bufferA->dtype) {
        case php_interop_polite_math_matrix_dtype_float32: {
            float *a = &(((float *)bufferA->data)[offsetA]);
            float *b = &(((float *)bufferB->data)[offsetB]);
            zend_long selector;
            size_t ldB = k;
            size_t ldIndex = k;
            for(zend_long j=0; j<n; j++,b+=ldB) {
                if(rindow_openblas_math_get_integer(
                            bufferX->dtype, bufferX->data, offsetX,1,
                            j, &selector)) {
                    zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type of label number.", 0);
                    return;
                }
                if(selector<0||selector>=numClass) {
                    zend_throw_exception(spl_ce_RuntimeException, "Label number is out of bounds.", 0);
                    return;
                }
                if(reverse) {
                    if(addMode) {
                        if(k==1) {
                            a[selector*ldIndex] += *b;
                        } else {
                            cblas_saxpy((blasint)k, 1.0,
                                b, (blasint)1,
                                &(a[selector*ldIndex]), (blasint)1);
                        }
                    } else {
                        if(k==1) {
                            a[selector*ldIndex] = *b;
                        } else {
                            cblas_scopy((blasint)k,
                                b, (blasint)1,
                                &(a[selector*ldIndex]), (blasint)1);
                        }
                    }
                } else {
                    if(addMode) {
                        if(k==1) {
                            *b += a[selector*ldIndex];
                        } else {
                            cblas_saxpy((blasint)k, 1.0,
                                &(a[selector*ldIndex]), (blasint)1,
                                b, (blasint)1);
                        }
                    } else {
                        if(k==1) {
                            *b = a[selector*ldIndex];
                        } else {
                            cblas_scopy((blasint)k,
                                &(a[selector*ldIndex]), (blasint)1,
                                b, (blasint)1);
                        }
                    }
                }
            }
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            double *a = &(((double *)bufferA->data)[offsetA]);
            double *b = &(((double *)bufferB->data)[offsetB]);
            zend_long selector;
            size_t ldB = k;
            size_t ldIndex = k;
            for(zend_long j=0; j<n; j++,b+=ldB) {
                if(rindow_openblas_math_get_integer(
                            bufferX->dtype, bufferX->data, offsetX,1,
                            j, &selector)) {
                    zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type of label number.", 0);
                    return;
                }
                if(selector<0||selector>=numClass) {
                    zend_throw_exception(spl_ce_RuntimeException, "Label number is out of bounds.", 0);
                    return;
                }
                if(reverse) {
                    if(addMode) {
                        if(k==1) {
                            a[selector*ldIndex] += *b;
                        } else {
                            cblas_daxpy((blasint)k, 1.0,
                                b, (blasint)1,
                                &(a[selector*ldIndex]), (blasint)1);
                        }
                    } else {
                        if(k==1) {
                            a[selector*ldIndex] = *b;
                        } else {
                            cblas_dcopy((blasint)k,
                                b, (blasint)1,
                                &(a[selector*ldIndex]), (blasint)1);
                        }
                    }
                } else {
                    if(addMode) {
                        if(k==1) {
                            *b += a[selector*ldIndex];
                        } else {
                            cblas_daxpy((blasint)k, 1.0,
                                &(a[selector*ldIndex]), (blasint)1,
                                b, (blasint)1);
                        }
                    } else {
                        if(k==1) {
                            *b = a[selector*ldIndex];
                        } else {
                            cblas_dcopy((blasint)k,
                                &(a[selector*ldIndex]), (blasint)1,
                                b, (blasint)1);
                        }
                    }
                }
            }
            break;
        }
        default: {
            if(!php_rindow_openblas_common_dtype_is_int(bufferA->dtype)&&
                !php_rindow_openblas_common_dtype_is_bool(bufferA->dtype)) {
                zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
                return;
            }

            int valueSize = php_rindow_openblas_common_dtype_to_valuesize(bufferA->dtype);
            zend_long selector;
            size_t ldB = k*valueSize;
            size_t ldIndex = k*valueSize;
            uint8_t *a, *b;
            int rc;
            a = php_rindow_openblas_get_address(bufferA,offsetA,valueSize);
            b = php_rindow_openblas_get_address(bufferB,offsetB,valueSize);

            for(zend_long j=0; j<n; j++,b+=ldB) {
                if(rindow_openblas_math_get_integer(
                            bufferX->dtype, bufferX->data, offsetX,1,
                            j, &selector)) {
                    zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type of label number.", 0);
                    return;
                }
                if(selector<0||selector>=numClass) {
                    zend_throw_exception(spl_ce_RuntimeException, "Label number is out of bounds.", 0);
                    return;
                }
                if(reverse) {
                    if(addMode) {
                        rc = php_rindow_openblas_math_add(
                            k, bufferA->dtype, b, 1, a+selector*ldIndex, 1);
                    } else {
                        rc = php_rindow_openblas_math_copy(
                            k, bufferA->dtype, b, 1, a+selector*ldIndex, 1);
                    }
                } else {
                    if(addMode) {
                        rc = php_rindow_openblas_math_add(
                            k, bufferA->dtype, a+selector*ldIndex, 1, b, 1);
                    } else {
                        rc = php_rindow_openblas_math_copy(
                            k, bufferA->dtype, a+selector*ldIndex, 1, b, 1);
                    }
                }
                if(rc)
                    return;
            }
            break;
        }
    }
}

/*
   B(m,n) := A(m,X(m,n))

   public function reduceGather(
       bool $reverse,
       bool $addMode,
       int $m,
       int $n,
       int $numClass,
       Buffer $X, int $offsetX,
       Buffer $A, int $offsetA,
       Buffer $B, int $offsetB
       ) : void
 {{{ */
static PHP_METHOD(Math, reduceGather)
{
    zend_bool reverse;
    zend_bool addMode;
    zend_long m;
    zend_long n;
    zend_long numClass;
    zval* x;
    zend_long offsetX;
    zval* a;
    zend_long offsetA;
    zval* b;
    zend_long offsetB;
    php_interop_polite_math_matrix_linear_buffer_t* bufferX;
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    php_interop_polite_math_matrix_linear_buffer_t* bufferB;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 11, 11)
        Z_PARAM_BOOL(reverse)
        Z_PARAM_BOOL(addMode)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_LONG(numClass)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_OBJECT(a) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetA)
        Z_PARAM_OBJECT(b) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetB)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "m", m)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    if(numClass<=0) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Argument numClass must be greater than or equal 0.", 0);
        return;
    }
    // Check Buffer X
    bufferX = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(bufferX,"x")) {
        return;
    }
    if(offsetX<0) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Argument offsetX must be greater than or equal 0.", 0);
        return;
    }
    if(offsetX+m*n > bufferX->size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Matrix X specification too large for buffer.", 0);
        return;
    }

    // Check Buffer A
    bufferA = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_buffer_type(bufferA,"a")) {
        return;
    }
    if(offsetA<0) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Argument offsetA must be greater than or equal 0.", 0);
        return;
    }
    if(offsetA+m*numClass > bufferA->size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Matrix A specification too large for buffer.", 0);
        return;
    }

    // Check Buffer B
    bufferB = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(b);
    if(php_rindow_openblas_assert_buffer_type(bufferB,"b")) {
        return;
    }
    if(offsetB<0) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Argument offsetB must be greater than or equal 0.", 0);
        return;
    }
    if(offsetB+m*n > bufferB->size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Matrix B specification too large for buffer.", 0);
        return;
    }

    // Check Buffer A and Y
    if(bufferA->dtype!=bufferB->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for A and B", 0);
        return;
    }
    if(bufferX->dtype==php_interop_polite_math_matrix_dtype_bool) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Data type of BufferX must not be bool", 0);
        return;
    }

    switch (bufferA->dtype) {
        case php_interop_polite_math_matrix_dtype_float32: {
            float *a = &(((float *)bufferA->data)[offsetA]);
            float *b = &(((float *)bufferB->data)[offsetB]);
            zend_long selector;
            zend_long idxX = offsetX;
            zend_long ldX = n;
            zend_long ldA = n*numClass;
            zend_long ldB = n;
            for(zend_long i=0; i<m; i++,idxX+=ldX,a+=ldA,b+=ldB) {
                for(zend_long j=0; j<n; j++) {
                    if(rindow_openblas_math_get_integer(
                                bufferX->dtype, bufferX->data, idxX, 1,
                                j, &selector)) {
                        zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type of label number.", 0);
                        return;
                    }
                    if(selector<0||selector>=numClass) {
                        zend_throw_exception(spl_ce_RuntimeException, "Label number is out of bounds.", 0);
                        return;
                    }
                    if(reverse) {
                        if(addMode) {
                            a[j+selector*ldB] += b[j];
                        } else {
                            a[j+selector*ldB] = b[j];
                        }
                    } else {
                        if(addMode) {
                            b[j] += a[j+selector*ldB];
                        } else {
                            b[j] = a[j+selector*ldB];
                        }
                    }
                }
            }
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            double *a = &(((double *)bufferA->data)[offsetA]);
            double *b = &(((double *)bufferB->data)[offsetB]);
            zend_long selector;
            zend_long idxX = offsetX;
            zend_long ldX = n;
            zend_long ldA = n*numClass;
            zend_long ldB = n;
            for(zend_long i=0; i<m; i++,idxX+=ldX,a+=ldA,b+=ldB) {
                for(zend_long j=0; j<n; j++) {
                    if(rindow_openblas_math_get_integer(
                                bufferX->dtype, bufferX->data, idxX, 1,
                                j, &selector)) {
                        zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type of label number.", 0);
                        return;
                    }
                    if(selector<0||selector>=numClass) {
                        zend_throw_exception(spl_ce_RuntimeException, "Label number is out of bounds.", 0);
                        return;
                    }
                    if(reverse) {
                        if(addMode) {
                            a[j+selector*ldB] += b[j];
                        } else {
                            a[j+selector*ldB] = b[j];
                        }
                    } else {
                        if(addMode) {
                            b[j] += a[j+selector*ldB];
                        } else {
                            b[j] = a[j+selector*ldB];
                        }
                    }
                }
            }
            break;
        }
        default: {
            if(!php_rindow_openblas_common_dtype_is_int(bufferA->dtype)&&
                !php_rindow_openblas_common_dtype_is_bool(bufferA->dtype)) {
                zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
                return;
            }

            int valueSize = php_rindow_openblas_common_dtype_to_valuesize(bufferA->dtype);
            zend_long selector;
            zend_long idxX = offsetX;
            zend_long ldX = n;
            zend_long ldA = n*numClass*valueSize;
            zend_long ldB = n*valueSize;
            zend_long ldIndex = n*valueSize;
            uint8_t *a, *b;
            int rc;
            a = php_rindow_openblas_get_address(bufferA,offsetA,valueSize);
            b = php_rindow_openblas_get_address(bufferB,offsetB,valueSize);

            for(zend_long i=0; i<m; i++,idxX+=ldX,a+=ldA,b+=ldB) {
                for(zend_long j=0; j<n; j++) {
                    if(rindow_openblas_math_get_integer(
                                bufferX->dtype, bufferX->data, idxX,1,
                                j, &selector)) {
                        zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type of label number.", 0);
                        return;
                    }
                    if(selector<0||selector>=numClass) {
                        zend_throw_exception(spl_ce_RuntimeException, "Label number is out of bounds.", 0);
                        return;
                    }
                    if(reverse) {
                        if(addMode) {
                            rc = php_rindow_openblas_math_add(
                                1, bufferA->dtype,
                                b+j*valueSize, 1,
                                a+j*valueSize+selector*ldIndex, 1);
                        } else {
                            rc = php_rindow_openblas_math_copy(
                                1, bufferA->dtype,
                                b+j*valueSize, 1,
                                a+j*valueSize+selector*ldIndex, 1);
                        }
                    } else {
                        if(addMode) {
                            rc = php_rindow_openblas_math_add(
                                1, bufferA->dtype,
                                a+j*valueSize+selector*ldIndex, 1,
                                b+j*valueSize, 1);
                        } else {
                            rc = php_rindow_openblas_math_copy(
                                1, bufferA->dtype,
                                a+j*valueSize+selector*ldIndex, 1,
                                b+j*valueSize, 1);
                        }
                    }
                    if(rc)
                        return;
                }
            }
            break;
        }
    }
}
/* }}} */
