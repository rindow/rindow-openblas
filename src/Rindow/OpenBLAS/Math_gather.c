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
            void *pDataX = rindow_matlib_common_get_address((dtype_t)bufferX->dtype, bufferX->data,(index_t)offsetX);
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataB,bufferB,offsetB)
            if(pDataX==NULL) {
                zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type of label number.", 0);
                return;
            }
            int32_t errcode = rindow_matlib_s_gather(reverse,addMode,(index_t)n,(index_t)k,(index_t)numClass,(dtype_t)bufferX->dtype,pDataX,pDataA,pDataB);
            if(errcode) {
                if(errcode == RINDOW_MATLIB_E_UNSUPPORTED_DATA_TYPE) {
                    zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type of label number.", 0);
                    return;
                } else if(errcode == RINDOW_MATLIB_E_PERM_OUT_OF_RANGE) {
                    zend_throw_exception(spl_ce_RuntimeException, "Label number is out of bounds.", 0);
                    return;
                } else {
                    zend_throw_exception_ex(spl_ce_InvalidArgumentException, 0, "Unknown error.(%d)", errcode);
                    return;
                }
            }
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            void *pDataX = rindow_matlib_common_get_address((dtype_t)bufferX->dtype, bufferX->data,(index_t)offsetX);
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataB,bufferB,offsetB)
            if(pDataX==NULL) {
                zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type of label number.", 0);
                return;
            }
            int32_t errcode = rindow_matlib_d_gather(reverse,addMode,(index_t)n,(index_t)k,(index_t)numClass,(dtype_t)bufferX->dtype,pDataX,pDataA,pDataB);
            if(errcode) {
                if(errcode == RINDOW_MATLIB_E_UNSUPPORTED_DATA_TYPE) {
                    zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type of label number.", 0);
                    return;
                } else if(errcode == RINDOW_MATLIB_E_PERM_OUT_OF_RANGE) {
                    zend_throw_exception(spl_ce_RuntimeException, "Label number is out of bounds.", 0);
                    return;
                } else {
                    zend_throw_exception_ex(spl_ce_InvalidArgumentException, 0, "Unknown error.: %d", errcode);
                    return;
                }
            }
            break;
        }
        default: {
            if(!php_rindow_openblas_common_dtype_is_int(bufferA->dtype)&&
                !php_rindow_openblas_common_dtype_is_bool(bufferA->dtype)) {
                zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
                return;
            }
            void *pDataX = rindow_matlib_common_get_address((dtype_t)bufferX->dtype, bufferX->data,(index_t)offsetX);
            void *pDataA = rindow_matlib_common_get_address((dtype_t)bufferA->dtype, bufferA->data,(index_t)offsetA);
            void *pDataB = rindow_matlib_common_get_address((dtype_t)bufferB->dtype, bufferB->data,(index_t)offsetB);
            if(pDataX==NULL) {
                zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type of label number.", 0);
                return;
            }
            int32_t errcode = rindow_matlib_i_gather(reverse,addMode,(index_t)n,(index_t)k,(index_t)numClass,(dtype_t)bufferX->dtype,pDataX,(dtype_t)bufferA->dtype,pDataA,pDataB);
            if(errcode) {
                if(errcode == RINDOW_MATLIB_E_UNSUPPORTED_DATA_TYPE) {
                    zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
                    return;
                } else if(errcode == RINDOW_MATLIB_E_PERM_OUT_OF_RANGE) {
                    zend_throw_exception(spl_ce_RuntimeException, "Label number is out of bounds.", 0);
                    return;
                } else {
                    zend_throw_exception(spl_ce_RuntimeException, "Unknown error.", 0);
                    return;
                }
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
            void *pDataX = rindow_matlib_common_get_address((dtype_t)bufferX->dtype, bufferX->data,(index_t)offsetX);
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataB,bufferB,offsetB)
            if(pDataX==NULL) {
                zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type of label number.", 0);
                return;
            }
            int32_t errcode = rindow_matlib_s_reducegather(reverse,addMode,(index_t)m,(index_t)n,(index_t)numClass,(dtype_t)bufferX->dtype,pDataX,pDataA,pDataB);
            if(errcode) {
                if(errcode == RINDOW_MATLIB_E_UNSUPPORTED_DATA_TYPE) {
                    zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type of label number.", 0);
                    return;
                } else if(errcode == RINDOW_MATLIB_E_PERM_OUT_OF_RANGE) {
                    zend_throw_exception(spl_ce_RuntimeException, "Label number is out of bounds.", 0);
                    return;
                } else {
                    zend_throw_exception_ex(spl_ce_InvalidArgumentException, 0, "Unknown error.(%d)", errcode);
                    return;
                }
            }
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            void *pDataX = rindow_matlib_common_get_address((dtype_t)bufferX->dtype, bufferX->data,(index_t)offsetX);
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataB,bufferB,offsetB)
            if(pDataX==NULL) {
                zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type of label number.", 0);
                return;
            }
            int32_t errcode = rindow_matlib_d_reducegather(reverse,addMode,(index_t)m,(index_t)n,(index_t)numClass,(dtype_t)bufferX->dtype,pDataX,pDataA,pDataB);
            if(errcode) {
                if(errcode == RINDOW_MATLIB_E_UNSUPPORTED_DATA_TYPE) {
                    zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type of label number.", 0);
                    return;
                } else if(errcode == RINDOW_MATLIB_E_PERM_OUT_OF_RANGE) {
                    zend_throw_exception(spl_ce_RuntimeException, "Label number is out of bounds.", 0);
                    return;
                } else {
                    zend_throw_exception_ex(spl_ce_InvalidArgumentException, 0, "Unknown error.: %d", errcode);
                    return;
                }
            }
            break;
        }
        default: {
            if(!php_rindow_openblas_common_dtype_is_int(bufferA->dtype)&&
                !php_rindow_openblas_common_dtype_is_bool(bufferA->dtype)) {
                zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
                return;
            }
            void *pDataX = rindow_matlib_common_get_address((dtype_t)bufferX->dtype, bufferX->data,(index_t)offsetX);
            void *pDataA = rindow_matlib_common_get_address((dtype_t)bufferA->dtype, bufferA->data,(index_t)offsetA);
            void *pDataB = rindow_matlib_common_get_address((dtype_t)bufferB->dtype, bufferB->data,(index_t)offsetB);
            if(pDataX==NULL) {
                zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type of label number.", 0);
                return;
            }
            int32_t errcode = rindow_matlib_i_reducegather(reverse,addMode,(index_t)m,(index_t)n,(index_t)numClass,(dtype_t)bufferX->dtype,pDataX,(dtype_t)bufferA->dtype,pDataA,pDataB);
            if(errcode) {
                if(errcode == RINDOW_MATLIB_E_UNSUPPORTED_DATA_TYPE) {
                    zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
                    return;
                } else if(errcode == RINDOW_MATLIB_E_PERM_OUT_OF_RANGE) {
                    zend_throw_exception(spl_ce_RuntimeException, "Label number is out of bounds.", 0);
                    return;
                } else {
                    zend_throw_exception(spl_ce_RuntimeException, "Unknown error.", 0);
                    return;
                }
            }
            break;
        }
    }
}
/* }}} */
