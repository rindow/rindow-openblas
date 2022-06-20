/*
   X(m) := sum( A(m,n) )

   Method Rindow\OpenBLAS\Math::
    public function reduceSum(
        int $m,
        int $n,
        int $k,
        Buffer $A, int $offsetA,
        Buffer $B, int $offsetB ) : void
 {{{ */
static PHP_METHOD(Math, reduceSum)
{
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    php_interop_polite_math_matrix_linear_buffer_t* bufferB;
    zend_long m;
    zend_long n;
    zend_long k;
    zval* a=NULL;
    zend_long offsetA;
    zval* b=NULL;
    zend_long offsetB;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 7, 7)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_LONG(k)
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
    if(php_rindow_openblas_assert_shape_parameter(
        "k", k)) {
        return;
    }
    // Check Buffer A
    bufferA = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_buffer_type(bufferA,"a")) {
        return;
    }
    if(offsetA<0) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Argument offsetA must be greater than or equals 0.", 0);
        return;
    }
    if(offsetA+m*n*k>bufferA->size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Matrix specification too large for bufferA.", 0);
        return;
    }

    // Check Buffer B
    bufferB = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(b);
    if(php_rindow_openblas_assert_buffer_type(bufferB,"b")) {
        return;
    }
    if(offsetB<0) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Argument offsetB must be greater than or equals 0.", 0);
        return;
    }
    if(offsetB+m*k>bufferB->size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Matrix specification too large for bufferB.", 0);
        return;
    }

    // Check Buffer A and B
    if(bufferA->dtype!=bufferB->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for A and B", 0);
        return;
    }

    switch (bufferA->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *a = &(((float *)bufferA->data)[offsetA]);
                float *b = &(((float *)bufferB->data)[offsetB]);
                zend_long idxA = 0;
                zend_long idxB = 0;
                zend_long ldA = n*k;
                zend_long ldB = k;
                for(zend_long i=0; i<m; i++,idxA+=ldA,idxB+=ldB) {
                    for(zend_long j=0; j<k; j++) {
                        b[idxB+j] = s_sum(n,&a[idxA+j],k);
                    }
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *a = &(((double *)bufferA->data)[offsetA]);
                double *b = &(((double *)bufferB->data)[offsetB]);
                zend_long idxA = 0;
                zend_long idxB = 0;
                zend_long ldA = n*k;
                zend_long ldB = k;
                for(zend_long i=0; i<m; i++,idxA+=ldA,idxB+=ldB) {
                    for(zend_long j=0; j<k; j++) {
                        b[idxB+j] = d_sum(n,&a[idxA+j],k);
                    }
                }
            }
            break;
        default:
            zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
            return;
    }
}
/* }}} */

/*
   X(m) := max( A(m,n) )

   Method Rindow\OpenBLAS\Math::
    public function reduceMax(
        int $m,
        int $n,
        int $k,
        Buffer $A, int $offsetA,
        Buffer $B, int $offsetB ) : void
 {{{ */
static PHP_METHOD(Math, reduceMax)
{
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    php_interop_polite_math_matrix_linear_buffer_t* bufferB;
    zend_long m;
    zend_long n;
    zend_long k;
    zval* a=NULL;
    zend_long offsetA;
    zval* b=NULL;
    zend_long offsetB;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 7, 7)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_LONG(k)
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
    if(php_rindow_openblas_assert_shape_parameter(
        "k", k)) {
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
    if(offsetA+m*n*k>bufferA->size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Matrix specification too large for bufferA.", 0);
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
    if(offsetB+m*k>bufferB->size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Matrix specification too large for bufferB.", 0);
        return;
    }

    // Check Buffer A and B
    if(bufferA->dtype!=bufferB->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for A and B", 0);
        return;
    }

    switch (bufferA->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *a = &(((float *)bufferA->data)[offsetA]);
                float *b = &(((float *)bufferB->data)[offsetB]);
                zend_long idxA = 0;
                zend_long idxB = 0;
                zend_long ldA = n*k;
                zend_long ldB = k;
                for(zend_long i=0; i<m; i++,idxA+=ldA,idxB+=ldB) {
                    for(zend_long j=0; j<k; j++) {
                        b[idxB+j] = s_max(n,&a[idxA+j],k);
                    }
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *a = &(((double *)bufferA->data)[offsetA]);
                double *b = &(((double *)bufferB->data)[offsetB]);
                zend_long idxA = 0;
                zend_long idxB = 0;
                zend_long ldA = n*k;
                zend_long ldB = k;
                for(zend_long i=0; i<m; i++,idxA+=ldA,idxB+=ldB) {
                    for(zend_long j=0; j<k; j++) {
                        b[idxB+j] = d_max(n,&a[idxA+j],k);
                    }
                }
            }
            break;
        default:
            zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
            return;
    }
}
/* }}} */

/*
   X(m) := argmax( A(m,n) )

   Method Rindow\OpenBLAS\Math::
    public function reduceArgMax(
        int $m,
        int $n,
        int $k,
        Buffer $A, int $offsetA,
        Buffer $B, int $offsetB ) : void
 {{{ */
static PHP_METHOD(Math, reduceArgMax)
{
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    php_interop_polite_math_matrix_linear_buffer_t* bufferB;
    zend_long m;
    zend_long n;
    zend_long k;
    zval* a=NULL;
    zend_long offsetA;
    zval* b=NULL;
    zend_long offsetB;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 7, 7)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_LONG(k)
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
    if(php_rindow_openblas_assert_shape_parameter(
        "k", k)) {
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
    if(offsetA+m*n*k>bufferA->size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Matrix specification too large for bufferA.", 0);
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
    if(offsetB+m*k>bufferB->size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Matrix specification too large for bufferB.", 0);
        return;
    }

    switch (bufferA->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *a = &(((float *)bufferA->data)[offsetA]);
                zend_long value;
                zend_long idxA = 0;
                zend_long idxB = 0;
                zend_long ldA = n*k;
                zend_long ldB = k;
                for(zend_long i=0; i<m; i++,idxA+=ldA,idxB+=ldB) {
                    for(zend_long j=0; j<k; j++) {
                        value = s_argmax(n,&a[idxA+j],k);
                        rindow_openblas_math_set_integer(bufferB->dtype, bufferB->data, offsetB, 1, idxB+j, value);
                    }
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *a = &(((double *)bufferA->data)[offsetA]);
                zend_long value;
                zend_long idxA = 0;
                zend_long idxB = 0;
                zend_long ldA = n*k;
                zend_long ldB = k;
                for(zend_long i=0; i<m; i++,idxA+=ldA,idxB+=ldB) {
                    for(zend_long j=0; j<k; j++) {
                        value = d_argmax(n,&a[idxA+j],k);
                        int rc = rindow_openblas_math_set_integer(bufferB->dtype, bufferB->data, offsetB, 1, idxB+j, value);
                        if(rc) {
                            zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
                            return;
                        }
                    }
                }
            }
            break;
        default:
            zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
            return;
    }
}
/* }}} */
