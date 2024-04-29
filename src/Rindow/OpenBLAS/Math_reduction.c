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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataB,bufferB,offsetB)
            rindow_matlib_s_reducesum((index_t)m,(index_t)n,(index_t)k,pDataA,pDataB);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataB,bufferB,offsetB)
            rindow_matlib_d_reducesum((index_t)m,(index_t)n,(index_t)k,pDataA,pDataB);
            break;
        }
        default: {
            zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
            return;
        }
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataB,bufferB,offsetB)
            rindow_matlib_s_reducemax((index_t)m,(index_t)n,(index_t)k,pDataA,pDataB);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataB,bufferB,offsetB)
            rindow_matlib_d_reducemax((index_t)m,(index_t)n,(index_t)k,pDataA,pDataB);
            break;
        }
        default: {
            zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
            return;
        }
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataA,bufferA,offsetA)
            void *pDataB = rindow_matlib_common_get_address((dtype_t)bufferB->dtype, bufferB->data,(index_t)offsetB);
            rindow_matlib_s_reduceargmax((index_t)m,(index_t)n,(index_t)k,pDataA,(dtype_t)bufferB->dtype,pDataB);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataA,bufferA,offsetA)
            void *pDataB = rindow_matlib_common_get_address((dtype_t)bufferB->dtype, bufferB->data,(index_t)offsetB);
            rindow_matlib_d_reduceargmax((index_t)m,(index_t)n,(index_t)k,pDataA,(dtype_t)bufferB->dtype,pDataB);
            break;
        }
        default: {
            zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
            return;
        }
    }
}
/* }}} */
