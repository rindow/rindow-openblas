/*
   B(n,repeats,k) := A(n,k)

   Method Rindow\OpenBLAS\Math::
   public function repeat(
       int $m,
       int $k,
       int $repeats,
       Buffer $A, int $offsetA,
       Buffer $B, int $offsetB
       ) : void
 {{{ */
static PHP_METHOD(Math, repeat)
{
    zend_long m;
    zend_long k;
    zend_long repeats;
    zval* a;
    zend_long offsetA;
    zval* b;
    zend_long offsetB;
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    php_interop_polite_math_matrix_linear_buffer_t* bufferB;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 7, 7)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(k)
        Z_PARAM_LONG(repeats)
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
        "k", k)) {
        return;
    }
    if(repeats<=0) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Argument repeats must be greater than or equal 0.", 0);
        return;
    }

    // Check Buffer A
    bufferA = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_buffer_type(bufferA,"a")) {
        return;
    }
    if(offsetA+m*k > bufferA->size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Matrix A specification too large for buffer.", 0);
        return;
    }

    // Check Buffer B
    bufferB = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(b);
    if(php_rindow_openblas_assert_buffer_type(bufferB,"b")) {
        return;
    }
    if(offsetB+m*repeats*k > bufferB->size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Matrix B specification too large for buffer.", 0);
        return;
    }

    // Check Buffer A and Y
    if(bufferA->dtype!=bufferB->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for A and B", 0);
        return;
    }

    switch (bufferA->dtype) {
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataB,bufferB,offsetB)
            rindow_matlib_s_repeat((index_t)m,(index_t)k,(index_t)repeats,pDataA,pDataB);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataB,bufferB,offsetB)
            rindow_matlib_d_repeat((index_t)m,(index_t)k,(index_t)repeats,pDataA,pDataB);
            break;
        }
        default:{
            if(!php_rindow_openblas_common_dtype_is_int(bufferA->dtype)&&
                !php_rindow_openblas_common_dtype_is_bool(bufferA->dtype)) {
                zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
                return;
            }
            void *pDataA = rindow_matlib_common_get_address((dtype_t)bufferA->dtype, bufferA->data,(index_t)offsetA);
            void *pDataB = rindow_matlib_common_get_address((dtype_t)bufferB->dtype, bufferB->data,(index_t)offsetB);
            rindow_matlib_i_repeat((index_t)m,(index_t)k,(index_t)repeats,(dtype_t)bufferA->dtype,pDataA,pDataB);
            break;
        }
    }
}
/* }}} */
