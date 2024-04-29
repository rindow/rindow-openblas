#include <ext/standard/php_rand.h>


/*
   X(i) := rand(seed)

   Method Rindow\OpenBLAS\Math::
    public function randomUniform(
        int $n,
        Buffer $X, int $offsetX, int $incX,
        float $low,
        float $high,
        int $seed
        ) : void
 {{{ */

static PHP_METHOD(Math, randomUniform)
{
    php_interop_polite_math_matrix_linear_buffer_t* bufferX;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zval* low_val=NULL;
    zval* high_val=NULL;
    zend_long low_int;
    zend_long high_int;
    double low_float;
    double high_float;
    zend_long seed;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 7, 7)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_ZVAL(low_val)
        Z_PARAM_ZVAL(high_val)
        Z_PARAM_LONG(seed)
    ZEND_PARSE_PARAMETERS_END();

    // Check Buffer X
    bufferX = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(bufferX,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", bufferX,n,offsetX,incX)) {
        return;
    }

    if(php_rindow_openblas_common_dtype_is_float(bufferX->dtype)) {
        if(php_rindow_openblas_val2float(low_val,&low_float,"low must be float or int")) {
            return;
        }
        if(php_rindow_openblas_val2float(high_val,&high_float,"high must be float or int")) {
            return;
        }
    } else if(php_rindow_openblas_common_dtype_is_int(bufferX->dtype)) {
        if(php_rindow_openblas_val2int(low_val,&low_int,"low must be float or int")) {
            return;
        }
        if(php_rindow_openblas_val2int(high_val,&high_int,"high must be float or int")) {
            return;
        }
    } else {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
        return;
    }

    switch(bufferX->dtype) {
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,bufferX,offsetX)
            rindow_matlib_s_randomuniform((index_t)n,pDataX,(index_t)incX,(float)low_float,(float)high_float,(int32_t)seed);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,bufferX,offsetX)
            rindow_matlib_d_randomuniform((index_t)n,pDataX,(index_t)incX,(double)low_float,(double)high_float,(int32_t)seed);
            break;
        }
        default: {
            void *pDataX = rindow_matlib_common_get_address((dtype_t)bufferX->dtype, bufferX->data,(index_t)offsetX);
            rindow_matlib_i_randomuniform((index_t)n,(dtype_t)bufferX->dtype,pDataX,(index_t)incX,(int32_t)low_int,(int32_t)high_int,(int32_t)seed);
            break;
        }
    }
}
/* }}} */

/*
   X(i) := rand(seed)

   Method Rindow\OpenBLAS\Math::
    public function randomNormal(
        int $n,
        Buffer $X, int $offsetX, int $incX,
        float $mean,
        float $scale,
        int $seed
        ) : void
 {{{ */

static PHP_METHOD(Math, randomNormal)
{
    php_interop_polite_math_matrix_linear_buffer_t* bufferX;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zval* low_val=NULL;
    zval* high_val=NULL;
    double mean;
    double scale;
    zend_long seed;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 7, 7)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_DOUBLE(mean)
        Z_PARAM_DOUBLE(scale)
        Z_PARAM_LONG(seed)
    ZEND_PARSE_PARAMETERS_END();

    // Check Buffer X
    bufferX = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(bufferX,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", bufferX,n,offsetX,incX)) {
        return;
    }

    switch(bufferX->dtype) {
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,bufferX,offsetX)
            rindow_matlib_s_randomnormal((index_t)n,pDataX,(index_t)incX,(float)mean,(float)scale,(int32_t)seed);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,bufferX,offsetX)
            rindow_matlib_d_randomnormal((index_t)n,pDataX,(index_t)incX,(double)mean,(double)scale,(int32_t)seed);
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
   X(i) := rand(seed)

   Method Rindow\OpenBLAS\Math::
    public function randomSequence(
        int $n,
        int $size,
        Buffer $X, int $offsetX, int $incX
        int $seed
        ) : void
 {{{ */

static PHP_METHOD(Math, randomSequence)
{
    php_interop_polite_math_matrix_linear_buffer_t* bufferX;
    zend_long n;
    zend_long size;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long seed;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 6, 6)
        Z_PARAM_LONG(n)
        Z_PARAM_LONG(size)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_LONG(seed)
    ZEND_PARSE_PARAMETERS_END();

    // Check Buffer X
    bufferX = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(bufferX,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", bufferX,n,offsetX,incX)) {
        return;
    }
    if(n<size||size<1) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "size must be smaller then n or equal.", 0);
        return;
    }
    if(bufferX->dtype!=php_interop_polite_math_matrix_dtype_int64&&
        bufferX->dtype!=php_interop_polite_math_matrix_dtype_int32) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "dtype must be int32 or int64.", 0);
        return;
    }

    void *pDataX = rindow_matlib_common_get_address((dtype_t)bufferX->dtype, bufferX->data,(index_t)offsetX);
    rindow_matlib_i_randomsequence((index_t)n,(index_t)size,(dtype_t)bufferX->dtype,pDataX,(index_t)incX,(int32_t)seed);
}   
/* }}} */
