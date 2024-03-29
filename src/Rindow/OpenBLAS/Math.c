#include <php.h>
#include <Zend/zend_interfaces.h>
#include <Zend/zend_exceptions.h>
#include <ext/spl/spl_iterators.h>
#include <ext/spl/spl_exceptions.h>
#include <ext/standard/php_rand.h>
#include <cblas.h>
#include <stdint.h>
#include <Interop/Polite/Math/Matrix.h>


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "php_rindow_openblas.h"

static float s_max(zend_long n,float *x,zend_long incX)
{
    zend_long i;
    float a;
    a = x[0];
    for(i=1;i<n;i++) {
        // if NaN set NaN
        // Compatible with reduce_max of tensorflow 2.6
        if(!(a>=x[i*incX])) {
            a = x[i*incX];
        }
    }
    return a;
}

static double d_max(zend_long n,double *x,zend_long incX)
{
    zend_long i;
    double a;
    a = x[0];
    for(i=1;i<n;i++) {
        if(a<x[i*incX]) {
            a = x[i*incX];
        }
    }
    return a;
}

static zend_long s_argmax(zend_long n,float *x,zend_long incX)
{
    zend_long i;
    zend_long idx;
    float a;
    idx = 0;
    a = x[0];
    for(i=1;i<n;i++) {
        if(a<x[i*incX]) {
            idx = i;
            a = x[i*incX];
        }
    }
    return idx;
}

static zend_long d_argmax(zend_long n,double *x,zend_long incX)
{
    zend_long i;
    zend_long idx;
    double a;
    idx = 0;
    a = x[0];
    for(i=1;i<n;i++) {
        if(a<x[i*incX]) {
            idx = i;
            a = x[i*incX];
        }
    }
    return idx;
}

static float s_sum(zend_long n,float *x,zend_long incX)
{
    zend_long i;
    float a=0;
    for(i=0; i<n; i++) {
        a += x[i*incX];
    }
    return a;
}

static double d_sum(zend_long n,double *x,zend_long incX)
{
    zend_long i;
    double a=0;
    for(i=0; i<n; i++) {
        a += x[i*incX];
    }
    return a;
}


static zend_object_handlers rindow_openblas_math_object_handlers;

// destractor
static void php_rindow_openblas_math_free_object(zend_object* object)
{
    zend_object_std_dtor(object);
}

// constructor
static zend_object* php_rindow_openblas_math_create_object(zend_class_entry* class_type) /* {{{ */
{
    zend_object* intern = NULL;

    intern = (zend_object*)ecalloc(1, sizeof(zend_object) + zend_object_properties_size(class_type));

    zend_object_std_init(intern, class_type);
    object_properties_init(intern, class_type);

    intern->handlers = &rindow_openblas_math_object_handlers;

    return intern;
} /* }}} */

#define PHP_RINDOW_OPENBLAS_MATH_SUM_TEMPLATE(data_type) { \
    data_type  *pDataX; \
    pDataX = &(((data_type *)buffer->data)[offsetX]); \
    result = 0.0; \
    for (i=0; i<n; i++,pDataX+=incX) { \
        result += (data_type)*pDataX; \
    } \
}

int php_rindow_openblas_val2int(
    zval* val_value,
    zend_long* integer_value,
    char* message)
{
	switch(Z_TYPE_P(val_value)) {
	    case IS_LONG:
	        *integer_value = Z_LVAL_P(val_value);
	        break;
	    case IS_DOUBLE:
	        *integer_value = (zend_long)Z_DVAL_P(val_value);
	        break;
		default:
            zend_throw_exception(spl_ce_InvalidArgumentException, message, 0);
			return -1;
	}
	return 0;
}

int php_rindow_openblas_val2float(
    zval* val_value,
    double* float_value,
    char* message)
{
	switch(Z_TYPE_P(val_value)) {
	    case IS_LONG:
	        *float_value = (double)Z_LVAL_P(val_value);
	        break;
	    case IS_DOUBLE:
	        *float_value = Z_DVAL_P(val_value);
	        break;
		default:
            zend_throw_exception(spl_ce_InvalidArgumentException, message, 0);
			return -1;
	}
	return 0;
}

#define PHP_RINDOW_OPENBLAS_MATH_ADD_TEMPLATE(data_type) { \
    data_type  *pDataX; \
    data_type  *pDataY; \
    pDataX = (data_type *)values; \
    pDataY = (data_type *)target; \
    for (i=0; i<n; i++) { \
        *pDataY += *pDataX; \
        pDataX+=incValue; \
        pDataY+=incTarget; \
    } \
}
int php_rindow_openblas_math_add(
    zend_long n,
    zend_long dtype,
    void* values,
    zend_long incValue,
    void* target,
    zend_long incTarget
    )
{
    switch (dtype) {
        zend_long i;
        case php_interop_polite_math_matrix_dtype_float32:
            PHP_RINDOW_OPENBLAS_MATH_ADD_TEMPLATE(float)
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            PHP_RINDOW_OPENBLAS_MATH_ADD_TEMPLATE(double)
            break;
        case php_interop_polite_math_matrix_dtype_int8:
            PHP_RINDOW_OPENBLAS_MATH_ADD_TEMPLATE(int8_t)
            break;
        case php_interop_polite_math_matrix_dtype_uint8:
            PHP_RINDOW_OPENBLAS_MATH_ADD_TEMPLATE(uint8_t)
            break;
        case php_interop_polite_math_matrix_dtype_int16:
            PHP_RINDOW_OPENBLAS_MATH_ADD_TEMPLATE(int16_t)
            break;
        case php_interop_polite_math_matrix_dtype_uint16:
            PHP_RINDOW_OPENBLAS_MATH_ADD_TEMPLATE(uint16_t)
            break;
        case php_interop_polite_math_matrix_dtype_int32:
            PHP_RINDOW_OPENBLAS_MATH_ADD_TEMPLATE(int32_t)
            break;
        case php_interop_polite_math_matrix_dtype_uint32:
            PHP_RINDOW_OPENBLAS_MATH_ADD_TEMPLATE(uint32_t)
            break;
        case php_interop_polite_math_matrix_dtype_int64:
            PHP_RINDOW_OPENBLAS_MATH_ADD_TEMPLATE(int64_t)
            break;
        case php_interop_polite_math_matrix_dtype_uint64:
            PHP_RINDOW_OPENBLAS_MATH_ADD_TEMPLATE(uint64_t)
            break;
        default:
            zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
            return -1;
    }
    return 0;
}

#define PHP_RINDOW_OPENBLAS_MATH_COPY_TEMPLATE(data_type) { \
    data_type  *pDataX; \
    data_type  *pDataY; \
    pDataX = (data_type *)source; \
    pDataY = (data_type *)dest; \
    for (i=0; i<n; i++) { \
        *pDataY += *pDataX; \
        pDataX+=incSource; \
        pDataY+=incDest; \
    } \
}
int php_rindow_openblas_math_copy(
    zend_long n,
    zend_long dtype,
    void* source,
    zend_long incSource,
    void* dest,
    zend_long incDest
    )
{
    switch (dtype) {
        zend_long i;
        case php_interop_polite_math_matrix_dtype_float32:
            PHP_RINDOW_OPENBLAS_MATH_COPY_TEMPLATE(float)
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            PHP_RINDOW_OPENBLAS_MATH_COPY_TEMPLATE(double)
            break;
        case php_interop_polite_math_matrix_dtype_int8:
            PHP_RINDOW_OPENBLAS_MATH_COPY_TEMPLATE(int8_t)
            break;
        case php_interop_polite_math_matrix_dtype_uint8:
            PHP_RINDOW_OPENBLAS_MATH_COPY_TEMPLATE(uint8_t)
            break;
        case php_interop_polite_math_matrix_dtype_int16:
            PHP_RINDOW_OPENBLAS_MATH_COPY_TEMPLATE(int16_t)
            break;
        case php_interop_polite_math_matrix_dtype_uint16:
            PHP_RINDOW_OPENBLAS_MATH_COPY_TEMPLATE(uint16_t)
            break;
        case php_interop_polite_math_matrix_dtype_int32:
            PHP_RINDOW_OPENBLAS_MATH_COPY_TEMPLATE(int32_t)
            break;
        case php_interop_polite_math_matrix_dtype_uint32:
            PHP_RINDOW_OPENBLAS_MATH_COPY_TEMPLATE(uint32_t)
            break;
        case php_interop_polite_math_matrix_dtype_int64:
            PHP_RINDOW_OPENBLAS_MATH_COPY_TEMPLATE(int64_t)
            break;
        case php_interop_polite_math_matrix_dtype_uint64:
            PHP_RINDOW_OPENBLAS_MATH_COPY_TEMPLATE(uint64_t)
            break;
        default:
            zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
            return -1;
    }
    return 0;
}

/* Method Rindow\OpenBLAS\Math::
    public function sum(
        int $n,
        Buffer $X, int $offsetX, int $incX ) : float
 {{{ */
static PHP_METHOD(Math, sum)
{
    php_interop_polite_math_matrix_linear_buffer_t* buffer;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    double result;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 4, 4)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    buffer = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(buffer,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        zend_long i;
        case php_interop_polite_math_matrix_dtype_float32:
            PHP_RINDOW_OPENBLAS_MATH_SUM_TEMPLATE(float)
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            PHP_RINDOW_OPENBLAS_MATH_SUM_TEMPLATE(double)
            break;
        case php_interop_polite_math_matrix_dtype_int8:
            PHP_RINDOW_OPENBLAS_MATH_SUM_TEMPLATE(int8_t)
            break;
        case php_interop_polite_math_matrix_dtype_uint8:
            PHP_RINDOW_OPENBLAS_MATH_SUM_TEMPLATE(uint8_t)
            break;
        case php_interop_polite_math_matrix_dtype_int16:
            PHP_RINDOW_OPENBLAS_MATH_SUM_TEMPLATE(int16_t)
            break;
        case php_interop_polite_math_matrix_dtype_uint16:
            PHP_RINDOW_OPENBLAS_MATH_SUM_TEMPLATE(uint16_t)
            break;
        case php_interop_polite_math_matrix_dtype_int32:
            PHP_RINDOW_OPENBLAS_MATH_SUM_TEMPLATE(int32_t)
            break;
        case php_interop_polite_math_matrix_dtype_uint32:
            PHP_RINDOW_OPENBLAS_MATH_SUM_TEMPLATE(uint32_t)
            break;
        case php_interop_polite_math_matrix_dtype_int64:
            PHP_RINDOW_OPENBLAS_MATH_SUM_TEMPLATE(int64_t)
            break;
        case php_interop_polite_math_matrix_dtype_uint64:
            PHP_RINDOW_OPENBLAS_MATH_SUM_TEMPLATE(uint64_t)
            break;
        case php_interop_polite_math_matrix_dtype_bool:
            {
                uint8_t *pBoolX;
                pBoolX = &(((uint8_t *)buffer->data)[offsetX]);
                result = 0.0;
                for (i=0; i<n; i++,pBoolX+=incX) {
                    if(*pBoolX!=0) {
                        result += 1;
                    }
                }
            }
            break;
        default:
            zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
            return;
    }
    RETURN_DOUBLE(result);
}
/* }}} */

/* Method Rindow\OpenBLAS\Math::
    public function imax(
        int $n,
        Buffer $X, int $offsetX, int $incX) : int
 {{{ */
static PHP_METHOD(Math, imax)
{
    php_interop_polite_math_matrix_linear_buffer_t* buffer;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    float  *pFloatX;
    double *pDoubleX;
    float  floatMax;
    double doubleMax;
    zend_long resultIdx;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 4, 4)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    buffer = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(buffer,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        zend_long i;
        case php_interop_polite_math_matrix_dtype_float32:
            pFloatX = &(((float *)buffer->data)[offsetX]);
            floatMax = *pFloatX;
            pFloatX += incX;
            resultIdx = 0;
            for (i=1; i<n; i++,pFloatX+=incX) {
                if(floatMax < *pFloatX || isnan(floatMax)) {
                    floatMax = *pFloatX;
                    resultIdx = i;
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            pDoubleX = &(((double *)buffer->data)[offsetX]);
            doubleMax = *pDoubleX;
            pDoubleX += incX;
            resultIdx = 0;
            for (i=1; i<n; i++,pDoubleX+=incX) {
                if(doubleMax < *pDoubleX || isnan(doubleMax)) {
                    doubleMax = *pDoubleX;
                    resultIdx = i;
                }
            }
            break;
        default:
            zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
            return;
    }
    RETURN_LONG(resultIdx);
}
/* }}} */

/* Method Rindow\OpenBLAS\Math::
    public function imin(
        int $n,
        Buffer $X, int $offsetX, int $incX) : int
 {{{ */
static PHP_METHOD(Math, imin)
{
    php_interop_polite_math_matrix_linear_buffer_t* buffer;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    float  *pFloatX;
    double *pDoubleX;
    float  floatMin;
    double doubleMin;
    zend_long resultIdx;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 4, 4)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    buffer = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(buffer,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        zend_long i;
        case php_interop_polite_math_matrix_dtype_float32:
            pFloatX = &(((float *)buffer->data)[offsetX]);
            floatMin = *pFloatX;
            pFloatX += incX;
            resultIdx = 0;
            for (i=1; i<n; i++,pFloatX+=incX) {
                if(floatMin > *pFloatX) {
                    floatMin = *pFloatX;
                    resultIdx = i;
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            pDoubleX = &(((double *)buffer->data)[offsetX]);
            doubleMin = *pDoubleX;
            pDoubleX += incX;
            resultIdx = 0;
            for (i=1; i<n; i++,pDoubleX+=incX) {
                if(doubleMin > *pDoubleX) {
                    doubleMin = *pDoubleX;
                    resultIdx = i;
                }
            }
            break;
        default:
            zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
            return;
    }
    RETURN_LONG(resultIdx);
}
/* }}} */

/*
   X := a*X + b

   Method Rindow\OpenBLAS\Math::
    public function increment(
        int $n,
        float $alpha,
        Buffer $X, int $offsetX, int $incX,
        float $beta) : void
 {{{ */
static PHP_METHOD(Math, increment)
{
    php_interop_polite_math_matrix_linear_buffer_t* buffer;
    zend_long n;
    double alpha;
    double beta;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long i;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 6, 6)
        Z_PARAM_LONG(n)
        Z_PARAM_DOUBLE(alpha)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_DOUBLE(beta)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    buffer = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(buffer,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    x[i*incX] = (float)alpha * x[i*incX] + (float)beta;
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *x = &(((double *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    x[i*incX] = (double)alpha * x[i*incX] + (double)beta;
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
   X := 1 / (a*X + b)

   Method Rindow\OpenBLAS\Math::
    public function reciprocal(
        int $n,
        float $alpha,
        Buffer $X, int $offsetX, int $incX,
        float $beta) : void
 {{{ */
static PHP_METHOD(Math, reciprocal)
{
    php_interop_polite_math_matrix_linear_buffer_t* buffer;
    zend_long n;
    double alpha;
    double beta;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long i;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 6, 6)
        Z_PARAM_LONG(n)
        Z_PARAM_DOUBLE(alpha)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_DOUBLE(beta)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    buffer = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(buffer,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    float t;
                    t = (float)alpha * x[i*incX] + (float)beta;
                    // *** CAUTION ***
                    // disable checking for INFINITY values
                    //if(t==0.0) {
                    //    zend_throw_exception(spl_ce_RuntimeException, "Zero divide.", 0);
                    //    return;
                    //}
                    x[i*incX] = 1 / t;
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *x = &(((double *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    double t;
                    t = (double)alpha * x[i*incX] + (double)beta;
                    // *** CAUTION ***
                    // disable checking for INFINITY values
                    //if(t==0.0) {
                    //    zend_throw_exception(spl_ce_RuntimeException, "Zero divide.", 0);
                    //    return;
                    //}
                    x[i*incX] = 1 / t;
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
    A[m,n] := A[m,n] (A[m,n] >  X[n])
    A[m,n] := X[n]   (A[m,n] <= X[n])

   Method Rindow\OpenBLAS\Math::
    public function maximum(
        int $m,
        int $n,
        Buffer $A, int $offsetA, int $ldA,
        Buffer $X, int $offsetX, int $incX,
        ) : void
 {{{ */
static PHP_METHOD(Math, maximum)
{
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    php_interop_polite_math_matrix_linear_buffer_t* bufferX;
    zend_long m;
    zend_long n;
    zval* a=NULL;
    zend_long offsetA;
    zend_long ldA;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 8, 8)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(a) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "m", m)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    bufferA = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_buffer_type(bufferA,"a")) {
        return;
    }
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        "A", bufferA,m,n,offsetA,ldA)) {
        return;
    }

    bufferX = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(bufferX,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", bufferX,n,offsetX,incX)) {
        return;
    }
    // Check Buffer X and A
    if(bufferA->dtype!=bufferX->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for A and X", 0);
        return;
    }
    switch (bufferA->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *a = &(((float *)bufferA->data)[offsetA]);
                float *x = &(((float *)bufferX->data)[offsetX]);
                float value;
                for(zend_long i=0;i<m;i++) {
                    for(zend_long j=0;j<n;j++) {
                        value = x[j*incX];
                        if(isnan(value)) {
                            a[i*ldA+j] = value;
                        } else {
                            // *** CAUTION ***
                            // if NaN then don't set alpha
                            if(a[i*ldA+j] < value) {
                                a[i*ldA+j] = value;
                            }
                        }
                    }
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *a = &(((double *)bufferA->data)[offsetA]);
                double *x = &(((double *)bufferX->data)[offsetX]);
                double value;
                for(zend_long i=0;i<m;i++) {
                    for(zend_long j=0;j<n;j++) {
                        value = x[j*incX];
                        if(isnan(value)) {
                            a[i*ldA+j] = value;
                        } else {
                            // *** CAUTION ***
                            // if NaN then don't set alpha
                            if(a[i*ldA+j] < value) {
                                a[i*ldA+j] = value;
                            }
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


/*
    A[m,n] := A[m,n] (A[m,n] <  X[n])
    A[m,n] := X[n]   (A[m,n] >= X[n])

   Method Rindow\OpenBLAS\Math::
    public function minimum(
        int $m,
        int $n,
        Buffer $A, int $offsetA, int $ldA,
        Buffer $X, int $offsetX, int $incX,
        ) : void
 {{{ */
static PHP_METHOD(Math, minimum)
{
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    php_interop_polite_math_matrix_linear_buffer_t* bufferX;
    zend_long m;
    zend_long n;
    zval* a=NULL;
    zend_long offsetA;
    zend_long ldA;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 8, 8)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(a) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "m", m)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    bufferA = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_buffer_type(bufferA,"a")) {
        return;
    }
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        "A", bufferA,m,n,offsetA,ldA)) {
        return;
    }

    bufferX = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(bufferX,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", bufferX,n,offsetX,incX)) {
        return;
    }
    // Check Buffer X and A
    if(bufferA->dtype!=bufferX->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for A and X", 0);
        return;
    }
    switch (bufferA->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *a = &(((float *)bufferA->data)[offsetA]);
                float *x = &(((float *)bufferX->data)[offsetX]);
                float value;
                for(zend_long i=0;i<m;i++) {
                    for(zend_long j=0;j<n;j++) {
                        value = x[j*incX];
                        if(isnan(value)) {
                            a[i*ldA+j] = value;
                        } else {
                            // *** CAUTION ***
                            // if NaN then don't set alpha
                            if(a[i*ldA+j] > value) {
                                a[i*ldA+j] = value;
                            }
                        }
                    }
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *a = &(((double *)bufferA->data)[offsetA]);
                double *x = &(((double *)bufferX->data)[offsetX]);
                double value;
                for(zend_long i=0;i<m;i++) {
                    for(zend_long j=0;j<n;j++) {
                        value = x[j*incX];
                        if(isnan(value)) {
                            a[i*ldA+j] = value;
                        } else {
                            // *** CAUTION ***
                            // if NaN then don't set alpha
                            if(a[i*ldA+j] > value) {
                                a[i*ldA+j] = value;
                            }
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

/*
    A[m,n] := 1 (A[m,n] >  X[n])
    A[m,n] := 0 (A[m,n] <= X[n])

   Method Rindow\OpenBLAS\Math::
    public function greater(
        int $m,
        int $n,
        Buffer $A, int $offsetA, int $ldA,
        Buffer $X, int $offsetX, int $incX,
        ) : void
 {{{ */
static PHP_METHOD(Math, greater)
{
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    php_interop_polite_math_matrix_linear_buffer_t* bufferX;
    zend_long m;
    zend_long n;
    zval* a=NULL;
    zend_long offsetA;
    zend_long ldA;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 8, 8)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(a) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "m", m)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    bufferA = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_buffer_type(bufferA,"a")) {
        return;
    }
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        "A", bufferA,m,n,offsetA,ldA)) {
        return;
    }

    bufferX = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(bufferX,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", bufferX,n,offsetX,incX)) {
        return;
    }
    // Check Buffer X and A
    if(bufferA->dtype!=bufferX->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for A and X", 0);
        return;
    }
    switch (bufferA->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *a = &(((float *)bufferA->data)[offsetA]);
                float *x = &(((float *)bufferX->data)[offsetX]);
                for(zend_long i=0;i<m;i++) {
                    for(zend_long j=0;j<n;j++) {
                        // *** CAUTION ***
                        // if NaN set 0.0
                        // if equal set 0.0
                        if(a[i*ldA+j] > x[j*incX]) {
                            a[i*ldA+j] = 1.0;
                        } else {
                            a[i*ldA+j] = 0.0;
                        }
                    }
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *a = &(((double *)bufferA->data)[offsetA]);
                double *x = &(((double *)bufferX->data)[offsetX]);
                for(zend_long i=0;i<m;i++) {
                    for(zend_long j=0;j<n;j++) {
                        // *** CAUTION ***
                        // if NaN set 0.0
                        // if equal set 0.0
                        if(a[i*ldA+j] > x[j*incX]) {
                            a[i*ldA+j] = 1.0;
                        } else {
                            a[i*ldA+j] = 0.0;
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

/*
    A[m,n] := 1 (A[m,n] >= X[n])
    A[m,n] := 0 (A[m,n] <  X[n])

   Method Rindow\OpenBLAS\Math::
    public function greaterEqual(
        int $m,
        int $n,
        Buffer $A, int $offsetA, int $ldA,
        Buffer $X, int $offsetX, int $incX,
        ) : void
 {{{ */
static PHP_METHOD(Math, greaterEqual)
{
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    php_interop_polite_math_matrix_linear_buffer_t* bufferX;
    zend_long m;
    zend_long n;
    zval* a=NULL;
    zend_long offsetA;
    zend_long ldA;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 8, 8)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(a) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "m", m)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    bufferA = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_buffer_type(bufferA,"a")) {
        return;
    }
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        "A", bufferA,m,n,offsetA,ldA)) {
        return;
    }

    bufferX = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(bufferX,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", bufferX,n,offsetX,incX)) {
        return;
    }
    // Check Buffer X and A
    if(bufferA->dtype!=bufferX->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for A and X", 0);
        return;
    }
    switch (bufferA->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *a = &(((float *)bufferA->data)[offsetA]);
                float *x = &(((float *)bufferX->data)[offsetX]);
                for(zend_long i=0;i<m;i++) {
                    for(zend_long j=0;j<n;j++) {
                        // *** CAUTION ***
                        // if NaN set 0.0
                        // if equal set 1.0
                        if(a[i*ldA+j] >= x[j*incX]) {
                            a[i*ldA+j] = 1.0;
                        } else {
                            a[i*ldA+j] = 0.0;
                        }
                    }
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *a = &(((double *)bufferA->data)[offsetA]);
                double *x = &(((double *)bufferX->data)[offsetX]);
                for(zend_long i=0;i<m;i++) {
                    for(zend_long j=0;j<n;j++) {
                        // *** CAUTION ***
                        // if NaN set 0.0
                        // if equal set 1.0
                        if(a[i*ldA+j] >= x[j*incX]) {
                            a[i*ldA+j] = 1.0;
                        } else {
                            a[i*ldA+j] = 0.0;
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

/*
    A[m,n] := 1 (A[m,n] <  X[n])
    A[m,n] := 0 (A[m,n] >= X[n])

   Method Rindow\OpenBLAS\Math::
    public function less(
        int $m,
        int $n,
        Buffer $A, int $offsetA, int $ldA,
        Buffer $X, int $offsetX, int $incX,
        ) : void
 {{{ */
static PHP_METHOD(Math, less)
{
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    php_interop_polite_math_matrix_linear_buffer_t* bufferX;
    zend_long m;
    zend_long n;
    zval* a=NULL;
    zend_long offsetA;
    zend_long ldA;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 8, 8)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(a) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "m", m)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    bufferA = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_buffer_type(bufferA,"a")) {
        return;
    }
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        "A", bufferA,m,n,offsetA,ldA)) {
        return;
    }

    bufferX = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(bufferX,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", bufferX,n,offsetX,incX)) {
        return;
    }
    // Check Buffer X and A
    if(bufferA->dtype!=bufferX->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for A and X", 0);
        return;
    }
    switch (bufferA->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *a = &(((float *)bufferA->data)[offsetA]);
                float *x = &(((float *)bufferX->data)[offsetX]);
                for(zend_long i=0;i<m;i++) {
                    for(zend_long j=0;j<n;j++) {
                        // *** CAUTION ***
                        // if NaN set 0.0
                        // if equal set 0.0
                        if(a[i*ldA+j] < x[j*incX]) {
                            a[i*ldA+j] = 1.0;
                        } else {
                            a[i*ldA+j] = 0.0;
                        }
                    }
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *a = &(((double *)bufferA->data)[offsetA]);
                double *x = &(((double *)bufferX->data)[offsetX]);
                for(zend_long i=0;i<m;i++) {
                    for(zend_long j=0;j<n;j++) {
                        // *** CAUTION ***
                        // if NaN set 0.0
                        // if equal set 1.0
                        if(a[i*ldA+j] < x[j*incX]) {
                            a[i*ldA+j] = 1.0;
                        } else {
                            a[i*ldA+j] = 0.0;
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


/*
    A[m,n] := 1 (A[m,n] <= X[n])
    A[m,n] := 0 (A[m,n] >  X[n])

   Method Rindow\OpenBLAS\Math::
    public function lessEqual(
        int $m,
        int $n,
        Buffer $A, int $offsetA, int $ldA,
        Buffer $X, int $offsetX, int $incX,
        ) : void
 {{{ */
static PHP_METHOD(Math, lessEqual)
{
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    php_interop_polite_math_matrix_linear_buffer_t* bufferX;
    zend_long m;
    zend_long n;
    zval* a=NULL;
    zend_long offsetA;
    zend_long ldA;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 8, 8)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(a) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "m", m)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    bufferA = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_buffer_type(bufferA,"a")) {
        return;
    }
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        "A", bufferA,m,n,offsetA,ldA)) {
        return;
    }

    bufferX = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(bufferX,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", bufferX,n,offsetX,incX)) {
        return;
    }
    // Check Buffer X and A
    if(bufferA->dtype!=bufferX->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for A and X", 0);
        return;
    }
    switch (bufferA->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *a = &(((float *)bufferA->data)[offsetA]);
                float *x = &(((float *)bufferX->data)[offsetX]);
                for(zend_long i=0;i<m;i++) {
                    for(zend_long j=0;j<n;j++) {
                        // *** CAUTION ***
                        // if NaN set 0.0
                        // if equal set 0.0
                        if(a[i*ldA+j] <= x[j*incX]) {
                            a[i*ldA+j] = 1.0;
                        } else {
                            a[i*ldA+j] = 0.0;
                        }
                    }
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *a = &(((double *)bufferA->data)[offsetA]);
                double *x = &(((double *)bufferX->data)[offsetX]);
                for(zend_long i=0;i<m;i++) {
                    for(zend_long j=0;j<n;j++) {
                        // *** CAUTION ***
                        // if NaN set 0.0
                        // if equal set 1.0
                        if(a[i*ldA+j] <= x[j*incX]) {
                            a[i*ldA+j] = 1.0;
                        } else {
                            a[i*ldA+j] = 0.0;
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


/*
   A(i) := X(i) * A(i)

   Method Rindow\OpenBLAS\Math::
    public function multiply(
        bool $trans,
        int $m,
        int $n,
        Buffer $X, int $offsetX, int $incX,
        Buffer $A, int $offsetA, int $ldA ) : void
 {{{ */
static PHP_METHOD(Math, multiply)
{
    php_interop_polite_math_matrix_linear_buffer_t* bufferX;
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    zend_bool trans;
    zend_long m;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zval* a=NULL;
    zend_long offsetA;
    zend_long ldA;
    zend_long rows,cols;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 9, 9)
        Z_PARAM_BOOL(trans)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_OBJECT(a) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "m", m)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    // Check Buffer X
    bufferX = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(bufferX,"x")) {
        return;
    }
    if(!trans) {
        rows = m; cols = n;
    } else {
        rows = n; cols = m;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", bufferX,cols,offsetX,incX)) {
        return;
    }

    // Check Buffer A
    bufferA = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_buffer_type(bufferA,"a")) {
        return;
    }
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        "A", bufferA,m,n,offsetA,ldA)) {
        return;
    }

    // Check Buffer X and A
    if(bufferX->dtype!=bufferA->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for X and A", 0);
        return;
    }

    switch (bufferX->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *x = &(((float *)bufferX->data)[offsetX]);
                float *a = &(((float *)bufferA->data)[offsetA]);
                zend_long i,j,incAj,incAi;
                if(!trans) { incAj = ldA; incAi = 1;}
                else       { incAj = 1;   incAi = ldA;}
                for(j=0; j<rows; j++) {
                    for(i=0; i<cols; i++) {
                        a[j*incAj+i*incAi] = x[i*incX] * a[j*incAj+i*incAi];
                    }
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *x = &(((double *)bufferX->data)[offsetX]);
                double *a = &(((double *)bufferA->data)[offsetA]);
                zend_long i,j,incAj,incAi;
                if(!trans) { incAj = ldA; incAi = 1;}
                else       { incAj = 1;   incAi = ldA;}
                for(j=0; j<rows; j++) {
                    for(i=0; i<cols; i++) {
                        a[j*incAj+i*incAi] = x[i*incX] * a[j*incAj+i*incAi];
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
   A(i) := alpha * X(i) + A(i)

   Method Rindow\OpenBLAS\Math::
    public function add(
        int $trans,
        int $m,
        int $n,
        float $alpha,
        Buffer $X, int $offsetX, int $incX,
        Buffer $A, int $offsetA, int $ldA ) : void
 {{{ */
static PHP_METHOD(Math, add)
{
    php_interop_polite_math_matrix_linear_buffer_t* bufferX;
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    zend_bool trans;
    zend_long m;
    zend_long n;
    double alpha;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zval* a=NULL;
    zend_long offsetA;
    zend_long ldA;
    zend_long rows,cols;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 10, 10)
        Z_PARAM_BOOL(trans)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_DOUBLE(alpha)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_OBJECT(a) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "m", m)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    // Check Buffer X
    bufferX = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(bufferX,"x")) {
        return;
    }
    if(!trans) {
        rows = m; cols = n;
    } else {
        rows = n; cols = m;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", bufferX,cols,offsetX,incX)) {
        return;
    }

    // Check Buffer A
    bufferA = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_buffer_type(bufferA,"a")) {
        return;
    }
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        "A", bufferA,m,n,offsetA,ldA)) {
        return;
    }

    // Check Buffer X and A
    if(bufferX->dtype!=bufferA->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for X and A", 0);
        return;
    }

    switch (bufferX->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *x = &(((float *)bufferX->data)[offsetX]);
                float *a = &(((float *)bufferA->data)[offsetA]);
                zend_long j,incAj,incAi;
                if(!trans) { incAj = ldA; incAi = 1;}
                else       { incAj = 1;   incAi = ldA;}
                for(j=0; j<rows; j++) {
                    cblas_saxpy((blasint)cols,(float)alpha,
                        (float*)x,(blasint)incX,
                        (float*)(&a[j*incAj]),(blasint)incAi);
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *x = &(((double *)bufferX->data)[offsetX]);
                double *a = &(((double *)bufferA->data)[offsetA]);
                zend_long j,incAj,incAi;
                if(!trans) { incAj = ldA; incAi = 1;}
                else       { incAj = 1;   incAi = ldA;}
                for(j=0; j<rows; j++) {
                    cblas_daxpy((blasint)cols,(double)alpha,
                        (double*)x,(blasint)incX,
                        (double*)(&a[j*incAj]),(blasint)incAi);
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
   A(m,n) := X(n)

   Method Rindow\OpenBLAS\Math::
    public function duplicate(
        bool $trans,
        int $m,
        int $n,
        Buffer $X, int $offsetX, int $incX,
        Buffer $A, int $offsetA, int $ldA ) : void
 {{{ */
static PHP_METHOD(Math, duplicate)
{
    php_interop_polite_math_matrix_linear_buffer_t* bufferX;
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    zend_bool trans;
    zend_long m;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zval* a=NULL;
    zend_long offsetA;
    zend_long ldA;
    zend_long rows,cols;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 9, 9)
        Z_PARAM_BOOL(trans)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_OBJECT(a) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "m", m)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    if(!trans) {
        rows = m; cols = n;
    } else {
        rows = n; cols = m;
    }
    // Check Buffer X
    bufferX = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(bufferX,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", bufferX,cols,offsetX,incX)) {
        return;
    }

    // Check Buffer A
    bufferA = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_buffer_type(bufferA,"a")) {
        return;
    }
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        "A", bufferA,m,n,offsetA,ldA)) {
        return;
    }

    // Check Buffer X and Y
    if(bufferX->dtype!=bufferA->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for X and Y", 0);
        return;
    }

    switch (bufferX->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *x = &(((float *)bufferX->data)[offsetX]);
                float *a = &(((float *)bufferA->data)[offsetA]);
                zend_long j,incAj,incAi;
                if(!trans) { incAj = ldA; incAi = 1;}
                else       { incAj = 1;   incAi = ldA;}
                for(j=0; j<rows; j++) {
                    cblas_scopy((blasint)cols,
                        x, (blasint)incX,
                        &(a[j*incAj]), (blasint)incAi);
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *x = &(((double *)bufferX->data)[offsetX]);
                double *a = &(((double *)bufferA->data)[offsetA]);
                zend_long j,incAj,incAi;
                if(!trans) { incAj = ldA; incAi = 1;}
                else       { incAj = 1;   incAi = ldA;}
                for(j=0; j<rows; j++) {
                    cblas_dcopy((blasint)n,
                        x, (blasint)incX,
                        &(a[j*incAj]), (blasint)incAi);
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
   X := X ^ 2

   Method Rindow\OpenBLAS\Math::
    public function square(
        int $n,
        Buffer $X, int $offsetX, int $incX) : void
 {{{ */
static PHP_METHOD(Math, square)
{
    php_interop_polite_math_matrix_linear_buffer_t* buffer;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long i;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 4, 4)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    buffer = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(buffer,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    x[i*incX] = x[i*incX] * x[i*incX];
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *x = &(((double *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    x[i*incX] = x[i*incX] * x[i*incX];
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
   X := sqrt(X)

   Method Rindow\OpenBLAS\Math::
    public function sqrt(
        int $n,
        Buffer $X, int $offsetX, int $incX) : void
 {{{ */
static PHP_METHOD(Math, sqrt)
{
    php_interop_polite_math_matrix_linear_buffer_t* buffer;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long i;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 4, 4)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    buffer = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(buffer,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    float t;
                    t = x[i*incX];
                    // *** CAUTION ***
                    // disable checking for NaN and INFINITY values
                    //if(t<0.0) {
                    //    zend_throw_exception(spl_ce_RuntimeException, "Invalid value in sqrt.", 0);
                    //    return;
                    //}
                    x[i*incX] = sqrtf(t);
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *x = &(((double *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    double t;
                    t = x[i*incX];
                    // *** CAUTION ***
                    // disable checking for NaN and INFINITY values
                    //if(t<0.0) {
                    //    zend_throw_exception(spl_ce_RuntimeException, "Invalid value in sqrt.", 0);
                    //    return;
                    //}
                    x[i*incX] = sqrt(t);
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
   X := 1 / (a * sqrt(X) + b)

   Method Rindow\OpenBLAS\Math::
    public function rsqrt(
        int $n,
        float $alpha,
        Buffer $X, int $offsetX, int $incX,
        float $beta) : void
 {{{ */
static PHP_METHOD(Math, rsqrt)
{
    php_interop_polite_math_matrix_linear_buffer_t* buffer;
    zend_long n;
    double alpha;
    double beta;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long i;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 6, 6)
        Z_PARAM_LONG(n)
        Z_PARAM_DOUBLE(alpha)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_DOUBLE(beta)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    buffer = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(buffer,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    float t;
                    // *** CAUTION ***
                    // disable checking for NaN and INFINITY values
                    //if(x[i*incX]<0.0) {
                    //    zend_throw_exception(spl_ce_RuntimeException, "Invalid value in sqrt.", 0);
                    //    return;
                    //}
                    t = (float)alpha * sqrtf(x[i*incX]) + (float)beta;
                    // *** CAUTION ***
                    // disable checking for NaN and INFINITY values
                    //if(t==0.0) {
                    //    zend_throw_exception(spl_ce_RuntimeException, "Zero divide.", 0);
                    //    return;
                    //}
                    x[i*incX] = 1 / t;
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *x = &(((double *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    double t;
                    // *** CAUTION ***
                    // disable checking for NaN and INFINITY values
                    //if(x[i*incX]<0.0) {
                    //    zend_throw_exception(spl_ce_RuntimeException, "Invalid value in sqrt.", 0);
                    //    return;
                    //}
                    t = (double)alpha * sqrt(x[i*incX]) + (double)beta;
                    // *** CAUTION ***
                    // disable checking for NaN and INFINITY values
                    //if(t==0.0) {
                    //    zend_throw_exception(spl_ce_RuntimeException, "Zero divide.", 0);
                    //    return;
                    //}
                    x[i*incX] = 1 / t;
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
   X := X ^ a

   Method Rindow\OpenBLAS\Math::
    public function pow(
        int $n,
        Buffer $X, int $offsetX, int $incX,
        float $alpha) : void
 {{{ */
static PHP_METHOD(Math, pow)
{
    php_interop_polite_math_matrix_linear_buffer_t* buffer;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    double alpha;
    zend_long i;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 5, 5)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_DOUBLE(alpha)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    buffer = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(buffer,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    x[i*incX] = powf(x[i*incX], (float)alpha);
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *x = &(((double *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    x[i*incX] = pow(x[i*incX], (double)alpha);
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
   X(i) := e ^ X(i)

   Method Rindow\OpenBLAS\Math::
    public function exp(
        int $n,
        Buffer $X, int $offsetX, int $incX) : void
 {{{ */
static PHP_METHOD(Math, exp)
{
    php_interop_polite_math_matrix_linear_buffer_t* buffer;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long i;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 4, 4)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    buffer = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(buffer,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    x[i*incX] = expf(x[i*incX]);
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *x = &(((double *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    x[i*incX] = exp(x[i*incX]);
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
   X := log(X)

   Method Rindow\OpenBLAS\Math::
    public function log(
        int $n,
        Buffer $X, int $offsetX, int $incX) : void
 {{{ */
static PHP_METHOD(Math, log)
{
    php_interop_polite_math_matrix_linear_buffer_t* buffer;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long i;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 4, 4)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    buffer = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(buffer,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    float t;
                    t = x[i*incX];
                    // *** CAUTION ***
                    // disable checking for NaN and INFINITY values
                    //if(t<=0.0) {
                    //    zend_throw_exception(spl_ce_RuntimeException, "Invalid value in log.", 0);
                    //    return;
                    //}
                    x[i*incX] = logf(t);
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *x = &(((double *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    double t;
                    t = x[i*incX];
                    // *** CAUTION ***
                    // disable checking for NaN and INFINITY values
                    //if(t<=0.0) {
                    //    zend_throw_exception(spl_ce_RuntimeException, "Invalid value in log.", 0);
                    //    return;
                    //}
                    x[i*incX] = log(t);
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
   X := tanh(X)

   Method Rindow\OpenBLAS\Math::
    public function tanh(
        int $n,
        Buffer $X, int $offsetX, int $incX) : void
 {{{ */
static PHP_METHOD(Math, tanh)
{
    php_interop_polite_math_matrix_linear_buffer_t* buffer;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long i;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 4, 4)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    buffer = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(buffer,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    float t;
                    t = x[i*incX];
                    x[i*incX] = tanhf(t);
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *x = &(((double *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    double t;
                    t = x[i*incX];
                    x[i*incX] = tanh(t);
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
   X := sin(X)

   Method Rindow\OpenBLAS\Math::
    public function sin(
        int $n,
        Buffer $X, int $offsetX, int $incX) : void
 {{{ */
static PHP_METHOD(Math, sin)
{
    php_interop_polite_math_matrix_linear_buffer_t* buffer;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long i;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 4, 4)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    buffer = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(buffer,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    float t;
                    t = x[i*incX];
                    x[i*incX] = sinf(t);
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *x = &(((double *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    double t;
                    t = x[i*incX];
                    x[i*incX] = sin(t);
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
   X := cos(X)

   Method Rindow\OpenBLAS\Math::
    public function cos(
        int $n,
        Buffer $X, int $offsetX, int $incX) : void
 {{{ */
static PHP_METHOD(Math, cos)
{
    php_interop_polite_math_matrix_linear_buffer_t* buffer;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long i;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 4, 4)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    buffer = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(buffer,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    float t;
                    t = x[i*incX];
                    x[i*incX] = cosf(t);
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *x = &(((double *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    double t;
                    t = x[i*incX];
                    x[i*incX] = cos(t);
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
   X := tan(X)

   Method Rindow\OpenBLAS\Math::
    public function tan(
        int $n,
        Buffer $X, int $offsetX, int $incX) : void
 {{{ */
static PHP_METHOD(Math, tan)
{
    php_interop_polite_math_matrix_linear_buffer_t* buffer;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long i;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 4, 4)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    buffer = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(buffer,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    float t;
                    t = x[i*incX];
                    x[i*incX] = tanf(t);
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *x = &(((double *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    double t;
                    t = x[i*incX];
                    x[i*incX] = tan(t);
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
   X := 0

   Method Rindow\OpenBLAS\Math::
    public function zeros(
        int $n,
        Buffer $X, int $offsetX, int $incX) : void
 {{{ */
static PHP_METHOD(Math, zeros)
{
    php_interop_polite_math_matrix_linear_buffer_t* buffer;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 4, 4)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    buffer = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(buffer,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                zend_long i;
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    x[i*incX] = 0;
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                zend_long i;
                double *x = &(((double *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    x[i*incX] = 0;
                }
            }
            break;
        default:
            {
                zend_long i;
                int valueSize;
                uint8_t *x;
                valueSize = php_rindow_openblas_common_dtype_to_valuesize(buffer->dtype);
                x = php_rindow_openblas_get_address(buffer,offsetX,valueSize);
                if(incX==1) {
                    memset(x,0,valueSize*n);
                } else {
                    for(i=0;i<n;i++) {
                        memset(&x[i*incX],0,valueSize);
                    }
                }
            }
            break;
    }
}
/* }}} */

#define RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(data_type,get_type) { \
    data_type *x = &(((data_type *)buffer)[offset]); \
    *value = (get_type)(x[index*incWidth]);  \
    return 0; \
}

#define RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(data_type) { \
    data_type *x = &(((data_type *)buffer)[offset]); \
    x[index*incWidth] = (data_type)value;  \
    return 0; \
}

static int rindow_openblas_math_get_integer(
    zend_long dtype,void *buffer, zend_long offset,zend_long incWidth,
    zend_long index, zend_long *value)
{
    switch (dtype) {
        case php_interop_polite_math_matrix_dtype_bool:
            {
                uint8_t *x = &(((uint8_t *)buffer)[offset]);
                if(x[index*incWidth]==0) { *value=0; return 0; }
                else                     { *value=1; return 0; }
            }
        case php_interop_polite_math_matrix_dtype_int8:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(int8_t,zend_long);
        case php_interop_polite_math_matrix_dtype_uint8:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(uint8_t,zend_long);
        case php_interop_polite_math_matrix_dtype_int16:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(int16_t,zend_long);
        case php_interop_polite_math_matrix_dtype_uint16:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(uint16_t,zend_long);
        case php_interop_polite_math_matrix_dtype_int32:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(int32_t,zend_long);
        case php_interop_polite_math_matrix_dtype_uint32:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(uint32_t,zend_long);
        case php_interop_polite_math_matrix_dtype_int64:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(int64_t,zend_long);
        case php_interop_polite_math_matrix_dtype_uint64:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(uint64_t,zend_long);
        case php_interop_polite_math_matrix_dtype_float32:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(float,zend_long);
        case php_interop_polite_math_matrix_dtype_float64:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(double,zend_long);
        default:
            return -1;
    }
}

static int rindow_openblas_math_set_integer(
    zend_long dtype,void *buffer, zend_long offset,zend_long incWidth,
    zend_long index, zend_long value)
{
    switch (dtype) {
        case php_interop_polite_math_matrix_dtype_bool:
        {
            uint8_t *x = &(((uint8_t *)buffer)[offset]);
            if(value==0) { x[index*incWidth]=0; return 0; }
            else         { x[index*incWidth]=1; return 0; }
        }
        case php_interop_polite_math_matrix_dtype_int8:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(int8_t);
        case php_interop_polite_math_matrix_dtype_uint8:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(uint8_t);
        case php_interop_polite_math_matrix_dtype_int16:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(int16_t);
        case php_interop_polite_math_matrix_dtype_uint16:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(uint16_t);
        case php_interop_polite_math_matrix_dtype_int32:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(int32_t);
        case php_interop_polite_math_matrix_dtype_uint32:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(uint32_t);
        case php_interop_polite_math_matrix_dtype_int64:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(int64_t);
        case php_interop_polite_math_matrix_dtype_uint64:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(uint64_t);
        case php_interop_polite_math_matrix_dtype_float32:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(float);
        case php_interop_polite_math_matrix_dtype_float64:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(double);
        default:
            return -1;
    }
}

static int rindow_openblas_math_get_float(
    zend_long dtype,void *buffer, zend_long offset,zend_long incWidth,
    zend_long index, double *value)
{
    switch (dtype) {
        case php_interop_polite_math_matrix_dtype_bool:
            {
                uint8_t *x = &(((uint8_t *)buffer)[offset]);
                if(x[index*incWidth]==0) { *value=0; return 0; }
                else                     { *value=1; return 0; }
            }
        case php_interop_polite_math_matrix_dtype_int8:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(int8_t,double);
        case php_interop_polite_math_matrix_dtype_uint8:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(uint8_t,double);
        case php_interop_polite_math_matrix_dtype_int16:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(int16_t,double);
        case php_interop_polite_math_matrix_dtype_uint16:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(uint16_t,double);
        case php_interop_polite_math_matrix_dtype_int32:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(int32_t,double);
        case php_interop_polite_math_matrix_dtype_uint32:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(uint32_t,double);
        case php_interop_polite_math_matrix_dtype_int64:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(int64_t,double);
        case php_interop_polite_math_matrix_dtype_uint64:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(uint64_t,double);
        case php_interop_polite_math_matrix_dtype_float32:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(float,double);
        case php_interop_polite_math_matrix_dtype_float64:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(double,double);
        default:
            return -1;
    }
}

static int rindow_openblas_math_set_float(
    zend_long dtype,void *buffer, zend_long offset,zend_long incWidth,
    zend_long index, double value)
{
    switch (dtype) {
        case php_interop_polite_math_matrix_dtype_bool:
        {
            uint8_t *x = &(((uint8_t *)buffer)[offset]);
            if(value==0) { x[index*incWidth]=0; return 0; }
            else         { x[index*incWidth]=1; return 0; }
        }
        case php_interop_polite_math_matrix_dtype_int8:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(int8_t);
        case php_interop_polite_math_matrix_dtype_uint8:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(uint8_t);
        case php_interop_polite_math_matrix_dtype_int16:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(int16_t);
        case php_interop_polite_math_matrix_dtype_uint16:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(uint16_t);
        case php_interop_polite_math_matrix_dtype_int32:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(int32_t);
        case php_interop_polite_math_matrix_dtype_uint32:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(uint32_t);
        case php_interop_polite_math_matrix_dtype_int64:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(int64_t);
        case php_interop_polite_math_matrix_dtype_uint64:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(uint64_t);
        case php_interop_polite_math_matrix_dtype_float32:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(float);
        case php_interop_polite_math_matrix_dtype_float64:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(double);
        default:
            return -1;
    }
}


/*
   Y := a * onehot(X) + Y

   Method Rindow\OpenBLAS\Math::
    public function updateAddOnehot(
        int $m,
        int $n,
        float $a,
        Buffer $X, int $offsetX, int $incX,
        Buffer $Y, int $offsetY, int $ldY ) : void
 {{{ */
static PHP_METHOD(Math, updateAddOnehot)
{
    php_interop_polite_math_matrix_linear_buffer_t* bufferX;
    php_interop_polite_math_matrix_linear_buffer_t* bufferY;
    zend_long m;
    zend_long n;
    double a;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zval* y=NULL;
    zend_long offsetY;
    zend_long ldY;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 9, 9)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_DOUBLE(a)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_OBJECT(y) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetY)
        Z_PARAM_LONG(ldY)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "m", m)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    // Check Buffer X
    bufferX = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(bufferX,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", bufferX,m,offsetX,incX)) {
        return;
    }

    // Check Buffer Y
    bufferY = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(y);
    if(php_rindow_openblas_assert_buffer_type(bufferY,"y")) {
        return;
    }
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        "Y", bufferY,m,n,offsetY,ldY)) {
        return;
    }

    // Check Buffer A and Y
    if(bufferX->dtype==php_interop_polite_math_matrix_dtype_bool) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Data type of BufferX must not be bool", 0);
        return;
    }

    switch (bufferY->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
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
                    y[i*ldY+selector] += (float)a;
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *y = &(((double *)bufferY->data)[offsetY]);
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
                    y[i*ldY+selector] += (double)a;
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
   A := softmax(A)

   Method Rindow\OpenBLAS\Math::
    public function softmax(
        int $m,
        int $n,
        Buffer $A, int $offsetA, int $ldA) : void
 {{{ */
static PHP_METHOD(Math, softmax)
{
    php_interop_polite_math_matrix_linear_buffer_t* buffer;
    zend_long m;
    zend_long n;
    zval* a=NULL;
    zend_long offsetA;
    zend_long ldA;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 5, 5)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(a) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "m", m)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    buffer = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_buffer_type(buffer,"a")) {
        return;
    }
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        "A", buffer,m,n,offsetA,ldA)) {
        return;
    }
    switch (buffer->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *a = &(((float *)buffer->data)[offsetA]);
                zend_long i,j;
                for(i=0;i<m;i++,a+=ldA) {
                    float t,max_a,sum_exp;
                    max_a = s_max(n,a,1);
                    sum_exp = 0;
                    for(j=0;j<n;j++) {
                        t = expf(a[j]-max_a);
                        sum_exp += t;
                        a[j] = t;
                    }
                    // *** CAUTION ***
                    // disable checking for NaN and INFINITY values
                    //if(sum_exp==0.0) {
                    //    zend_throw_exception(spl_ce_RuntimeException, "Zero divide in softmax.", 0);
                    //    return;
                    //}
                    for(j=0;j<n;j++) {
                        a[j] = a[j] / sum_exp;
                    }
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *a = &(((double *)buffer->data)[offsetA]);
                zend_long i,j;
                for(i=0;i<m;i++,a+=ldA) {
                    double t,max_a,sum_exp;
                    max_a = d_max(n,a,1);
                    sum_exp = 0;
                    for(j=0;j<n;j++) {
                        t = exp(a[j]-max_a);
                        sum_exp += t;
                        a[j] = t;
                    }
                    // *** CAUTION ***
                    // disable checking for NaN and INFINITY values
                    //if(sum_exp==0.0) {
                    //    zend_throw_exception(spl_ce_RuntimeException, "Zero divide in softmax.", 0);
                    //    return;
                    //}
                    for(j=0;j<n;j++) {
                        a[j] = a[j] / sum_exp;
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
   Y(i) := 1  ( X(i) == Y(i) )
   Y(i) := 0  ( X(i) != Y(i) )

   Method Rindow\OpenBLAS\Math::
    public function equal(
        int $n,
        Buffer $X, int $offsetX, int $incX,
        Buffer $Y, int $offsetY, int $incY ) : void
 {{{ */
static PHP_METHOD(Math, equal)
{
    php_interop_polite_math_matrix_linear_buffer_t* bufferX;
    php_interop_polite_math_matrix_linear_buffer_t* bufferY;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zval* y=NULL;
    zend_long offsetY;
    zend_long incY;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 7, 7)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_OBJECT(y) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetY)
        Z_PARAM_LONG(incY)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    // Check Buffer X
    bufferX = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(bufferX,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", bufferX,n,offsetX,incX)) {
        return;
    }

    // Check Buffer Y
    bufferY = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(y);
    if(php_rindow_openblas_assert_buffer_type(bufferX,"y")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "Y", bufferY,n,offsetY,incY)) {
        return;
    }

    // Check Buffer X and Y
    if(bufferX->dtype!=bufferY->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for X and Y", 0);
        return;
    }

    switch (bufferX->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *x = &(((float *)bufferX->data)[offsetX]);
                float *y = &(((float *)bufferY->data)[offsetY]);
                zend_long i;
                for(i=0; i<n; i++) {
                    if(x[i*incX] == y[i*incY])
                        y[i*incY] = 1;
                    else
                        y[i*incY] = 0;
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *x = &(((double *)bufferX->data)[offsetX]);
                double *y = &(((double *)bufferY->data)[offsetY]);
                zend_long i;
                for(i=0; i<n; i++) {
                    if(x[i*incX] == y[i*incY])
                        y[i*incY] = 1;
                    else
                        y[i*incY] = 0;
                }
            }
            break;
        default:
            if(!php_rindow_openblas_common_dtype_is_int(bufferX->dtype)&&
                !php_rindow_openblas_common_dtype_is_bool(bufferX->dtype)) {
                zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
                return;
            }
            {
                int valueSize;
                uint8_t *x, *y;
                zend_long i,zero,one;
                zero = 0;
                one = 1;
                valueSize = php_rindow_openblas_common_dtype_to_valuesize(bufferX->dtype);
                x = php_rindow_openblas_get_address(bufferX,offsetX,valueSize);
                y = php_rindow_openblas_get_address(bufferY,offsetY,valueSize);
                for(i=0; i<n; i++) {
                    if(memcmp(&x[i*incX*valueSize],&y[i*incY*valueSize],valueSize)==0) {
                        memcpy(&y[i*incY*valueSize],&one,valueSize);
                    } else {
                        memcpy(&y[i*incY*valueSize],&zero,valueSize);
                    }
                }
            }
            break;
    }
}
/* }}} */

/*
   Y := cast<dtype> X

   Method Rindow\OpenBLAS\Math::
    public function astype(
        int $n,
        int $dtype,
        Buffer $X, int $offsetX, int $incX,
        Buffer $Y, int $offsetY, int $incY) : void
 {{{ */
static PHP_METHOD(Math, astype)
{
    php_interop_polite_math_matrix_linear_buffer_t* bufferX;
    php_interop_polite_math_matrix_linear_buffer_t* bufferY;
    zend_long n;
    zend_long dtype;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zval* y=NULL;
    zend_long offsetY;
    zend_long incY;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 8, 8)
        Z_PARAM_LONG(n)
        Z_PARAM_LONG(dtype)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_OBJECT(y) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetY)
        Z_PARAM_LONG(incY)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    bufferX = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(bufferX,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", bufferX,n,offsetX,incX)) {
        return;
    }
    bufferY = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(y);
    if(php_rindow_openblas_assert_buffer_type(bufferY,"y")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "Y", bufferY,n,offsetY,incY)) {
        return;
    }
    // Check dtype and Buffer Y
    if(dtype!=bufferY->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for Y", 0);
        return;
    }

    if(php_rindow_openblas_common_dtype_is_int(dtype) || php_rindow_openblas_common_dtype_is_bool(dtype)) {
        zend_long i,value;
        for(i=0;i<n;i++) {
            if(rindow_openblas_math_get_integer(
                        bufferX->dtype, bufferX->data, offsetX,incX,
                        i, &value)) {
                zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type of X.", 0);
                return;
            }
            rindow_openblas_math_set_integer(bufferY->dtype, bufferY->data, offsetY, incY, i, value);
        }
    } else if(php_rindow_openblas_common_dtype_is_float(dtype)) {
        zend_long i;
        double value;
        for(i=0;i<n;i++) {
            if(rindow_openblas_math_get_float(
                        bufferX->dtype, bufferX->data, offsetX,incX,
                        i, &value)) {
                zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type of X.", 0);
                return;
            }
            rindow_openblas_math_set_float(bufferY->dtype, bufferY->data, offsetY, incY, i, value);
        }
    } else {
        zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
        return;
    }
}
/* }}} */

/*
   B(m,n) :=  A(m,n) : trans=false
   B(n,m) :=  A(m,n) : trans=true

   Method Rindow\OpenBLAS\Math::
    public function matrixcopy(
        bool $trans,
        int $m,
        int $n,
        float $alpha,
        Buffer $A, int $offsetA, int $ldA,
        Buffer $B, int $offsetB, int $ldB ) : void
 {{{ */
static PHP_METHOD(Math, matrixcopy)
{
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    php_interop_polite_math_matrix_linear_buffer_t* bufferB;
    zend_bool trans;
    zend_long m;
    zend_long n;
    double alpha;
    zval* a=NULL;
    zend_long offsetA;
    zend_long ldA;
    zval* b=NULL;
    zend_long offsetB;
    zend_long ldB;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 10, 10)
        Z_PARAM_BOOL(trans)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_DOUBLE(alpha)
        Z_PARAM_OBJECT(a) // Interop\Polite\Math\Matrix\LinearBuffer

        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)
        Z_PARAM_OBJECT(b) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetB)
        Z_PARAM_LONG(ldB)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "m", m)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    // Check Buffer A
    bufferA = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_buffer_type(bufferA,"a")) {
        return;
    }
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        "A", bufferA,m,n,offsetA,ldA)) {
        return;
    }

    {
        zend_long rows,cols;
        if(!trans) {
            rows = m;
            cols = n;
        } else {
            rows = n;
            cols = m;
        }
        // Check Buffer B
        bufferB = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(b);
        if(php_rindow_openblas_assert_buffer_type(bufferB,"b")) {
            return;
        }
        if(php_rindow_openblas_assert_matrix_buffer_spec(
            "B", bufferB,rows,cols,offsetB,ldB)) {
            return;
        }
    }

    if(bufferA->dtype!=bufferB->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type A and B", 0);
        return;
    }

    switch (bufferA->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *a = &(((float *)bufferA->data)[offsetA]);
                float *b = &(((float *)bufferB->data)[offsetB]);
                zend_long i,j;
                if(!trans) {
                    for(i=0;i<m;i++) {
                        for(j=0;j<n;j++) {
                            b[i*ldB+j] = (float)alpha * a[i*ldA+j];
                        }
                    }
                } else {
                    for(i=0;i<m;i++) {
                        for(j=0;j<n;j++) {
                            b[j*ldB+i] = (float)alpha * a[i*ldA+j];
                        }
                    }
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *a = &(((double *)bufferA->data)[offsetA]);
                double *b = &(((double *)bufferB->data)[offsetB]);
                zend_long i,j;
                if(!trans) {
                    for(i=0;i<m;i++) {
                        for(j=0;j<n;j++) {
                            b[i*ldB+j] = (double)alpha * a[i*ldA+j];
                        }
                    }
                } else {
                    for(i=0;i<m;i++) {
                        for(j=0;j<n;j++) {
                            b[j*ldB+i] = (double)alpha * a[i*ldA+j];
                        }
                    }
                }
            }
            break;
        default:
            zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type of A.", 0);
            return;
    }
}
/* }}} */

/*
   Method Rindow\OpenBLAS\Math::
    public function imagecopy(
        int $height,
        int $width,
        int $channels,
        Buffer $A, int $offsetA,
        Buffer $B, int $offsetB,
        bool $channelsFirst,
        int $heightShift,
        int $widthShift,
        bool $verticalFlip,
        bool $horizontalFlip,
        bool $rgbFlip
    ) : void
 {{{ */
static PHP_METHOD(Math, imagecopy)
{
    zend_long height;
    zend_long width;
    zend_long channels;
    zval* a;
    zend_long offsetA;
    zval* b;
    zend_long offsetB;
    zend_bool channelsFirst;
    zend_long heightShift;
    zend_long widthShift;
    zend_bool verticalFlip;
    zend_bool horizontalFlip;
    zend_bool rgbFlip;
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    php_interop_polite_math_matrix_linear_buffer_t* bufferB;
    zend_long ldC,ldY,ldX;
    zend_long directionY,directionX,biasY,biasX;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 13, 13)
        Z_PARAM_LONG(height)
        Z_PARAM_LONG(width)
        Z_PARAM_LONG(channels)
        Z_PARAM_OBJECT(a) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetA)

        Z_PARAM_OBJECT(b) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetB)
        Z_PARAM_BOOL(channelsFirst)
        Z_PARAM_LONG(heightShift)
        Z_PARAM_LONG(widthShift)

        Z_PARAM_BOOL(verticalFlip)
        Z_PARAM_BOOL(horizontalFlip)
        Z_PARAM_BOOL(rgbFlip)
    ZEND_PARSE_PARAMETERS_END();

    if(height<1) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "height must be greater then 0", 0);
        return;
    }
    if(width<1) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "width must be greater then 0", 0);
        return;
    }
    if(channels<1) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "channels must be greater then 0", 0);
        return;
    }
    // Check Buffer A
    bufferA = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_buffer_type(bufferA,"a")) {
        return;
    }
    if(bufferA->size < height*width*channels+offsetA) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Matrix specification too large for bufferA", 0);
        return;
    }
    // Check Buffer B
    bufferB = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(b);
    if(php_rindow_openblas_assert_buffer_type(bufferA,"b")) {
        return;
    }
    if(bufferB->size < height*width*channels+offsetB) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Matrix specification too large for bufferB", 0);
        return;
    }

    if(bufferA->dtype!=bufferB->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type A and B", 0);
        return;
    }

    if(channelsFirst) {
        ldC = width*height;
        ldY = width;
        ldX = 1;
    } else {
        ldY = width*channels;
        ldX = channels;
        ldC = 1;
    }
    directionY = directionX = 1;
    biasY = biasX = 0;
    if(verticalFlip) {
        directionY = -directionY;
        biasY = height-1;
    }
    if(horizontalFlip) {
        directionX = -directionX;
        biasX = width-1;
    }
    biasY -= heightShift*directionY;
    biasX -= widthShift*directionX;

    switch (bufferA->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *a = &(((float *)bufferA->data)[offsetA]);
                float *b = &(((float *)bufferB->data)[offsetB]);
                for(zend_long y=0;y<height;y++) {
                    for(zend_long x=0;x<width;x++) {
                        for(zend_long c=0;c<channels;c++) {
                            zend_long sy = y*directionY+biasY;
                            zend_long sx = x*directionX+biasX;
                            if(sy<0) {
                                sy = 0;
                            } else if(sy>=height) {
                                sy = height-1;
                            }
                            if(sx<0) {
                                sx = 0;
                            } else if(sx>=width) {
                                sx = width-1;
                            }
                            zend_long srcc = (rgbFlip&&c<3)?(2-c):c;
                            b[y*ldY+x*ldX+c*ldC] =
                                a[sy*ldY+sx*ldX+srcc*ldC];
                        }
                    }
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *a = &(((double *)bufferA->data)[offsetA]);
                double *b = &(((double *)bufferB->data)[offsetB]);
                for(zend_long y=0;y<height;y++) {
                    for(zend_long x=0;x<width;x++) {
                        for(zend_long c=0;c<channels;c++) {
                            zend_long sy = y*directionY+biasY;
                            zend_long sx = x*directionX+biasX;
                            if(sy<0) {
                                sy = 0;
                            } else if(sy>=height) {
                                sy = height-1;
                            }
                            if(sx<0) {
                                sx = 0;
                            } else if(sx>=width) {
                                sx = width-1;
                            }
                            zend_long srcc = (rgbFlip&&c<3)?(2-c):c;
                            b[y*ldY+x*ldX+c*ldC] =
                                a[sy*ldY+sx*ldX+srcc*ldC];
                        }
                    }
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_uint8:
            {
                uint8_t *a = &(((uint8_t *)bufferA->data)[offsetA]);
                uint8_t *b = &(((uint8_t *)bufferB->data)[offsetB]);
                for(zend_long y=0;y<height;y++) {
                    for(zend_long x=0;x<width;x++) {
                        for(zend_long c=0;c<channels;c++) {
                            zend_long sy = y*directionY+biasY;
                            zend_long sx = x*directionX+biasX;
                            if(sy<0) {
                                sy = 0;
                            } else if(sy>=height) {
                                sy = height-1;
                            }
                            if(sx<0) {
                                sx = 0;
                            } else if(sx>=width) {
                                sx = width-1;
                            }
                            zend_long srcc = (rgbFlip&&c<3)?(2-c):c;
                            b[y*ldY+x*ldX+c*ldC] =
                                a[sy*ldY+sx*ldX+srcc*ldC];
                        }
                    }
                }
            }
            break;
        default:
            zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type of A.", 0);
            return;
    }
}
/* }}} */

/*
   X(n) :=  P

   Method Rindow\OpenBLAS\Math::
    public function fill(
        int $n,
        Buffer $value, int $offsetV,
        Buffer $X, int $offsetX, int $incX ) : void
 {{{ */
static PHP_METHOD(Math, fill)
{
    php_interop_polite_math_matrix_linear_buffer_t* bufferV;
    php_interop_polite_math_matrix_linear_buffer_t* bufferX;
    zend_long n;
    zval* value=NULL;
    zend_long offsetV;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 6, 6)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(value) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetV)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    // Check Buffer V
    bufferV = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(value);
    if(php_rindow_openblas_assert_buffer_type(bufferV,"value")) {
        return;
    }
    if(offsetV >= bufferV->size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "value buffer size is too small", 0);
        return;
    }
    // Check Buffer X
    bufferX = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(bufferX,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", bufferX,n,offsetX,incX)) {
        return;
    }

    if(bufferV->dtype!=bufferX->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type X and value", 0);
        return;
    }

    {
        size_t value_size = php_rindow_openblas_common_dtype_to_valuesize(bufferV->dtype);
        char *value = &(((char *)(bufferV->data))[offsetV*value_size]);
        char *x = &(((char *)(bufferX->data))[offsetX*value_size]);
        zend_long i;
        size_t step = incX*value_size;
        for(i=0;i<n;i++,x+=step) {
            memcpy(x,value,value_size);
        }
    }
}
/* }}} */

/*
      X := a  (X == NaN)
      X := X  (X != NaN)

   Method Rindow\OpenBLAS\Math::
    public function nan2num(
        int $n,
        Buffer $X, int $offsetX, int $incX,
        float $alpha) : void
 {{{ */
static PHP_METHOD(Math, nan2num)
{
    php_interop_polite_math_matrix_linear_buffer_t* buffer;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    double alpha;
    zend_long i;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 5, 5)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_DOUBLE(alpha)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    buffer = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(buffer,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    if(isnan(x[i*incX])) {
                        x[i*incX] = (float)alpha;
                    }
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *x = &(((double *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    if(isnan(x[i*incX])) {
                        x[i*incX] = (double)alpha;
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
   X := isnan(X)

   Method Rindow\OpenBLAS\Math::
    public function isnan(
        int $n,
        Buffer $X, int $offsetX, int $incX) : void
 {{{ */
static PHP_METHOD(Math, isnan)
{
    php_interop_polite_math_matrix_linear_buffer_t* buffer;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long i;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 4, 4)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    buffer = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(buffer,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    if(isnan(x[i*incX])) {
                        x[i*incX] = 1.0;
                    } else {
                        x[i*incX] = 0.0;
                    }
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *x = &(((double *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    if(isnan(x[i*incX])) {
                        x[i*incX] = 1.0;
                    } else {
                        x[i*incX] = 0.0;
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
   Y(n) := searchsorted( A(m), X(n) )

   Method Rindow\OpenBLAS\Math::
    public function searchsorted(
        int $m,
        Buffer $A, int $offsetA, int $incA,
        int $n,
        Buffer $X, int $offsetX, int $incX,
        bool $right,
        Buffer $Y, int $offsetY int $incY, ) : void
 {{{ */

static PHP_METHOD(Math, searchsorted)
{
    zend_long m;
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    zend_long offsetA;
    zend_long incA;
    zend_long n;
    php_interop_polite_math_matrix_linear_buffer_t* bufferX;
    zend_long offsetX;
    zend_long incX;
    zend_bool right;
    php_interop_polite_math_matrix_linear_buffer_t* bufferY;
    zend_long offsetY;
    zend_long incY;
    zval* a=NULL;
    zval* x=NULL;
    zval* y=NULL;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 12, 12)
        Z_PARAM_LONG(m)
        Z_PARAM_OBJECT(a) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(incA)
        Z_PARAM_LONG(n)

        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_BOOL(right)
        Z_PARAM_OBJECT(y) // Interop\Polite\Math\Matrix\LinearBuffer

        Z_PARAM_LONG(offsetY)
        Z_PARAM_LONG(incY)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "m", m)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }
    // Check Buffer A
    bufferA = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_buffer_type(bufferA,"A")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "A", bufferA,m,offsetA,incA)) {
        return;
    }

    // Check Buffer X
    bufferX = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(bufferX,"X")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", bufferX,n,offsetX,incX)) {
        return;
    }

    // Check Buffer Y
    bufferY = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(y);
    if(php_rindow_openblas_assert_buffer_type(bufferY,"Y")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "Y", bufferY,n,offsetY,incY)) {
        return;
    }

    // Check Buffer A and X
    if(bufferA->dtype!=bufferX->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for A and X", 0);
        return;
    }

    switch (bufferA->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *a = &(((float *)bufferA->data)[offsetA]);
                float *x = &(((float *)bufferX->data)[offsetX]);
                float value;
                zend_long idxA;
                zend_long idxX = 0;
                zend_long idxY = 0;
                for(zend_long i=0;i<n;i++,idxX+=incX,idxY+=incY) {
                    zend_long j;
                    value = x[idxX];
                    idxA = 0;
                    if(right) {
                        for(j=0;j<m;j++,idxA+=incA) {
                            if(!(value>=a[idxA])) {
                                break;
                            }
                        }
                    } else {
                        for(j=0;j<m;j++,idxA+=incA) {
                            if(!(value>a[idxA])) {
                                break;
                            }
                        }
                    }
                    rindow_openblas_math_set_integer(bufferY->dtype, bufferY->data, offsetY, 1, idxY, j);
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *a = &(((double *)bufferA->data)[offsetA]);
                double *x = &(((double *)bufferX->data)[offsetX]);
                double value;
                zend_long idxA;
                zend_long idxX = 0;
                zend_long idxY = 0;
                for(zend_long i=0;i<n;i++,idxX+=incX,idxY+=incY) {
                    zend_long j;
                    value = x[idxX];
                    idxA = offsetA;
                    if(right) {
                        for(j=0;j<m;j++,idxA+=incA) {
                            if(!(value>=a[idxA])) {
                                break;
                            }
                        }
                    } else {
                        for(j=0;j<m;j++,idxA+=incA) {
                            if(!(value>a[idxA])) {
                                break;
                            }
                        }
                    }
                    rindow_openblas_math_set_integer(bufferY->dtype, bufferY->data, offsetY, 1, idxY, j);
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
   Y(n) := X(n) + Y(n-1)

   Method Rindow\OpenBLAS\Math::
    public function cumsum(
        int $n,
        Buffer $X, int $offsetX, int $incX,
        bool $exclusive,
        bool $reverse,
        Buffer $Y, int $offsetY int $incY, ) : void
 {{{ */

static PHP_METHOD(Math, cumsum)
{
    zend_long n;
    php_interop_polite_math_matrix_linear_buffer_t* bufferX;
    zend_long offsetX;
    zend_long incX;
    zend_bool exclusive;
    zend_bool reverse;
    php_interop_polite_math_matrix_linear_buffer_t* bufferY;
    zend_long offsetY;
    zend_long incY;
    zval* x=NULL;
    zval* y=NULL;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 9, 9)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(x) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_BOOL(exclusive)

        Z_PARAM_BOOL(reverse)
        Z_PARAM_OBJECT(y) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetY)
        Z_PARAM_LONG(incY)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        "n", n)) {
        return;
    }

    // Check Buffer X
    bufferX = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(bufferX,"X")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", bufferX,n,offsetX,incX)) {
        return;
    }

    // Check Buffer Y
    bufferY = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(y);
    if(php_rindow_openblas_assert_buffer_type(bufferY,"Y")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "Y", bufferY,n,offsetY,incY)) {
        return;
    }

    // Check Buffer A and X
    if(bufferX->dtype!=bufferY->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for X and Y", 0);
        return;
    }

    switch (bufferX->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            {
                float *x = &(((float *)bufferX->data)[offsetX]);
                float *y = &(((float *)bufferY->data)[offsetY]);
                zend_long idxX,idxY;
                float value = 0.0;
                if(reverse) {
                    idxX = 0;
                    idxY = incY*(n-1);
                    incY = -incY;
                } else {
                    idxX = 0;
                    idxY = 0;
                }
                if(exclusive) {
                    for(zend_long i=0;i<n;i++,idxX+=incX,idxY+=incY) {
                        y[idxY] = value;
                        value += x[idxX];
                    }
                } else {
                    for(zend_long i=0;i<n;i++,idxX+=incX,idxY+=incY) {
                        value += x[idxX];
                        y[idxY] = value;
                    }
                }
            }
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            {
                double *x = &(((double *)bufferX->data)[offsetX]);
                double *y = &(((double *)bufferY->data)[offsetY]);
                zend_long idxX,idxY;
                double value = 0.0;
                if(reverse) {
                    idxX = 0;
                    idxY = incY*(n-1);
                    incY = -incY;
                } else {
                    idxX = 0;
                    idxY = 0;
                }
                if(exclusive) {
                    for(zend_long i=0;i<n;i++,idxX+=incX,idxY+=incY) {
                        y[idxY] = value;
                        value += x[idxX];
                    }
                } else {
                    for(zend_long i=0;i<n;i++,idxX+=incX,idxY+=incY) {
                        value += x[idxX];
                        y[idxY] = value;
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

#include "Math_gather.c"
#include "Math_repeat.c"
#include "Math_slice.c"
#include "Math_reduction.c"
#include "Math_im2col1d.c"
#include "Math_im2col2d.c"
#include "Math_im2col3d.c"
#include "Math_random.c"

ZEND_BEGIN_ARG_INFO_EX(ai_Math_sum, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_imax, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_imin, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_increment, 0, 0, 6)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, alpha)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_INFO(0, beta)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_reciprocal, 0, 0, 6)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, alpha)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_INFO(0, beta)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_maximum, 0, 0, 8)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_INFO(0, ldA)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_minimum, 0, 0, 8)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_INFO(0, ldA)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_greater, 0, 0, 8)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_INFO(0, ldA)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_greaterEqual, 0, 0, 8)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_INFO(0, ldA)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_less, 0, 0, 8)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_INFO(0, ldA)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_lessEqual, 0, 0, 8)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_INFO(0, ldA)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_multiply, 0, 0, 9)
    ZEND_ARG_INFO(0, trans)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_OBJ_INFO(0, y, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, incY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_add, 0, 0, 10)
    ZEND_ARG_INFO(0, trans)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, alpha)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_OBJ_INFO(0, y, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, incY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_duplicate, 0, 0, 9)
    ZEND_ARG_INFO(0, trans)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_OBJ_INFO(0, y, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, incY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_square, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_sqrt, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_rsqrt, 0, 0, 6)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, alpha)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_INFO(0, beta)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_pow, 0, 0, 5)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_INFO(0, alpha)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_exp, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_log, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_tanh, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_sin, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_cos, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_tan, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_zeros, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_gather, 0, 0, 11)
    ZEND_ARG_INFO(0, reverse)
    ZEND_ARG_INFO(0, addMode)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, k)
    ZEND_ARG_INFO(0, numClass)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_OBJ_INFO(0, a, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_OBJ_INFO(0, b, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetB)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_reduceGather, 0, 0, 11)
    ZEND_ARG_INFO(0, reverse)
    ZEND_ARG_INFO(0, addMode)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, numClass)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_OBJ_INFO(0, a, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_OBJ_INFO(0, b, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetB)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_repeat, 0, 0, 7)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, k)
    ZEND_ARG_INFO(0, repeats)
    ZEND_ARG_OBJ_INFO(0, a, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_OBJ_INFO(0, b, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetB)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_slice, 0, 0, 18)
    ZEND_ARG_INFO(0, reverse)
    ZEND_ARG_INFO(0, addMode)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, k)
    ZEND_ARG_INFO(0, size)
    ZEND_ARG_OBJ_INFO(0, a, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_INFO(0, incA)
    ZEND_ARG_OBJ_INFO(0, y, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, incY)
    ZEND_ARG_INFO(0, startAxis0)
    ZEND_ARG_INFO(0, sizeAxis0)
    ZEND_ARG_INFO(0, startAxis1)
    ZEND_ARG_INFO(0, sizeAxis1)
    ZEND_ARG_INFO(0, startAxis2)
    ZEND_ARG_INFO(0, sizeAxis2)
ZEND_END_ARG_INFO()


ZEND_BEGIN_ARG_INFO_EX(ai_Math_updateAddOnehot, 0, 0, 9)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, a)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_OBJ_INFO(0, y, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, ldY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_softmax, 0, 0, 5)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_equal, 0, 0, 7)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_OBJ_INFO(0, y, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, incY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_reduceSum, 0, 0, 7)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, k)
    ZEND_ARG_OBJ_INFO(0, a, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_OBJ_INFO(0, b, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetB)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_reduceMax, 0, 0, 7)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, k)
    ZEND_ARG_OBJ_INFO(0, a, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_OBJ_INFO(0, b, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetB)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_reduceArgMax, 0, 0, 7)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, k)
    ZEND_ARG_OBJ_INFO(0, a, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_OBJ_INFO(0, b, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetB)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_astype, 0, 0, 8)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, dtype)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_OBJ_INFO(0, y, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, incY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_matrixcopy, 0, 0, 10)
    ZEND_ARG_INFO(0, trans)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, alpha)
    ZEND_ARG_OBJ_INFO(0, a, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)

    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_INFO(0, ldA)
    ZEND_ARG_OBJ_INFO(0, b, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetB)
    ZEND_ARG_INFO(0, ldB)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_imagecopy, 0, 0, 13)
    ZEND_ARG_INFO(0, height)
    ZEND_ARG_INFO(0, width)
    ZEND_ARG_INFO(0, channels)
    ZEND_ARG_OBJ_INFO(0, a, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_OBJ_INFO(0, b, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetB)
    ZEND_ARG_INFO(0, channelsFirst)
    ZEND_ARG_INFO(0, heightShift)
    ZEND_ARG_INFO(0, widthShift)
    ZEND_ARG_INFO(0, verticalFlip)
    ZEND_ARG_INFO(0, horizontalFlip)
    ZEND_ARG_INFO(0, rgbFlip)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_fill, 0, 0, 6)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, value, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetV)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_nan2num, 0, 0, 5)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_INFO(0, alpha)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_isnan, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_searchsorted, 0, 0, 12)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_OBJ_INFO(0, a, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_INFO(0, incA)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_INFO(0, right)
    ZEND_ARG_OBJ_INFO(0, y, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, incY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_cumsum, 0, 0, 9)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_INFO(0, exclusive)
    ZEND_ARG_INFO(0, reverse)
    ZEND_ARG_OBJ_INFO(0, y, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, incY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_im2col1d, 0, 0, 16)
    ZEND_ARG_INFO(0, reverse)
    ZEND_ARG_OBJ_INFO(0, images_obj, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, images_offset)
    ZEND_ARG_INFO(0, images_size)
    ZEND_ARG_INFO(0, batches)

    ZEND_ARG_INFO(0, im_w)
    ZEND_ARG_INFO(0, channels)
    ZEND_ARG_INFO(0, filter_w)
    ZEND_ARG_INFO(0, stride_w)
    ZEND_ARG_INFO(0, padding)

    ZEND_ARG_INFO(0, channels_first)
    ZEND_ARG_INFO(0, dilation_w)
    ZEND_ARG_INFO(0, cols_channels_first)
    ZEND_ARG_OBJ_INFO(0, cols_obj,Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, cols_offset)

    ZEND_ARG_INFO(0, cols_size)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_im2col2d, 0, 0, 18)
    ZEND_ARG_INFO(0, reverse)
    ZEND_ARG_OBJ_INFO(0, images_obj, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, images_offset)
    ZEND_ARG_INFO(0, images_size)
    ZEND_ARG_INFO(0, batches)

    ZEND_ARG_INFO(0, im_h)
    ZEND_ARG_INFO(0, im_w)
    ZEND_ARG_INFO(0, channels)
    ZEND_ARG_INFO(0, filter_h)
    ZEND_ARG_INFO(0, filter_w)

    ZEND_ARG_INFO(0, stride_h)
    ZEND_ARG_INFO(0, stride_w)
    ZEND_ARG_INFO(0, padding)
    ZEND_ARG_INFO(0, channels_first)
    ZEND_ARG_INFO(0, dilation_h)

    ZEND_ARG_INFO(0, dilation_w)
    ZEND_ARG_INFO(0, cols_channels_first)
    ZEND_ARG_OBJ_INFO(0, cols_obj,Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, cols_offset)
    ZEND_ARG_INFO(0, cols_size)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_im2col3d, 0, 0, 24)
    ZEND_ARG_INFO(0, reverse)
    ZEND_ARG_OBJ_INFO(0, images_obj, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, images_offset)
    ZEND_ARG_INFO(0, images_size)
    ZEND_ARG_INFO(0, batches)

    ZEND_ARG_INFO(0, im_d)
    ZEND_ARG_INFO(0, im_h)
    ZEND_ARG_INFO(0, im_w)
    ZEND_ARG_INFO(0, channels)
    ZEND_ARG_INFO(0, filter_d)

    ZEND_ARG_INFO(0, filter_h)
    ZEND_ARG_INFO(0, filter_w)
    ZEND_ARG_INFO(0, stride_d)
    ZEND_ARG_INFO(0, stride_h)
    ZEND_ARG_INFO(0, stride_w)

    ZEND_ARG_INFO(0, padding)
    ZEND_ARG_INFO(0, channels_first)
    ZEND_ARG_INFO(0, dilation_d)
    ZEND_ARG_INFO(0, dilation_h)
    ZEND_ARG_INFO(0, dilation_w)

    ZEND_ARG_INFO(0, cols_channels_first)
    ZEND_ARG_OBJ_INFO(0, cols_obj,Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, cols_offset)
    ZEND_ARG_INFO(0, cols_size)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_randomUniform, 0, 0, 7)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_INFO(0, low)
    ZEND_ARG_INFO(0, high)
    ZEND_ARG_INFO(0, seed)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_randomNormal, 0, 0, 7)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_INFO(0, mean)
    ZEND_ARG_INFO(0, scale)
    ZEND_ARG_INFO(0, seed)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_randomSequence, 0, 0, 6)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, size)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_INFO(0, seed)
ZEND_END_ARG_INFO()

/* {{{ Rindow\OpenBLAS\Blas function entries */
static zend_function_entry php_rindow_openblas_math_me[] = {
    /* clang-format off */
    PHP_ME(Math, sum,            ai_Math_sum,            ZEND_ACC_PUBLIC)
    PHP_ME(Math, imax,           ai_Math_imax,           ZEND_ACC_PUBLIC)
    PHP_ME(Math, imin,           ai_Math_imin,           ZEND_ACC_PUBLIC)
    PHP_ME(Math, increment,      ai_Math_increment,      ZEND_ACC_PUBLIC)
    PHP_ME(Math, reciprocal,     ai_Math_reciprocal,     ZEND_ACC_PUBLIC)
    PHP_ME(Math, maximum,        ai_Math_maximum,        ZEND_ACC_PUBLIC)
    PHP_ME(Math, minimum,        ai_Math_minimum,        ZEND_ACC_PUBLIC)
    PHP_ME(Math, greater,        ai_Math_greater,        ZEND_ACC_PUBLIC)
    PHP_ME(Math, greaterEqual,   ai_Math_greaterEqual,   ZEND_ACC_PUBLIC)
    PHP_ME(Math, less,           ai_Math_less,           ZEND_ACC_PUBLIC)
    PHP_ME(Math, lessEqual,      ai_Math_lessEqual,      ZEND_ACC_PUBLIC)
    PHP_ME(Math, multiply,       ai_Math_multiply,       ZEND_ACC_PUBLIC)
    PHP_ME(Math, add,            ai_Math_add,            ZEND_ACC_PUBLIC)
    PHP_ME(Math, duplicate,      ai_Math_duplicate,      ZEND_ACC_PUBLIC)
    PHP_ME(Math, square,         ai_Math_square,         ZEND_ACC_PUBLIC)
    PHP_ME(Math, sqrt,           ai_Math_sqrt,           ZEND_ACC_PUBLIC)
    PHP_ME(Math, rsqrt,          ai_Math_rsqrt,          ZEND_ACC_PUBLIC)
    PHP_ME(Math, pow,            ai_Math_pow,            ZEND_ACC_PUBLIC)
    PHP_ME(Math, exp,            ai_Math_exp,            ZEND_ACC_PUBLIC)
    PHP_ME(Math, log,            ai_Math_log,            ZEND_ACC_PUBLIC)
    PHP_ME(Math, tanh,           ai_Math_tanh,           ZEND_ACC_PUBLIC)
    PHP_ME(Math, sin,            ai_Math_sin,            ZEND_ACC_PUBLIC)
    PHP_ME(Math, cos,            ai_Math_cos,            ZEND_ACC_PUBLIC)
    PHP_ME(Math, tan,            ai_Math_tan,            ZEND_ACC_PUBLIC)
    PHP_ME(Math, zeros,          ai_Math_zeros,          ZEND_ACC_PUBLIC)
    PHP_ME(Math, gather,         ai_Math_gather,         ZEND_ACC_PUBLIC)
    PHP_ME(Math, reduceGather,   ai_Math_reduceGather,   ZEND_ACC_PUBLIC)
    PHP_ME(Math, repeat,         ai_Math_repeat,         ZEND_ACC_PUBLIC)
    PHP_ME(Math, slice,          ai_Math_slice,          ZEND_ACC_PUBLIC)
    PHP_ME(Math, updateAddOnehot,ai_Math_updateAddOnehot,ZEND_ACC_PUBLIC)
    PHP_ME(Math, softmax,        ai_Math_softmax,        ZEND_ACC_PUBLIC)
    PHP_ME(Math, equal,          ai_Math_equal,          ZEND_ACC_PUBLIC)
    PHP_ME(Math, astype,         ai_Math_astype,         ZEND_ACC_PUBLIC)
    PHP_ME(Math, matrixcopy,     ai_Math_matrixcopy,     ZEND_ACC_PUBLIC)
    PHP_ME(Math, imagecopy,      ai_Math_imagecopy,      ZEND_ACC_PUBLIC)
    PHP_ME(Math, fill,           ai_Math_fill,           ZEND_ACC_PUBLIC)
    PHP_ME(Math, nan2num,        ai_Math_nan2num,        ZEND_ACC_PUBLIC)
    PHP_ME(Math, isnan,          ai_Math_isnan,          ZEND_ACC_PUBLIC)
    PHP_ME(Math, searchsorted,   ai_Math_searchsorted,   ZEND_ACC_PUBLIC)
    PHP_ME(Math, cumsum,         ai_Math_cumsum,         ZEND_ACC_PUBLIC)
    PHP_ME(Math, reduceSum,      ai_Math_reduceSum,      ZEND_ACC_PUBLIC)
    PHP_ME(Math, reduceMax,      ai_Math_reduceMax,      ZEND_ACC_PUBLIC)
    PHP_ME(Math, reduceArgMax,   ai_Math_reduceArgMax,   ZEND_ACC_PUBLIC)
    PHP_ME(Math, im2col1d,       ai_Math_im2col1d,       ZEND_ACC_PUBLIC)
    PHP_ME(Math, im2col2d,       ai_Math_im2col2d,       ZEND_ACC_PUBLIC)
    PHP_ME(Math, im2col3d,       ai_Math_im2col3d,       ZEND_ACC_PUBLIC)
    PHP_ME(Math, randomUniform,  ai_Math_randomUniform,  ZEND_ACC_PUBLIC)
    PHP_ME(Math, randomNormal,   ai_Math_randomNormal,   ZEND_ACC_PUBLIC)
    PHP_ME(Math, randomSequence, ai_Math_randomSequence, ZEND_ACC_PUBLIC)
    PHP_FE_END
    /* clang-format on */
};
/* }}} */

/* Class Rindow\OpenBLAS\Math {{{ */
static zend_class_entry* rindow_openblas_math_ce;

void php_rindow_openblas_math_init_ce(INIT_FUNC_ARGS)
{
    zend_class_entry ce;

    INIT_NS_CLASS_ENTRY(ce, "Rindow\\OpenBLAS", "Math", php_rindow_openblas_math_me);
    rindow_openblas_math_ce = zend_register_internal_class(&ce);
    rindow_openblas_math_ce->create_object = php_rindow_openblas_math_create_object;

    memcpy(&rindow_openblas_math_object_handlers, zend_get_std_object_handlers(), sizeof(zend_object_handlers));
    rindow_openblas_math_object_handlers.offset    = 0;
    rindow_openblas_math_object_handlers.free_obj  = php_rindow_openblas_math_free_object;
    rindow_openblas_math_object_handlers.clone_obj = NULL;

    //zend_class_implements(rindow_openblas_math_ce, 2, spl_ce_ArrayAccess, spl_ce_Countable);
}
/* }}} */
