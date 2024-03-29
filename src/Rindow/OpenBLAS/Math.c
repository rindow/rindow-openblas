#include <php.h>
#include <Zend/zend_interfaces.h>
#include <Zend/zend_exceptions.h>
#include <ext/spl/spl_iterators.h>
#include <ext/spl/spl_exceptions.h>
#include <ext/standard/php_rand.h>
#include <cblas.h>
#include <stdint.h>
#include <math.h>
#include <Interop/Polite/Math/Matrix.h>


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "php_rindow_openblas.h"

#define RINDOW_MATLIB_INCLUDING_SOURCE 1

#include "../../../matlib.c"

//static float s_max(long n, float *x, long incX)
//{
//    long i;
//    float a;
//    a = x[0];
//    for(i=1;i<n;i++) {
//        // if NaN set NaN
//        // Compatible with reduce_max of tensorflow 2.6
//        if(!(a>=x[i*incX])) {
//            a = x[i*incX];
//        }
//    }
//    return a;
//}
//
//static double d_max(long n, double *x, long incX)
//{
//    long i;
//    double a;
//    a = x[0];
//    for(i=1;i<n;i++) {
//        if(a<x[i*incX]) {
//            a = x[i*incX];
//        }
//    }
//    return a;
//}
//
//static long s_argmax(long n, float *x, long incX)
//{
//    long i;
//    long idx;
//    float a;
//    idx = 0;
//    a = x[0];
//    for(i=1;i<n;i++) {
//        if(a<x[i*incX]) {
//            idx = i;
//            a = x[i*incX];
//        }
//    }
//    return idx;
//}
//
//static long d_argmax(long n, double *x, long incX)
//{
//    long i;
//    long idx;
//    double a;
//    idx = 0;
//    a = x[0];
//    for(i=1;i<n;i++) {
//        if(a<x[i*incX]) {
//            idx = i;
//            a = x[i*incX];
//        }
//    }
//    return idx;
//}
//
//static float s_sum_sb(long n, float *x, long incX)
//{
//    long i;
//    float a=0;
//    for(i=0; i<n; i++) {
//        a += x[i*incX];
//    }
//    return a;
//}
//
//static double d_sum_sb(long n, double *x, long incX)
//{
//    long i;
//    double a=0;
//    for(i=0; i<n; i++) {
//        a += x[i*incX];
//    }
//    return a;
//}
//

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

#define PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(data_type,variable,buffer,offset) \
    data_type  *variable; \
    variable = &(((data_type *)buffer->data)[offset]); \


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

//#define PHP_RINDOW_OPENBLAS_MATH_ADD_TEMPLATE(data_type) { \
//    data_type  *pDataX; \
//    data_type  *pDataY; \
//    pDataX = (data_type *)values; \
//    pDataY = (data_type *)target; \
//    for (i=0; i<n; i++) { \
//        *pDataY += *pDataX; \
//        pDataX+=incValue; \
//        pDataY+=incTarget; \
//    } \
//}
//int php_rindow_openblas_math_add(
//    zend_long n,
//    zend_long dtype,
//    void* values,
//    zend_long incValue,
//    void* target,
//    zend_long incTarget
//    )
//{
//    switch (dtype) {
//        zend_long i;
//        case php_interop_polite_math_matrix_dtype_float32:
//            PHP_RINDOW_OPENBLAS_MATH_ADD_TEMPLATE(float)
//            break;
//        case php_interop_polite_math_matrix_dtype_float64:
//            PHP_RINDOW_OPENBLAS_MATH_ADD_TEMPLATE(double)
//            break;
//        case php_interop_polite_math_matrix_dtype_int8:
//            PHP_RINDOW_OPENBLAS_MATH_ADD_TEMPLATE(int8_t)
//            break;
//        case php_interop_polite_math_matrix_dtype_uint8:
//            PHP_RINDOW_OPENBLAS_MATH_ADD_TEMPLATE(uint8_t)
//            break;
//        case php_interop_polite_math_matrix_dtype_int16:
//            PHP_RINDOW_OPENBLAS_MATH_ADD_TEMPLATE(int16_t)
//            break;
//        case php_interop_polite_math_matrix_dtype_uint16:
//            PHP_RINDOW_OPENBLAS_MATH_ADD_TEMPLATE(uint16_t)
//            break;
//        case php_interop_polite_math_matrix_dtype_int32:
//            PHP_RINDOW_OPENBLAS_MATH_ADD_TEMPLATE(int32_t)
//            break;
//        case php_interop_polite_math_matrix_dtype_uint32:
//            PHP_RINDOW_OPENBLAS_MATH_ADD_TEMPLATE(uint32_t)
//            break;
//        case php_interop_polite_math_matrix_dtype_int64:
//            PHP_RINDOW_OPENBLAS_MATH_ADD_TEMPLATE(int64_t)
//            break;
//        case php_interop_polite_math_matrix_dtype_uint64:
//            PHP_RINDOW_OPENBLAS_MATH_ADD_TEMPLATE(uint64_t)
//            break;
//        default:
//            zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
//            return -1;
//    }
//    return 0;
//}

//#define PHP_RINDOW_OPENBLAS_MATH_COPY_TEMPLATE(data_type) { \
//    data_type  *pDataX; \
//    data_type  *pDataY; \
//    pDataX = (data_type *)source; \
//    pDataY = (data_type *)dest; \
//    for (i=0; i<n; i++) { \
//        *pDataY += *pDataX; \
//        pDataX+=incSource; \
//        pDataY+=incDest; \
//    } \
//}
//int php_rindow_openblas_math_copy(
//    zend_long n,
//    zend_long dtype,
//    void* source,
//    zend_long incSource,
//    void* dest,
//    zend_long incDest
//    )
//{
//    switch (dtype) {
//        zend_long i;
//        case php_interop_polite_math_matrix_dtype_float32:
//            PHP_RINDOW_OPENBLAS_MATH_COPY_TEMPLATE(float)
//            break;
//        case php_interop_polite_math_matrix_dtype_float64:
//            PHP_RINDOW_OPENBLAS_MATH_COPY_TEMPLATE(double)
//            break;
//        case php_interop_polite_math_matrix_dtype_int8:
//            PHP_RINDOW_OPENBLAS_MATH_COPY_TEMPLATE(int8_t)
//            break;
//        case php_interop_polite_math_matrix_dtype_uint8:
//            PHP_RINDOW_OPENBLAS_MATH_COPY_TEMPLATE(uint8_t)
//            break;
//        case php_interop_polite_math_matrix_dtype_int16:
//            PHP_RINDOW_OPENBLAS_MATH_COPY_TEMPLATE(int16_t)
//            break;
//        case php_interop_polite_math_matrix_dtype_uint16:
//            PHP_RINDOW_OPENBLAS_MATH_COPY_TEMPLATE(uint16_t)
//            break;
//        case php_interop_polite_math_matrix_dtype_int32:
//            PHP_RINDOW_OPENBLAS_MATH_COPY_TEMPLATE(int32_t)
//            break;
//        case php_interop_polite_math_matrix_dtype_uint32:
//            PHP_RINDOW_OPENBLAS_MATH_COPY_TEMPLATE(uint32_t)
//            break;
//        case php_interop_polite_math_matrix_dtype_int64:
//            PHP_RINDOW_OPENBLAS_MATH_COPY_TEMPLATE(int64_t)
//            break;
//        case php_interop_polite_math_matrix_dtype_uint64:
//            PHP_RINDOW_OPENBLAS_MATH_COPY_TEMPLATE(uint64_t)
//            break;
//        default:
//            zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
//            return -1;
//    }
//    return 0;
//}

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
        case php_interop_polite_math_matrix_dtype_float32:{
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,buffer,offsetX)
            result = rindow_matlib_s_sum(n,pDataX,incX);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64:{
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,buffer,offsetX)
            result = rindow_matlib_d_sum(n,pDataX,incX);
            break;
        }
        case php_interop_polite_math_matrix_dtype_int8:
        case php_interop_polite_math_matrix_dtype_uint8:
        case php_interop_polite_math_matrix_dtype_int16:
        case php_interop_polite_math_matrix_dtype_uint16:
        case php_interop_polite_math_matrix_dtype_int32:
        case php_interop_polite_math_matrix_dtype_uint32:
        case php_interop_polite_math_matrix_dtype_int64:
        case php_interop_polite_math_matrix_dtype_uint64:
        case php_interop_polite_math_matrix_dtype_bool: {
            void *pDataX = rindow_matlib_common_get_address(buffer->dtype,buffer->data,offsetX);
            result = (double)rindow_matlib_i_sum(buffer->dtype, n, pDataX, incX);
            break;
        }
        default:{
            zend_throw_exception_ex(spl_ce_InvalidArgumentException, 0, "Unsupported data type.");
            return;
        }
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
        case php_interop_polite_math_matrix_dtype_float32:{
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,buffer,offsetX)
            resultIdx = rindow_matlib_s_imax(n,pDataX,incX);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64:{
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,buffer,offsetX)
            resultIdx = rindow_matlib_d_imax(n,pDataX,incX);
            break;
        }
        default: {
            if(!rindow_matlib_common_dtype_is_int(buffer->dtype)) {
                zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
                return;
            }
            void *pDataX = rindow_matlib_common_get_address(buffer->dtype,buffer->data,offsetX);
            resultIdx = rindow_matlib_i_imax(buffer->dtype, n, pDataX, incX);
            break;
        }
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
        case php_interop_polite_math_matrix_dtype_float32:{
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,buffer,offsetX)
            resultIdx = rindow_matlib_s_imin(n,pDataX,incX);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64:{
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,buffer,offsetX)
            resultIdx = rindow_matlib_d_imin(n,pDataX,incX);
            break;
        }
        default: {
            if(!rindow_matlib_common_dtype_is_int(buffer->dtype)) {
                zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
                return;
            }
            void *pDataX = rindow_matlib_common_get_address(buffer->dtype,buffer->data,offsetX);
            resultIdx = rindow_matlib_i_imin(buffer->dtype, n, pDataX, incX);
            break;
        }
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
        case php_interop_polite_math_matrix_dtype_float32:{
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,buffer,offsetX)
            rindow_matlib_s_increment(n, pDataX, incX, alpha, beta);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64:{
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,buffer,offsetX)
            rindow_matlib_d_increment(n, pDataX, incX, alpha, beta);
            break;
        }
        default:{
            zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
            return;
        }
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,buffer,offsetX)
            rindow_matlib_s_reciprocal(n, pDataX, incX, alpha, beta);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,buffer,offsetX)
            rindow_matlib_d_reciprocal(n, pDataX, incX, alpha, beta);
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,bufferX,offsetX)
            rindow_matlib_s_maximum(m, n, pDataA, ldA, pDataX, incX);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,bufferX,offsetX)
            rindow_matlib_d_maximum(m, n, pDataA, ldA, pDataX, incX);
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,bufferX,offsetX)
            rindow_matlib_s_minimum(m, n, pDataA, ldA, pDataX, incX);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,bufferX,offsetX)
            rindow_matlib_d_minimum(m, n, pDataA, ldA, pDataX, incX);
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,bufferX,offsetX)
            rindow_matlib_s_greater(m, n, pDataA, ldA, pDataX, incX);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,bufferX,offsetX)
            rindow_matlib_d_greater(m, n, pDataA, ldA, pDataX, incX);
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,bufferX,offsetX)
            rindow_matlib_s_greater_equal(m, n, pDataA, ldA, pDataX, incX);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,bufferX,offsetX)
            rindow_matlib_d_greater_equal(m, n, pDataA, ldA, pDataX, incX);
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,bufferX,offsetX)
            rindow_matlib_s_less(m, n, pDataA, ldA, pDataX, incX);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,bufferX,offsetX)
            rindow_matlib_d_less(m, n, pDataA, ldA, pDataX, incX);
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,bufferX,offsetX)
            rindow_matlib_s_less_equal(m, n, pDataA, ldA, pDataX, incX);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,bufferX,offsetX)
            rindow_matlib_d_less_equal(m, n, pDataA, ldA, pDataX, incX);
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
    zend_long cols;
    zend_long transCode;

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

    if(trans) {
        transCode = RINDOW_MATLIB_TRANS;
        cols = m;
    } else {
        transCode = RINDOW_MATLIB_NO_TRANS;
        cols = n;
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

    // Check Buffer X and A
    if(bufferX->dtype!=bufferA->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for X and A", 0);
        return;
    }

    switch (bufferX->dtype) {
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,bufferX,offsetX)
            rindow_matlib_s_multiply(transCode, m, n, pDataX, incX, pDataA, ldA);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,bufferX,offsetX)
            rindow_matlib_d_multiply(transCode, m, n, pDataX, incX, pDataA, ldA);
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
    zend_long cols;
    zend_long transCode;

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
    if(trans) {
        transCode = RINDOW_MATLIB_TRANS;
        cols = m;
    } else {
        transCode = RINDOW_MATLIB_NO_TRANS;
        cols = n;
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,bufferX,offsetX)
            rindow_matlib_s_add(transCode, m, n, (float)alpha, pDataX, incX, pDataA, ldA);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,bufferX,offsetX)
            rindow_matlib_d_add(transCode, m, n, (double)alpha, pDataX, incX, pDataA, ldA);
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
    zend_long cols;
    zend_long transCode;

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
    if(trans) {
        transCode = RINDOW_MATLIB_TRANS;
        cols = m;
    } else {
        transCode = RINDOW_MATLIB_NO_TRANS;
        cols = n;
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,bufferX,offsetX)
            rindow_matlib_s_duplicate(transCode, m, n, pDataX, incX, pDataA, ldA);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,bufferX,offsetX)
            rindow_matlib_d_duplicate(transCode, m, n, pDataX, incX, pDataA, ldA);
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,buffer,offsetX)
            rindow_matlib_s_square(n, pDataX, incX);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,buffer,offsetX)
            rindow_matlib_d_square(n, pDataX, incX);
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,buffer,offsetX)
            rindow_matlib_s_sqrt(n, pDataX, incX);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,buffer,offsetX)
            rindow_matlib_d_sqrt(n, pDataX, incX);
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,buffer,offsetX)
            rindow_matlib_s_rsqrt(n, (float)alpha, pDataX, incX, (float)beta);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,buffer,offsetX)
            rindow_matlib_d_rsqrt(n, (double)alpha, pDataX, incX, (double)beta);
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
   A(m,n) := A(m,n) ** X(n)

   Method Rindow\OpenBLAS\Math::
    public function pow(
        int $m,
        int $n,
        Buffer $A, int $offsetA, int $ldA,
        Buffer $X, int $offsetX, int $incX,
        ) : void
 {{{ */
static PHP_METHOD(Math, pow)
{
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    php_interop_polite_math_matrix_linear_buffer_t* bufferX;
    zend_bool trans;
    zend_long m;
    zend_long n;
    zval* a=NULL;
    zend_long offsetA;
    zend_long ldA;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long i;
    zend_long cols;
    zend_long transCode;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 9, 9)
        Z_PARAM_BOOL(trans)
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
    if(trans) {
        transCode = RINDOW_MATLIB_TRANS;
        cols = m;
    } else {
        transCode = RINDOW_MATLIB_NO_TRANS;
        cols = n;
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
        "X", bufferX,cols,offsetX,incX)) {
        return;
    }
    // Check Buffer X and A
    if(bufferA->dtype!=bufferX->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for A and X", 0);
        return;
    }
    switch (bufferX->dtype) {
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,bufferX,offsetX)
            rindow_matlib_s_pow(transCode, m, n, pDataA, ldA, pDataX, incX);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,bufferX,offsetX)
            rindow_matlib_d_pow(transCode, m, n, pDataA, ldA, pDataX, incX);
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,buffer,offsetX)
            rindow_matlib_s_exp(n, pDataX, incX);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,buffer,offsetX)
            rindow_matlib_d_exp(n, pDataX, incX);
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,buffer,offsetX)
            rindow_matlib_s_log(n, pDataX, incX);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,buffer,offsetX)
            rindow_matlib_d_log(n, pDataX, incX);
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,buffer,offsetX)
            rindow_matlib_s_tanh(n, pDataX, incX);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,buffer,offsetX)
            rindow_matlib_d_tanh(n, pDataX, incX);
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,buffer,offsetX)
            rindow_matlib_s_sin(n, pDataX, incX);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,buffer,offsetX)
            rindow_matlib_d_sin(n, pDataX, incX);
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,buffer,offsetX)
            rindow_matlib_s_cos(n, pDataX, incX);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,buffer,offsetX)
            rindow_matlib_d_cos(n, pDataX, incX);
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,buffer,offsetX)
            rindow_matlib_s_tan(n, pDataX, incX);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,buffer,offsetX)
            rindow_matlib_d_tan(n, pDataX, incX);
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,buffer,offsetX)
            rindow_matlib_s_zeros(n, pDataX, incX);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,buffer,offsetX)
            rindow_matlib_d_zeros(n, pDataX, incX);
            break;
        }
        case php_interop_polite_math_matrix_dtype_int8:
        case php_interop_polite_math_matrix_dtype_uint8:
        case php_interop_polite_math_matrix_dtype_int16:
        case php_interop_polite_math_matrix_dtype_uint16:
        case php_interop_polite_math_matrix_dtype_int32:
        case php_interop_polite_math_matrix_dtype_uint32:
        case php_interop_polite_math_matrix_dtype_int64:
        case php_interop_polite_math_matrix_dtype_uint64:
        case php_interop_polite_math_matrix_dtype_bool: {
            int valueSize;
            void *pDataX;
            valueSize = php_rindow_openblas_common_dtype_to_valuesize(buffer->dtype);
            pDataX = php_rindow_openblas_get_address(buffer,offsetX,valueSize);
            rindow_matlib_i_zeros(buffer->dtype, n, pDataX, incX);
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

    // Check Buffer X
    if(bufferX->dtype==php_interop_polite_math_matrix_dtype_bool) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Data type of BufferX must not be bool", 0);
        return;
    }
    if(rindow_matlib_common_dtype_to_valuesize(bufferX->dtype)==0) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type of label number.", 0);
        return;
    }

    switch (bufferY->dtype) {
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataY,bufferY,offsetY)
            void* pDataX = rindow_matlib_common_get_address(bufferX->dtype, bufferX->data, offsetX);
            if(rindow_matlib_s_onehot(bufferX->dtype, m, n, pDataX, incX, (float)a, pDataY, ldY)) {
                zend_throw_exception(spl_ce_RuntimeException, "Label number is out of bounds.", 0);
                return;
            }
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataY,bufferY,offsetY)
            void* pDataX = rindow_matlib_common_get_address(bufferX->dtype, bufferX->data, offsetX);
            if(rindow_matlib_d_onehot(bufferX->dtype, m, n, pDataX, incX, (double)a, pDataY, ldY)) {
                zend_throw_exception(spl_ce_RuntimeException, "Label number is out of bounds.", 0);
                return;
            }
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataA,buffer,offsetA)
            rindow_matlib_s_softmax(m, n, pDataA, ldA);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataA,buffer,offsetA)
            rindow_matlib_d_softmax(m, n, pDataA, ldA);
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,bufferX,offsetX)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataY,bufferY,offsetY)
            rindow_matlib_s_equal(n, pDataX, incX, pDataY, incY);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,bufferX,offsetX)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataY,bufferY,offsetY)
            rindow_matlib_d_equal(n, pDataX, incX, pDataY, incY);
            break;
        }
        default: {
            if(!rindow_matlib_common_dtype_is_int(bufferX->dtype)&&
                !rindow_matlib_common_dtype_is_bool(bufferX->dtype)) {
                zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
                return;
            }
            void *pDataX = rindow_matlib_common_get_address(bufferX->dtype,bufferX->data,offsetX);
            void *pDataY = rindow_matlib_common_get_address(bufferY->dtype,bufferY->data,offsetY);
            rindow_matlib_i_equal(bufferX->dtype, n, pDataX, incX, pDataY, incY);
            break;
        }
    }
}
/* }}} */

/*
   Y(i) := 1  ( X(i) != Y(i) )
   Y(i) := 0  ( X(i) == Y(i) )

   Method Rindow\OpenBLAS\Math::
    public function notEqual(
        int $n,
        Buffer $X, int $offsetX, int $incX,
        Buffer $Y, int $offsetY, int $incY ) : void
 {{{ */
static PHP_METHOD(Math, notEqual)
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,bufferX,offsetX)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataY,bufferY,offsetY)
            rindow_matlib_s_notequal(n, pDataX, incX, pDataY, incY);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,bufferX,offsetX)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataY,bufferY,offsetY)
            rindow_matlib_d_notequal(n, pDataX, incX, pDataY, incY);
            break;
        }
        default: {
            if(!rindow_matlib_common_dtype_is_int(bufferX->dtype)&&
                !rindow_matlib_common_dtype_is_bool(bufferX->dtype)) {
                zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
                return;
            }
            void *pDataX = rindow_matlib_common_get_address(bufferX->dtype,bufferX->data,offsetX);
            void *pDataY = rindow_matlib_common_get_address(bufferY->dtype,bufferY->data,offsetY);
            rindow_matlib_i_notequal(bufferX->dtype, n, pDataX, incX, pDataY, incY);
            break;
        }
    }
}
/* }}} */

/*
   X(i) := 1  ( X(i) == 0 )
   X(i) := 0  ( X(i) != 0 )

   Method Rindow\OpenBLAS\Math::
    public function not(
        int $n,
        Buffer $X, int $offsetX, int $incX ) : void
 {{{ */
static PHP_METHOD(Math, not)
{
    php_interop_polite_math_matrix_linear_buffer_t* bufferX;
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
    // Check Buffer X
    bufferX = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(bufferX,"x")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", bufferX,n,offsetX,incX)) {
        return;
    }

    switch (bufferX->dtype) {
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,bufferX,offsetX)
            rindow_matlib_s_not(n, pDataX, incX);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,bufferX,offsetX)
            rindow_matlib_d_not(n, pDataX, incX);
            break;
        }
        default: {
            if(!rindow_matlib_common_dtype_is_int(bufferX->dtype)&&
                !rindow_matlib_common_dtype_is_bool(bufferX->dtype)) {
                zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
                return;
            }
            void *pDataX = rindow_matlib_common_get_address(bufferX->dtype,bufferX->data,offsetX);
            rindow_matlib_i_not(bufferX->dtype, n, pDataX, incX);
            break;
        }
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


    {
        void *pDataX = rindow_matlib_common_get_address(bufferX->dtype,bufferX->data,offsetX);
        void *pDataY = rindow_matlib_common_get_address(bufferY->dtype,bufferY->data,offsetY);

        if(rindow_matlib_astype(n, bufferX->dtype, pDataX, incX, bufferY->dtype, pDataY, incY)) {
            zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type of X or Y.", 0);
            return;
        }
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
    zend_long transCode;

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
    if(trans) {
        transCode = RINDOW_MATLIB_TRANS;
    } else {
        transCode = RINDOW_MATLIB_NO_TRANS;
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataB,bufferB,offsetB)
            rindow_matlib_s_matrixcopy(transCode, m, n, alpha, pDataA, ldA, pDataB, ldB);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataB,bufferB,offsetB)
            rindow_matlib_d_matrixcopy(transCode, m, n, alpha, pDataA, ldA, pDataB, ldB);
            break;
        }
        default: {
            zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type of A.", 0);
            return;
        }
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataB,bufferB,offsetB)
            rindow_matlib_s_imagecopy(height,width,channels,pDataA,pDataB,
                channelsFirst,heightShift,widthShift,verticalFlip,horizontalFlip,rgbFlip);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataB,bufferB,offsetB)
            rindow_matlib_d_imagecopy(height,width,channels,pDataA,pDataB,
                channelsFirst,heightShift,widthShift,verticalFlip,horizontalFlip,rgbFlip);
            break;
        }
        case php_interop_polite_math_matrix_dtype_uint8: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(uint8_t,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(uint8_t,pDataB,bufferB,offsetB)
            rindow_matlib_i8_imagecopy(height,width,channels,pDataA,pDataB,
                channelsFirst,heightShift,widthShift,verticalFlip,horizontalFlip,rgbFlip);
            break;
        }
        default: {
            zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type of A.", 0);
            return;
        }
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
        void *value = rindow_matlib_common_get_address(bufferV->dtype, bufferV->data, offsetV);
        void *x     = rindow_matlib_common_get_address(bufferX->dtype, bufferX->data, offsetX);
        rindow_matlib_fill(bufferX->dtype, n, value, x, incX);
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,buffer,offsetX)
            rindow_matlib_s_nan2num(n, pDataX, incX, (float)alpha);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,buffer,offsetX)
            rindow_matlib_d_nan2num(n, pDataX, incX, (double)alpha);
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,buffer,offsetX)
            rindow_matlib_s_isnan(n, pDataX, incX);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,buffer,offsetX)
            rindow_matlib_d_isnan(n, pDataX, incX);
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
   Y(n) := searchsorted( A(m,n), X(m) )

   Method Rindow\OpenBLAS\Math::
    public function searchsorted(
        int $m,
        int $n,
        Buffer $A, int $offsetA, int $ldA,
        Buffer $X, int $offsetX, int $incX,
        bool $right,
        Buffer $Y, int $offsetY int $incY, ) : void
 {{{ */

static PHP_METHOD(Math, searchsorted)
{
    zend_long m;
    zend_long n;
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    zend_long offsetA;
    zend_long ldA;
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
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT(a) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)

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
    if(bufferA->data==NULL) {
        zend_throw_exception_ex(spl_ce_InvalidArgumentException, 0, "uninitialized array: A");
        return;
    }
    if(offsetA<0) {
        zend_throw_exception_ex(spl_ce_InvalidArgumentException, 0, "Argument offsetA must be greater than or equals 0.");
        return;
    }
    if(ldA<0) {
        zend_throw_exception_ex(spl_ce_InvalidArgumentException, 0, "Argument ldA must be greater than or equals 0.");
        return;
    }
    if(offsetA+(m-1)*ldA+(n-1) >= bufferA->size) {
        zend_throw_exception_ex(spl_ce_InvalidArgumentException, 0, "Matrix specification too large for bufferA.");
        return;
    }

    // Check Buffer X
    bufferX = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_buffer_type(bufferX,"X")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "X", bufferX,m,offsetX,incX)) {
        return;
    }

    // Check Buffer Y
    bufferY = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(y);
    if(php_rindow_openblas_assert_buffer_type(bufferY,"Y")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "Y", bufferY,m,offsetY,incY)) {
        return;
    }

    // Check Buffer A and X
    if(bufferA->dtype!=bufferX->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for A and X", 0);
        return;
    }

    switch (bufferA->dtype) {
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,bufferX,offsetX)
            void *pDataY = rindow_matlib_common_get_address(bufferY->dtype, bufferY->data,offsetY);
            rindow_matlib_s_searchsorted(m,n,pDataA,ldA,pDataX,incX,right,bufferY->dtype,pDataY,incY);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,bufferX,offsetX)
            void *pDataY = rindow_matlib_common_get_address(bufferY->dtype, bufferY->data,offsetY);
            rindow_matlib_d_searchsorted(m,n,pDataA,ldA,pDataX,incX,right,bufferY->dtype,pDataY,incY);
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
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataX,bufferX,offsetX)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataY,bufferY,offsetY)
            rindow_matlib_s_cumsum(n,pDataX,incX,exclusive,reverse,pDataY,incY);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataX,bufferX,offsetX)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataY,bufferY,offsetY)
            rindow_matlib_d_cumsum(n,pDataX,incX,exclusive,reverse,pDataY,incY);
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
   Y(n) := X(n) + Y(n-1)

   Method Rindow\OpenBLAS\Math::
    public function transpose(
        Buffer $shape,
        Buffer $perm,
        Buffer $A, int $offsetA,
        Buffer $B, int $offsetB, 
        ) : void
 {{{ */

static PHP_METHOD(Math, transpose)
{
    zend_long ndim;
    php_interop_polite_math_matrix_linear_buffer_t* bufferShape;
    int32_t* shapevals;
    php_interop_polite_math_matrix_linear_buffer_t* bufferPerm;
    int32_t* permvals;
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    zend_long offsetA;
    php_interop_polite_math_matrix_linear_buffer_t* bufferB;
    zend_long offsetB;
    int32_t size;
    zval* shape=NULL;
    zval* perm=NULL;
    zval* a=NULL;
    zval* b=NULL;
    int32_t status;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 6, 6)
        Z_PARAM_OBJECT(shape) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_OBJECT(perm)  // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_OBJECT(a)     // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetA)
        Z_PARAM_OBJECT(b)     // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetB)
    ZEND_PARSE_PARAMETERS_END();

    // Check Buffer Shape
    bufferShape = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(shape);
    if(php_rindow_openblas_assert_buffer_type(bufferShape,"shape")) {
        return;
    }
    if(bufferShape->data==NULL) {
        zend_throw_exception_ex(spl_ce_DomainException, 0, "shapeBuffer is not initialized");
        return;
    }
    ndim = bufferShape->size;
    if(ndim<=0) {
        zend_throw_exception_ex(spl_ce_DomainException, 0, "ndim must be greater than 0.");
        return;
    }
    if(bufferShape->dtype!=php_interop_polite_math_matrix_dtype_int32) {
        zend_throw_exception_ex(spl_ce_InvalidArgumentException, 0, "data type of shape buffer must be int32.");
        return;
    }

    size = 1;
    shapevals=bufferShape->data; 
    for(int i=0;i<ndim;i++) {
        if(shapevals[i]<=0) {
            zend_throw_exception_ex(spl_ce_InvalidArgumentException, 0, "shape values must be greater than 0.");
            return;
        }
        size *= shapevals[i];
    }

    // Check Buffer perm
    bufferPerm = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(perm);
    if(php_rindow_openblas_assert_buffer_type(bufferPerm,"perm")) {
        return;
    }
    if(bufferPerm->data==NULL) {
        zend_throw_exception_ex(spl_ce_DomainException, 0, "bufferPerm is not initialized");
        return;
    }
    if(ndim != bufferPerm->size) {
        zend_throw_exception_ex(spl_ce_InvalidArgumentException, 0, "matrix shape and perm must be same size.");
        return;
    }
    if(bufferPerm->dtype!=php_interop_polite_math_matrix_dtype_int32) {
        zend_throw_exception_ex(spl_ce_InvalidArgumentException, 0, "data type of perm buffer must be int32.");
        return;
    }
    permvals=bufferPerm->data; 

    // Check Buffer A
    bufferA = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_buffer_type(bufferA,"A")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "A", bufferA,size,offsetA,1)) {
        return;
    }

    // Check Buffer B
    bufferB = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(b);
    if(php_rindow_openblas_assert_buffer_type(bufferB,"B")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "B", bufferB,size,offsetB,1)) {
        return;
    }

    // Check Buffer A and B
    if(bufferA->dtype!=bufferB->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for A and B.", 0);
        return;
    }

    switch (bufferA->dtype) {
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataB,bufferB,offsetB)
            status = rindow_matlib_s_transpose(ndim, shapevals, permvals, pDataA, pDataB);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataA,bufferA,offsetA)
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataB,bufferB,offsetB)
            status = rindow_matlib_d_transpose(ndim, shapevals, permvals, pDataA, pDataB);
            break;
        }
        case php_interop_polite_math_matrix_dtype_bool:
        case php_interop_polite_math_matrix_dtype_int8:
        case php_interop_polite_math_matrix_dtype_uint8:
        case php_interop_polite_math_matrix_dtype_int16:
        case php_interop_polite_math_matrix_dtype_uint16:
        case php_interop_polite_math_matrix_dtype_int32:
        case php_interop_polite_math_matrix_dtype_uint32:
        case php_interop_polite_math_matrix_dtype_int64:
        case php_interop_polite_math_matrix_dtype_uint64: {
            size_t value_bytes = rindow_matlib_common_dtype_to_valuesize(bufferA->dtype);
            void* pDataA = (int8_t*)(bufferA->data)+(value_bytes*offsetA);
            void* pDataB = (int8_t*)(bufferB->data)+(value_bytes*offsetB);
            status = rindow_matlib_i_transpose(bufferA->dtype, ndim, shapevals, permvals, pDataA, pDataB);
            break;
        }
        default: {
            zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
            return;
        }
    }
    if(status==RINDOW_MATLIB_SUCCESS) {
        return;
    }

    switch(status) {
        case RINDOW_MATLIB_E_MEM_ALLOC_FAILURE: {
            zend_throw_exception(spl_ce_RuntimeException, "memory allocation failure", 0);
            return;
        }
        case RINDOW_MATLIB_E_PERM_OUT_OF_RANGE: {
            zend_throw_exception(spl_ce_InvalidArgumentException, "perm contained an out-of-bounds axis", 0);
            return;
        }
        case RINDOW_MATLIB_E_DUP_AXIS: {
            zend_throw_exception(spl_ce_InvalidArgumentException, "Perm contained duplicate axis", 0);
            return;
        }
        default: {
            zend_throw_exception(spl_ce_RuntimeException, "Unknown error.", 0);
            return;
        }

    }
}

/* }}} */

/*
   A(m,n,k) :=  {in_band(n,k)==1: A(m,n,k) , in_band(n,k)==0: 0 }
    # in_band(n, k) = (lower < 0 || (n-k) <= lower)) && (upper < 0 || (k-n) <= upper)

   Method Rindow\OpenBLAS\Math::
    public function bandpart(
        int m,
        int n,
        int k,
        Buffer $A, int $offset,
        int $lower,
        int $upper,
        ) : void
 {{{ */

static PHP_METHOD(Math, bandpart)
{
    zend_long m,n,k;
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    zend_long offsetA;
    zend_long lower;
    zend_long upper;
    zval* a=NULL;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 7, 7)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_LONG(k)
        Z_PARAM_OBJECT(a) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetA)

        Z_PARAM_LONG(lower)
        Z_PARAM_LONG(upper)
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
    if(php_rindow_openblas_assert_buffer_type(bufferA,"A")) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        "A", bufferA,m*n*k,offsetA,1)) {
        return;
    }

    switch (bufferA->dtype) {
        case php_interop_polite_math_matrix_dtype_float32: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(float,pDataA,bufferA,offsetA)
            rindow_matlib_s_bandpart(m,n,k,pDataA,lower,upper);
            break;
        }
        case php_interop_polite_math_matrix_dtype_float64: {
            PHP_RINDOW_OPENBLAS_MATH_DEFDATA_TEMPLATE(double,pDataA,bufferA,offsetA)
            rindow_matlib_d_bandpart(m,n,k,pDataA,lower,upper);
            break;
        }
        default: {
            zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
            return;
        }
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

ZEND_BEGIN_ARG_INFO_EX(ai_Math_pow, 0, 0, 9)
    ZEND_ARG_INFO(0, trans)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, a, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_INFO(0, ldA)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
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

ZEND_BEGIN_ARG_INFO_EX(ai_Math_notEqual, 0, 0, 7)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_OBJ_INFO(0, y, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, incY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_not, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
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

ZEND_BEGIN_ARG_INFO_EX(ai_Math_transpose, 0, 0, 6)
    ZEND_ARG_OBJ_INFO(0, shape, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_OBJ_INFO(0, perm, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_OBJ_INFO(0, a, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_OBJ_INFO(0, b, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetB)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_bandpart, 0, 0, 7)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, k)
    ZEND_ARG_OBJ_INFO(0, a, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offset)
    ZEND_ARG_INFO(0, lower)
    ZEND_ARG_INFO(0, upper)
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
    PHP_ME(Math, updateAddOnehot,ai_Math_updateAddOnehot,ZEND_ACC_PUBLIC)
    PHP_ME(Math, softmax,        ai_Math_softmax,        ZEND_ACC_PUBLIC)
    PHP_ME(Math, equal,          ai_Math_equal,          ZEND_ACC_PUBLIC)
    PHP_ME(Math, notEqual,       ai_Math_notEqual,       ZEND_ACC_PUBLIC)
    PHP_ME(Math, not,            ai_Math_not,            ZEND_ACC_PUBLIC)
    PHP_ME(Math, astype,         ai_Math_astype,         ZEND_ACC_PUBLIC)
    PHP_ME(Math, matrixcopy,     ai_Math_matrixcopy,     ZEND_ACC_PUBLIC)
    PHP_ME(Math, imagecopy,      ai_Math_imagecopy,      ZEND_ACC_PUBLIC)
    PHP_ME(Math, fill,           ai_Math_fill,           ZEND_ACC_PUBLIC)
    PHP_ME(Math, nan2num,        ai_Math_nan2num,        ZEND_ACC_PUBLIC)
    PHP_ME(Math, isnan,          ai_Math_isnan,          ZEND_ACC_PUBLIC)
    PHP_ME(Math, searchsorted,   ai_Math_searchsorted,   ZEND_ACC_PUBLIC)
    PHP_ME(Math, cumsum,         ai_Math_cumsum,         ZEND_ACC_PUBLIC)
    PHP_ME(Math, transpose,      ai_Math_transpose,      ZEND_ACC_PUBLIC)
    PHP_ME(Math, bandpart,       ai_Math_bandpart,       ZEND_ACC_PUBLIC)
    PHP_ME(Math, gather,         ai_Math_gather,         ZEND_ACC_PUBLIC)
    PHP_ME(Math, reduceGather,   ai_Math_reduceGather,   ZEND_ACC_PUBLIC)
    PHP_ME(Math, repeat,         ai_Math_repeat,         ZEND_ACC_PUBLIC)
    PHP_ME(Math, slice,          ai_Math_slice,          ZEND_ACC_PUBLIC)
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
