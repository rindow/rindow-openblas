#include <php.h>
#include <Zend/zend_interfaces.h>
#include <Zend/zend_exceptions.h>
#include <ext/spl/spl_iterators.h>
#include <ext/spl/spl_exceptions.h>
#include <cblas.h>
#include <stdint.h>
#include <Rindow/OpenBLAS/Buffer.h>


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
        if(a<x[i*incX]) {
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

/* Method Rindow\OpenBLAS\Math::
    public function sum(
        int $n,
        Buffer $X, int $offsetX, int $incX ) : float
 {{{ */
static PHP_METHOD(Math, sum)
{
    php_rindow_openblas_buffer_t* buffer;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    double result;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "lOll",
    //        &n,&x,php_rindow_openblas_buffer_ce,&offsetX,&incX) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 4, 4)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    buffer = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        zend_long i;
        case php_rindow_openblas_dtype_float32:
            PHP_RINDOW_OPENBLAS_MATH_SUM_TEMPLATE(float)
            break;
        case php_rindow_openblas_dtype_float64:
            PHP_RINDOW_OPENBLAS_MATH_SUM_TEMPLATE(double)
            break;
        case php_rindow_openblas_dtype_int8:
            PHP_RINDOW_OPENBLAS_MATH_SUM_TEMPLATE(int8_t)
            break;
        case php_rindow_openblas_dtype_uint8:
            PHP_RINDOW_OPENBLAS_MATH_SUM_TEMPLATE(uint8_t)
            break;
        case php_rindow_openblas_dtype_int16:
            PHP_RINDOW_OPENBLAS_MATH_SUM_TEMPLATE(int16_t)
            break;
        case php_rindow_openblas_dtype_uint16:
            PHP_RINDOW_OPENBLAS_MATH_SUM_TEMPLATE(uint16_t)
            break;
        case php_rindow_openblas_dtype_int32:
            PHP_RINDOW_OPENBLAS_MATH_SUM_TEMPLATE(int32_t)
            break;
        case php_rindow_openblas_dtype_uint32:
            PHP_RINDOW_OPENBLAS_MATH_SUM_TEMPLATE(uint32_t)
            break;
        case php_rindow_openblas_dtype_int64:
            PHP_RINDOW_OPENBLAS_MATH_SUM_TEMPLATE(int64_t)
            break;
        case php_rindow_openblas_dtype_uint64:
            PHP_RINDOW_OPENBLAS_MATH_SUM_TEMPLATE(uint64_t)
            break;
        case php_rindow_openblas_dtype_bool:
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
            zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
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
    php_rindow_openblas_buffer_t* buffer;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    float  *pFloatX;
    double *pDoubleX;
    float  floatMax;
    double doubleMax;
    zend_long resultIdx;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "lOll",
    //        &n,&x,php_rindow_openblas_buffer_ce,&offsetX,&incX) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 4, 4)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    buffer = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        zend_long i;
        case php_rindow_openblas_dtype_float32:
            pFloatX = &(((float *)buffer->data)[offsetX]);
            floatMax = *pFloatX;
            pFloatX += incX;
            resultIdx = 0;
            for (i=1; i<n; i++,pFloatX+=incX) {
                if(floatMax < *pFloatX) {
                    floatMax = *pFloatX;
                    resultIdx = i;
                }
            }
            break;
        case php_rindow_openblas_dtype_float64:
            pDoubleX = &(((double *)buffer->data)[offsetX]);
            doubleMax = *pDoubleX;
            pDoubleX += incX;
            resultIdx = 0;
            for (i=1; i<n; i++,pDoubleX+=incX) {
                if(doubleMax < *pDoubleX) {
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
    php_rindow_openblas_buffer_t* buffer;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    float  *pFloatX;
    double *pDoubleX;
    float  floatMin;
    double doubleMin;
    zend_long resultIdx;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "lOll",
    //        &n,&x,php_rindow_openblas_buffer_ce,&offsetX,&incX) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 4, 4)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    buffer = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        zend_long i;
        case php_rindow_openblas_dtype_float32:
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
        case php_rindow_openblas_dtype_float64:
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
    php_rindow_openblas_buffer_t* buffer;
    zend_long n;
    double alpha;
    double beta;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long i;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "ldOlld",
    //        &n,&alpha,&x,php_rindow_openblas_buffer_ce,&offsetX,&incX,&beta) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 6, 6)
        Z_PARAM_LONG(n)
        Z_PARAM_DOUBLE(alpha)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_DOUBLE(beta)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    buffer = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_rindow_openblas_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    x[i*incX] = (float)alpha * x[i*incX] + (float)beta;
                }
            }
            break;
        case php_rindow_openblas_dtype_float64:
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
    php_rindow_openblas_buffer_t* buffer;
    zend_long n;
    double alpha;
    double beta;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long i;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "ldOlld",
    //        &n,&alpha,&x,php_rindow_openblas_buffer_ce,&offsetX,&incX,&beta) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 6, 6)
        Z_PARAM_LONG(n)
        Z_PARAM_DOUBLE(alpha)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_DOUBLE(beta)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    buffer = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_rindow_openblas_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    float t;
                    t = (float)alpha * x[i*incX] + (float)beta;
                    if(t==0.0) {
                        zend_throw_exception(spl_ce_RuntimeException, "Zero divide.", 0);
                        return;
                    }
                    x[i*incX] = 1 / t;
                }
            }
            break;
        case php_rindow_openblas_dtype_float64:
            {
                double *x = &(((double *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    double t;
                    t = (double)alpha * x[i*incX] + (double)beta;
                    if(t==0.0) {
                        zend_throw_exception(spl_ce_RuntimeException, "Zero divide.", 0);
                        return;
                    }
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
   X := X  (X > a)
   X := a  (X <= a)

   Method Rindow\OpenBLAS\Math::
    public function maximum(
        int $n,
        float $alpha,
        Buffer $X, int $offsetX, int $incX) : void
 {{{ */
static PHP_METHOD(Math, maximum)
{
    php_rindow_openblas_buffer_t* buffer;
    zend_long n;
    double alpha;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long i;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "ldOll",
    //        &n,&alpha,&x,php_rindow_openblas_buffer_ce,&offsetX,&incX) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 5, 5)
        Z_PARAM_LONG(n)
        Z_PARAM_DOUBLE(alpha)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    buffer = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_rindow_openblas_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    if((float)alpha > x[i*incX]) {
                        x[i*incX] = (float)alpha;
                    }
                }
            }
            break;
        case php_rindow_openblas_dtype_float64:
            {
                double *x = &(((double *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    if((double)alpha > x[i*incX]) {
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
    X := X  (X < a)
    X := a  (X >= a)

   Method Rindow\OpenBLAS\Math::
    public function minimum(
        int $n,
        float $alpha,
        Buffer $X, int $offsetX, int $incX) : void
 {{{ */
static PHP_METHOD(Math, minimum)
{
    php_rindow_openblas_buffer_t* buffer;
    zend_long n;
    double alpha;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long i;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "ldOll",
    //        &n,&alpha,&x,php_rindow_openblas_buffer_ce,&offsetX,&incX) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 5, 5)
        Z_PARAM_LONG(n)
        Z_PARAM_DOUBLE(alpha)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    buffer = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_rindow_openblas_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    if((float)alpha < x[i*incX]) {
                        x[i*incX] = (float)alpha;
                    }
                }
            }
            break;
        case php_rindow_openblas_dtype_float64:
            {
                double *x = &(((double *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    if((double)alpha < x[i*incX]) {
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
      X := 1  (X > a)
      X := 0  (X <= a)

   Method Rindow\OpenBLAS\Math::
    public function greater(
        int $n,
        float $alpha,
        Buffer $X, int $offsetX, int $incX) : void
 {{{ */
static PHP_METHOD(Math, greater)
{
    php_rindow_openblas_buffer_t* buffer;
    zend_long n;
    double alpha;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long i;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "ldOll",
    //        &n,&alpha,&x,php_rindow_openblas_buffer_ce,&offsetX,&incX) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 5, 5)
        Z_PARAM_LONG(n)
        Z_PARAM_DOUBLE(alpha)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    buffer = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_rindow_openblas_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    if(x[i*incX] > (float)alpha) {
                        x[i*incX] = 1.0;
                    } else {
                        x[i*incX] = 0.0;
                    }
                }
            }
            break;
        case php_rindow_openblas_dtype_float64:
            {
                double *x = &(((double *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    if(x[i*incX] > (double)alpha) {
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
        X := 1  (X < a)
        X := 0  (X >= a)

   Method Rindow\OpenBLAS\Math::
    public function less(
        int $n,
        float $alpha,
        Buffer $X, int $offsetX, int $incX) : void
 {{{ */
static PHP_METHOD(Math, less)
{
    php_rindow_openblas_buffer_t* buffer;
    zend_long n;
    double alpha;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long i;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "ldOll",
    //        &n,&alpha,&x,php_rindow_openblas_buffer_ce,&offsetX,&incX) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 5, 5)
        Z_PARAM_LONG(n)
        Z_PARAM_DOUBLE(alpha)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    buffer = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_rindow_openblas_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    if(x[i*incX] < (float)alpha) {
                        x[i*incX] = 1.0;
                    } else {
                        x[i*incX] = 0.0;
                    }
                }
            }
            break;
        case php_rindow_openblas_dtype_float64:
            {
                double *x = &(((double *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    if(x[i*incX] < (double)alpha) {
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
    php_rindow_openblas_buffer_t* bufferX;
    php_rindow_openblas_buffer_t* bufferA;
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

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "bllOllOll",
    //        &broadcast,
    //        &m,
    //        &n,
    //        &x,php_rindow_openblas_buffer_ce,&offsetX,&incX,
    //        &a,php_rindow_openblas_buffer_ce,&offsetA,&ldA) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 9, 9)
        Z_PARAM_BOOL(trans)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_OBJECT_OF_CLASS(a,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_M, m)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    // Check Buffer X
    bufferX = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(!trans) {
        rows = m; cols = n;
    } else {
        rows = n; cols = m;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, bufferX,cols,offsetX,incX)) {
        return;
    }

    // Check Buffer A
    bufferA = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_A, bufferA,m,n,offsetA,ldA)) {
        return;
    }

    // Check Buffer X and A
    if(bufferX->dtype!=bufferA->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for X and A", 0);
        return;
    }

    switch (bufferX->dtype) {
        case php_rindow_openblas_dtype_float32:
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
        case php_rindow_openblas_dtype_float64:
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
    php_rindow_openblas_buffer_t* bufferX;
    php_rindow_openblas_buffer_t* bufferA;
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

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "blldOllOll",
    //        &broadcast,
    //        &m,
    //        &n,
    //        &alpha,
    //        &x,php_rindow_openblas_buffer_ce,&offsetX,&incX,
    //        &a,php_rindow_openblas_buffer_ce,&offsetA,&ldA) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 10, 10)
        Z_PARAM_BOOL(trans)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_DOUBLE(alpha)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_OBJECT_OF_CLASS(a,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_M, m)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    // Check Buffer X
    bufferX = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(!trans) {
        rows = m; cols = n;
    } else {
        rows = n; cols = m;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, bufferX,cols,offsetX,incX)) {
        return;
    }

    // Check Buffer A
    bufferA = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_A, bufferA,m,n,offsetA,ldA)) {
        return;
    }

    // Check Buffer X and A
    if(bufferX->dtype!=bufferA->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for X and A", 0);
        return;
    }

    switch (bufferX->dtype) {
        case php_rindow_openblas_dtype_float32:
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
        case php_rindow_openblas_dtype_float64:
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
    php_rindow_openblas_buffer_t* bufferX;
    php_rindow_openblas_buffer_t* bufferA;
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

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "llOllOll",
    //        &m,
    //        &n,
    //        &x,php_rindow_openblas_buffer_ce,&offsetX,&incX,
    //        &y,php_rindow_openblas_buffer_ce,&offsetY,&incY) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 9, 9)
        Z_PARAM_BOOL(trans)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_OBJECT_OF_CLASS(a,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_M, m)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    if(!trans) {
        rows = m; cols = n;
    } else {
        rows = n; cols = m;
    }
    // Check Buffer X
    bufferX = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, bufferX,cols,offsetX,incX)) {
        return;
    }

    // Check Buffer A
    bufferA = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_A, bufferA,m,n,offsetA,ldA)) {
        return;
    }

    // Check Buffer X and Y
    if(bufferX->dtype!=bufferA->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for X and Y", 0);
        return;
    }

    switch (bufferX->dtype) {
        case php_rindow_openblas_dtype_float32:
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
        case php_rindow_openblas_dtype_float64:
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
    php_rindow_openblas_buffer_t* buffer;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long i;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "lOll",
    //        &n,&x,php_rindow_openblas_buffer_ce,&offsetX,&incX) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 4, 4)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    buffer = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_rindow_openblas_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    x[i*incX] = x[i*incX] * x[i*incX];
                }
            }
            break;
        case php_rindow_openblas_dtype_float64:
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
    php_rindow_openblas_buffer_t* buffer;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long i;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "lOll",
    //        &n,&x,php_rindow_openblas_buffer_ce,&offsetX,&incX) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 4, 4)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    buffer = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_rindow_openblas_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    float t;
                    t = x[i*incX];
                    if(t<0.0) {
                        zend_throw_exception(spl_ce_RuntimeException, "Invalid value in sqrt.", 0);
                        return;
                    }
                    x[i*incX] = sqrtf(t);
                }
            }
            break;
        case php_rindow_openblas_dtype_float64:
            {
                double *x = &(((double *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    double t;
                    t = x[i*incX];
                    if(t<0.0) {
                        zend_throw_exception(spl_ce_RuntimeException, "Invalid value in sqrt.", 0);
                        return;
                    }
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
    php_rindow_openblas_buffer_t* buffer;
    zend_long n;
    double alpha;
    double beta;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long i;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "ldOlld",
    //        &n,&alpha,&x,php_rindow_openblas_buffer_ce,&offsetX,&incX,&beta) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 6, 6)
        Z_PARAM_LONG(n)
        Z_PARAM_DOUBLE(alpha)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_DOUBLE(beta)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    buffer = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_rindow_openblas_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    float t;
                    if(x[i*incX]<0.0) {
                        zend_throw_exception(spl_ce_RuntimeException, "Invalid value in sqrt.", 0);
                        return;
                    }
                    t = (float)alpha * sqrtf(x[i*incX]) + (float)beta;
                    if(t==0.0) {
                        zend_throw_exception(spl_ce_RuntimeException, "Zero divide.", 0);
                        return;
                    }
                    x[i*incX] = 1 / t;
                }
            }
            break;
        case php_rindow_openblas_dtype_float64:
            {
                double *x = &(((double *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    double t;
                    if(x[i*incX]<0.0) {
                        zend_throw_exception(spl_ce_RuntimeException, "Invalid value in sqrt.", 0);
                        return;
                    }
                    t = (double)alpha * sqrt(x[i*incX]) + (double)beta;
                    if(t==0.0) {
                        zend_throw_exception(spl_ce_RuntimeException, "Zero divide.", 0);
                        return;
                    }
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
        float $alpha,
        Buffer $X, int $offsetX, int $incX) : void
 {{{ */
static PHP_METHOD(Math, pow)
{
    php_rindow_openblas_buffer_t* buffer;
    zend_long n;
    double alpha;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long i;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "ldOll",
    //        &n,&alpha,&x,php_rindow_openblas_buffer_ce,&offsetX,&incX) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 5, 5)
        Z_PARAM_LONG(n)
        Z_PARAM_DOUBLE(alpha)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    buffer = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_rindow_openblas_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    x[i*incX] = powf(x[i*incX], (float)alpha);
                }
            }
            break;
        case php_rindow_openblas_dtype_float64:
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
    php_rindow_openblas_buffer_t* buffer;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long i;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "lOll",
    //        &n,&x,php_rindow_openblas_buffer_ce,&offsetX,&incX) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 4, 4)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    buffer = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_rindow_openblas_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    x[i*incX] = expf(x[i*incX]);
                }
            }
            break;
        case php_rindow_openblas_dtype_float64:
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
    php_rindow_openblas_buffer_t* buffer;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long i;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "lOll",
    //        &n,&x,php_rindow_openblas_buffer_ce,&offsetX,&incX) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 4, 4)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    buffer = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_rindow_openblas_dtype_float32:
            {
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    float t;
                    t = x[i*incX];
                    if(t<=0.0) {
                        zend_throw_exception(spl_ce_RuntimeException, "Invalid value in log.", 0);
                        return;
                    }
                    x[i*incX] = logf(t);
                }
            }
            break;
        case php_rindow_openblas_dtype_float64:
            {
                double *x = &(((double *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    double t;
                    t = x[i*incX];
                    if(t<=0.0) {
                        zend_throw_exception(spl_ce_RuntimeException, "Invalid value in log.", 0);
                        return;
                    }
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
   X := 0

   Method Rindow\OpenBLAS\Math::
    public function zeros(
        int $n,
        Buffer $X, int $offsetX, int $incX) : void
 {{{ */
static PHP_METHOD(Math, zeros)
{
    php_rindow_openblas_buffer_t* buffer;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "lOll",
    //        &n,&x,php_rindow_openblas_buffer_ce,&offsetX,&incX) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 4, 4)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    buffer = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_rindow_openblas_dtype_float32:
            {
                zend_long i;
                float *x = &(((float *)buffer->data)[offsetX]);
                for(i=0;i<n;i++) {
                    x[i*incX] = 0;
                }
            }
            break;
        case php_rindow_openblas_dtype_float64:
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
                valueSize = php_rindow_openblas_dtype_to_valuesize(buffer->dtype);
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
        case php_rindow_openblas_dtype_bool:
            {
                uint8_t *x = &(((uint8_t *)buffer)[offset]);
                if(x[index*incWidth]==0) { *value=0; return 0; }
                else                     { *value=1; return 0; }
            }
        case php_rindow_openblas_dtype_int8:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(int8_t,zend_long);
        case php_rindow_openblas_dtype_uint8:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(uint8_t,zend_long);
        case php_rindow_openblas_dtype_int16:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(int16_t,zend_long);
        case php_rindow_openblas_dtype_uint16:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(uint16_t,zend_long);
        case php_rindow_openblas_dtype_int32:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(int32_t,zend_long);
        case php_rindow_openblas_dtype_uint32:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(uint32_t,zend_long);
        case php_rindow_openblas_dtype_int64:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(int64_t,zend_long);
        case php_rindow_openblas_dtype_uint64:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(uint64_t,zend_long);
        case php_rindow_openblas_dtype_float32:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(float,zend_long);
        case php_rindow_openblas_dtype_float64:
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
        case php_rindow_openblas_dtype_bool:
        {
            uint8_t *x = &(((uint8_t *)buffer)[offset]);
            if(value==0) { x[index*incWidth]=0; return 0; }
            else         { x[index*incWidth]=1; return 0; }
        }
        case php_rindow_openblas_dtype_int8:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(int8_t);
        case php_rindow_openblas_dtype_uint8:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(uint8_t);
        case php_rindow_openblas_dtype_int16:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(int16_t);
        case php_rindow_openblas_dtype_uint16:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(uint16_t);
        case php_rindow_openblas_dtype_int32:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(int32_t);
        case php_rindow_openblas_dtype_uint32:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(uint32_t);
        case php_rindow_openblas_dtype_int64:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(int64_t);
        case php_rindow_openblas_dtype_uint64:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(uint64_t);
        case php_rindow_openblas_dtype_float32:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(float);
        case php_rindow_openblas_dtype_float64:
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
        case php_rindow_openblas_dtype_bool:
            {
                uint8_t *x = &(((uint8_t *)buffer)[offset]);
                if(x[index*incWidth]==0) { *value=0; return 0; }
                else                     { *value=1; return 0; }
            }
        case php_rindow_openblas_dtype_int8:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(int8_t,double);
        case php_rindow_openblas_dtype_uint8:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(uint8_t,double);
        case php_rindow_openblas_dtype_int16:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(int16_t,double);
        case php_rindow_openblas_dtype_uint16:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(uint16_t,double);
        case php_rindow_openblas_dtype_int32:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(int32_t,double);
        case php_rindow_openblas_dtype_uint32:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(uint32_t,double);
        case php_rindow_openblas_dtype_int64:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(int64_t,double);
        case php_rindow_openblas_dtype_uint64:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(uint64_t,double);
        case php_rindow_openblas_dtype_float32:
            RINDOW_OPENBLAS_MATH_GET_CAST_TEMPLATE(float,double);
        case php_rindow_openblas_dtype_float64:
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
        case php_rindow_openblas_dtype_bool:
        {
            uint8_t *x = &(((uint8_t *)buffer)[offset]);
            if(value==0) { x[index*incWidth]=0; return 0; }
            else         { x[index*incWidth]=1; return 0; }
        }
        case php_rindow_openblas_dtype_int8:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(int8_t);
        case php_rindow_openblas_dtype_uint8:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(uint8_t);
        case php_rindow_openblas_dtype_int16:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(int16_t);
        case php_rindow_openblas_dtype_uint16:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(uint16_t);
        case php_rindow_openblas_dtype_int32:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(int32_t);
        case php_rindow_openblas_dtype_uint32:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(uint32_t);
        case php_rindow_openblas_dtype_int64:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(int64_t);
        case php_rindow_openblas_dtype_uint64:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(uint64_t);
        case php_rindow_openblas_dtype_float32:
            RINDOW_OPENBLAS_MATH_SET_CAST_TEMPLATE(float);
        case php_rindow_openblas_dtype_float64:
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
    php_rindow_openblas_buffer_t* bufferX;
    php_rindow_openblas_buffer_t* bufferY;
    zend_long m;
    zend_long n;
    double a;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zval* y=NULL;
    zend_long offsetY;
    zend_long ldY;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "lldOllOll",
    //        &m,
    //        &n,
    //        &a,
    //        &x,php_rindow_openblas_buffer_ce,&offsetX,&incX,
    //        &y,php_rindow_openblas_buffer_ce,&offsetY,&ldY) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 9, 9)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_DOUBLE(a)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_OBJECT_OF_CLASS(y,php_rindow_openblas_buffer_ce)
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
    // Check Buffer X
    bufferX = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, bufferX,m,offsetX,incX)) {
        return;
    }

    // Check Buffer Y
    bufferY = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(y);
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_Y, bufferY,m,n,offsetY,ldY)) {
        return;
    }

    // Check Buffer A and Y
    if(bufferX->dtype==php_rindow_openblas_dtype_bool) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Data type of BufferX must not be bool", 0);
        return;
    }

    switch (bufferY->dtype) {
        case php_rindow_openblas_dtype_float32:
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
        case php_rindow_openblas_dtype_float64:
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
    php_rindow_openblas_buffer_t* buffer;
    zend_long m;
    zend_long n;
    zval* a=NULL;
    zend_long offsetA;
    zend_long ldA;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 5, 5)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT_OF_CLASS(a,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_M, m)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    buffer = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_A, buffer,m,n,offsetA,ldA)) {
        return;
    }
    switch (buffer->dtype) {
        case php_rindow_openblas_dtype_float32:
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
                    if(sum_exp==0.0) {
                        zend_throw_exception(spl_ce_RuntimeException, "Zero divide in softmax.", 0);
                        return;
                    }
                    for(j=0;j<n;j++) {
                        a[j] = a[j] / sum_exp;
                    }
                }
            }
            break;
        case php_rindow_openblas_dtype_float64:
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
                    if(sum_exp==0.0) {
                        zend_throw_exception(spl_ce_RuntimeException, "Zero divide in softmax.", 0);
                        return;
                    }
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
    php_rindow_openblas_buffer_t* bufferX;
    php_rindow_openblas_buffer_t* bufferY;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zval* y=NULL;
    zend_long offsetY;
    zend_long incY;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "bllOllOll",
    //        &broadcast,
    //        &m,
    //        &n,
    //        &x,php_rindow_openblas_buffer_ce,&offsetX,&incX,
    //        &y,php_rindow_openblas_buffer_ce,&offsetY,&incY) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 7, 7)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_OBJECT_OF_CLASS(y,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetY)
        Z_PARAM_LONG(incY)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    // Check Buffer X
    bufferX = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, bufferX,n,offsetX,incX)) {
        return;
    }

    // Check Buffer Y
    bufferY = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(y);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_Y, bufferY,n,offsetY,incY)) {
        return;
    }

    // Check Buffer X and Y
    if(bufferX->dtype!=bufferY->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for X and Y", 0);
        return;
    }

    switch (bufferX->dtype) {
        case php_rindow_openblas_dtype_float32:
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
        case php_rindow_openblas_dtype_float64:
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
            if(!php_rindow_openblas_dtype_is_int(bufferX->dtype)&&
                !php_rindow_openblas_dtype_is_bool(bufferX->dtype)) {
                zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
                return;
            }
            {
                int valueSize;
                uint8_t *x, *y;
                zend_long i,zero,one;
                zero = 0;
                one = 1;
                valueSize = php_rindow_openblas_dtype_to_valuesize(bufferX->dtype);
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
   X(m) := sum( A(m,n) )

   Method Rindow\OpenBLAS\Math::
    public function reduceSum(
        bool $trans,
        int $m,
        int $n,
        Buffer $A, int $offsetA, int $ldA,
        Buffer $X, int $offsetX, int $incX ) : void
 {{{ */
static PHP_METHOD(Math, reduceSum)
{
    php_rindow_openblas_buffer_t* bufferA;
    php_rindow_openblas_buffer_t* bufferX;
    zend_bool trans;
    zend_long m;
    zend_long n;
    zval* a=NULL;
    zend_long offsetA;
    zend_long ldA;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long rows,cols;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "bllOllOll",
    //        &broadcast,
    //        &m,
    //        &n,
    //        &x,php_rindow_openblas_buffer_ce,&offsetX,&incX,
    //        &a,php_rindow_openblas_buffer_ce,&offsetA,&ldA) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 9, 9)
        Z_PARAM_BOOL(trans)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT_OF_CLASS(a,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_M, m)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    // Check Buffer X
    bufferX = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(!trans) {
        rows = m; cols = n;
    } else {
        rows = n; cols = m;
    }

    // Check Buffer A
    bufferA = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_A, bufferA,m,n,offsetA,ldA)) {
        return;
    }

    // Check Buffer X
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, bufferX,rows,offsetX,incX)) {
        return;
    }

    // Check Buffer X and A
    if(bufferX->dtype!=bufferA->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for X and A", 0);
        return;
    }

    switch (bufferX->dtype) {
        case php_rindow_openblas_dtype_float32:
            {
                float *x = &(((float *)bufferX->data)[offsetX]);
                float *a = &(((float *)bufferA->data)[offsetA]);
                zend_long j,incAj,incAi;
                if(!trans) { incAj = ldA; incAi = 1;}
                else       { incAj = 1;   incAi = ldA;}
                for(j=0; j<rows; j++) {
                    x[j*incX] = s_sum(cols,&a[j*incAj],incAi);
                }
            }
            break;
        case php_rindow_openblas_dtype_float64:
            {
                double *x = &(((double *)bufferX->data)[offsetX]);
                double *a = &(((double *)bufferA->data)[offsetA]);
                zend_long j,incAj,incAi;
                if(!trans) { incAj = ldA; incAi = 1;}
                else       { incAj = 1;   incAi = ldA;}
                for(j=0; j<rows; j++) {
                    x[j*incX] = d_sum(cols,&a[j*incAj],incAi);
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
        bool $trans,
        int $m,
        int $n,
        Buffer $A, int $offsetA, int $ldA,
        Buffer $X, int $offsetX, int $incX ) : void
 {{{ */
static PHP_METHOD(Math, reduceMax)
{
    php_rindow_openblas_buffer_t* bufferA;
    php_rindow_openblas_buffer_t* bufferX;
    zend_bool trans;
    zend_long m;
    zend_long n;
    zval* a=NULL;
    zend_long offsetA;
    zend_long ldA;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long rows,cols;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "bllOllOll",
    //        &broadcast,
    //        &m,
    //        &n,
    //        &x,php_rindow_openblas_buffer_ce,&offsetX,&incX,
    //        &a,php_rindow_openblas_buffer_ce,&offsetA,&ldA) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 9, 9)
        Z_PARAM_BOOL(trans)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT_OF_CLASS(a,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_M, m)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    // Check Buffer X
    bufferX = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(!trans) {
        rows = m; cols = n;
    } else {
        rows = n; cols = m;
    }

    // Check Buffer A
    bufferA = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_A, bufferA,m,n,offsetA,ldA)) {
        return;
    }

    // Check Buffer X
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, bufferX,rows,offsetX,incX)) {
        return;
    }

    // Check Buffer X and A
    if(bufferX->dtype!=bufferA->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for X and A", 0);
        return;
    }

    switch (bufferX->dtype) {
        case php_rindow_openblas_dtype_float32:
            {
                float *x = &(((float *)bufferX->data)[offsetX]);
                float *a = &(((float *)bufferA->data)[offsetA]);
                zend_long j,incAj,incAi;
                if(!trans) { incAj = ldA; incAi = 1;}
                else       { incAj = 1;   incAi = ldA;}
                for(j=0; j<rows; j++) {
                    x[j*incX] = s_max(cols,&a[j*incAj],incAi);
                }
            }
            break;
        case php_rindow_openblas_dtype_float64:
            {
                double *x = &(((double *)bufferX->data)[offsetX]);
                double *a = &(((double *)bufferA->data)[offsetA]);
                zend_long j,incAj,incAi;
                if(!trans) { incAj = ldA; incAi = 1;}
                else       { incAj = 1;   incAi = ldA;}
                for(j=0; j<rows; j++) {
                    x[j*incX] = d_max(cols,&a[j*incAj],incAi);
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
    public function reduceMax(
        bool $trans,
        int $m,
        int $n,
        Buffer $A, int $offsetA, int $ldA,
        Buffer $X, int $offsetX, int $incX ) : void
 {{{ */
static PHP_METHOD(Math, reduceArgMax)
{
    php_rindow_openblas_buffer_t* bufferA;
    php_rindow_openblas_buffer_t* bufferX;
    zend_bool trans;
    zend_long m;
    zend_long n;
    zval* a=NULL;
    zend_long offsetA;
    zend_long ldA;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zend_long rows,cols;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "bllOllOll",
    //        &broadcast,
    //        &m,
    //        &n,
    //        &x,php_rindow_openblas_buffer_ce,&offsetX,&incX,
    //        &a,php_rindow_openblas_buffer_ce,&offsetA,&ldA) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 9, 9)
        Z_PARAM_BOOL(trans)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT_OF_CLASS(a,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_M, m)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    // Check Buffer X
    bufferX = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(!trans) {
        rows = m; cols = n;
    } else {
        rows = n; cols = m;
    }

    // Check Buffer A
    bufferA = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_A, bufferA,m,n,offsetA,ldA)) {
        return;
    }

    // Check Buffer X
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, bufferX,rows,offsetX,incX)) {
        return;
    }

    switch (bufferA->dtype) {
        case php_rindow_openblas_dtype_float32:
            {
                zend_long value;
                float *a = &(((float *)bufferA->data)[offsetA]);
                zend_long j,incAj,incAi;
                if(!trans) { incAj = ldA; incAi = 1;}
                else       { incAj = 1;   incAi = ldA;}
                for(j=0; j<rows; j++) {
                    value = s_argmax(cols,&a[j*incAj],incAi);
                    rindow_openblas_math_set_integer(bufferX->dtype, bufferX->data, offsetX, incX, j, value);
                }
            }
            break;
        case php_rindow_openblas_dtype_float64:
            {
                zend_long value;
                double *a = &(((double *)bufferA->data)[offsetA]);
                zend_long j,incAj,incAi;
                if(!trans) { incAj = ldA; incAi = 1;}
                else       { incAj = 1;   incAi = ldA;}
                for(j=0; j<rows; j++) {
                    value = d_argmax(cols,&a[j*incAj],incAi);
                    rindow_openblas_math_set_integer(bufferX->dtype, bufferX->data, offsetX, incX, j, value);
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
    php_rindow_openblas_buffer_t* bufferX;
    php_rindow_openblas_buffer_t* bufferY;
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
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_OBJECT_OF_CLASS(y,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetY)
        Z_PARAM_LONG(incY)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    bufferX = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, bufferX,n,offsetX,incX)) {
        return;
    }
    bufferY = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(y);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_Y, bufferY,n,offsetY,incY)) {
        return;
    }
    // Check dtype and Buffer Y
    if(dtype!=bufferY->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for Y", 0);
        return;
    }

    if(php_rindow_openblas_dtype_is_int(dtype) || php_rindow_openblas_dtype_is_bool(dtype)) {
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
    } else if(php_rindow_openblas_dtype_is_float(dtype)) {
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


#include "Math_select.c"
#include "Math_scatter.c"
#include "Math_im2col1d.c"
#include "Math_im2col2d.c"
#include "Math_im2col3d.c"

ZEND_BEGIN_ARG_INFO_EX(ai_Math_sum, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_imax, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_imin, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_increment, 0, 0, 6)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, alpha)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_INFO(0, beta)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_reciprocal, 0, 0, 6)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, alpha)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_INFO(0, beta)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_maximum, 0, 0, 5)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, alpha)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_minimum, 0, 0, 5)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, alpha)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_greater, 0, 0, 5)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, alpha)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_less, 0, 0, 5)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, alpha)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_multiply, 0, 0, 9)
    ZEND_ARG_INFO(0, trans)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_OBJ_INFO(0, y, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, incY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_add, 0, 0, 10)
    ZEND_ARG_INFO(0, trans)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, alpha)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_OBJ_INFO(0, y, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, incY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_duplicate, 0, 0, 9)
    ZEND_ARG_INFO(0, trans)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_OBJ_INFO(0, y, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, incY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_square, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_sqrt, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_rsqrt, 0, 0, 6)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, alpha)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_INFO(0, beta)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_pow, 0, 0, 5)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, alpha)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_exp, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_log, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_zeros, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_selectAxis0, 0, 0, 12)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, k)
    ZEND_ARG_OBJ_INFO(0, a, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_INFO(0, ldA)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_OBJ_INFO(0, y, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, ldY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_selectAxis1, 0, 0, 11)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, a, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_INFO(0, ldA)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_OBJ_INFO(0, y, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, incY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_scatterAxis0, 0, 0, 12)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, k)
    ZEND_ARG_OBJ_INFO(0, a, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_INFO(0, ldA)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_OBJ_INFO(0, y, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, ldY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_scatterAxis1, 0, 0, 11)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, a, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_INFO(0, ldA)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_OBJ_INFO(0, y, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, incY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_updateAddOnehot, 0, 0, 9)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, a)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_OBJ_INFO(0, y, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, ldY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_softmax, 0, 0, 5)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_equal, 0, 0, 7)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_OBJ_INFO(0, y, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, incY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_reduceSum, 0, 0, 9)
    ZEND_ARG_INFO(0, trans)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_OBJ_INFO(0, y, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, incY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_reduceMax, 0, 0, 9)
    ZEND_ARG_INFO(0, trans)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_OBJ_INFO(0, y, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, incY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_reduceArgMax, 0, 0, 9)
    ZEND_ARG_INFO(0, trans)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_OBJ_INFO(0, y, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, incY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_astype, 0, 0, 8)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, dtype)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_OBJ_INFO(0, y, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, incY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_im2col1d, 0, 0, 15)
    ZEND_ARG_INFO(0, reverse)
    ZEND_ARG_OBJ_INFO(0, images_obj, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, images_offset)
    ZEND_ARG_INFO(0, images_size)
    ZEND_ARG_INFO(0, batches)

    ZEND_ARG_INFO(0, im_w)
    ZEND_ARG_INFO(0, channels)
    ZEND_ARG_INFO(0, filter_w)
    ZEND_ARG_INFO(0, stride_w)
    ZEND_ARG_INFO(0, padding)

    ZEND_ARG_INFO(0, channels_first)
    ZEND_ARG_INFO(0, cols_channels_first)
    ZEND_ARG_OBJ_INFO(0, cols_obj,Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, cols_offset)
    ZEND_ARG_INFO(0, cols_size)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_im2col2d, 0, 0, 18)
    ZEND_ARG_INFO(0, reverse)
    ZEND_ARG_OBJ_INFO(0, images_obj, Rindow\\OpenBLAS\\Buffer, 0)
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
    ZEND_ARG_INFO(0, cols_channels_first)

    ZEND_ARG_OBJ_INFO(0, cols_obj,Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, cols_offset)
    ZEND_ARG_INFO(0, cols_size)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Math_im2col3d, 0, 0, 21)
    ZEND_ARG_INFO(0, reverse)
    ZEND_ARG_OBJ_INFO(0, images_obj, Rindow\\OpenBLAS\\Buffer, 0)
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
    ZEND_ARG_INFO(0, cols_channels_first)
    ZEND_ARG_OBJ_INFO(0, cols_obj,Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, cols_offset)

    ZEND_ARG_INFO(0, cols_size)
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
    PHP_ME(Math, less,           ai_Math_less,           ZEND_ACC_PUBLIC)
    PHP_ME(Math, multiply,       ai_Math_multiply,       ZEND_ACC_PUBLIC)
    PHP_ME(Math, add,            ai_Math_add,            ZEND_ACC_PUBLIC)
    PHP_ME(Math, duplicate,      ai_Math_duplicate,      ZEND_ACC_PUBLIC)
    PHP_ME(Math, square,         ai_Math_square,         ZEND_ACC_PUBLIC)
    PHP_ME(Math, sqrt,           ai_Math_sqrt,           ZEND_ACC_PUBLIC)
    PHP_ME(Math, rsqrt,          ai_Math_rsqrt,          ZEND_ACC_PUBLIC)
    PHP_ME(Math, pow,            ai_Math_pow,            ZEND_ACC_PUBLIC)
    PHP_ME(Math, exp,            ai_Math_exp,            ZEND_ACC_PUBLIC)
    PHP_ME(Math, log,            ai_Math_log,            ZEND_ACC_PUBLIC)
    PHP_ME(Math, zeros,          ai_Math_zeros,          ZEND_ACC_PUBLIC)
    PHP_ME(Math, selectAxis0,    ai_Math_selectAxis0,    ZEND_ACC_PUBLIC)
    PHP_ME(Math, selectAxis1,    ai_Math_selectAxis1,    ZEND_ACC_PUBLIC)
    PHP_ME(Math, scatterAxis0,    ai_Math_scatterAxis0,    ZEND_ACC_PUBLIC)
    PHP_ME(Math, scatterAxis1,    ai_Math_scatterAxis1,    ZEND_ACC_PUBLIC)
    PHP_ME(Math, updateAddOnehot,ai_Math_updateAddOnehot,ZEND_ACC_PUBLIC)
    PHP_ME(Math, softmax,        ai_Math_softmax,        ZEND_ACC_PUBLIC)
    PHP_ME(Math, equal,          ai_Math_equal,          ZEND_ACC_PUBLIC)
    PHP_ME(Math, reduceSum,      ai_Math_reduceSum,      ZEND_ACC_PUBLIC)
    PHP_ME(Math, reduceMax,      ai_Math_reduceMax,      ZEND_ACC_PUBLIC)
    PHP_ME(Math, reduceArgMax,      ai_Math_reduceArgMax,      ZEND_ACC_PUBLIC)
    PHP_ME(Math, astype,         ai_Math_astype,         ZEND_ACC_PUBLIC)
    PHP_ME(Math, im2col1d,         ai_Math_im2col1d,         ZEND_ACC_PUBLIC)
    PHP_ME(Math, im2col2d,         ai_Math_im2col2d,         ZEND_ACC_PUBLIC)
    PHP_ME(Math, im2col3d,         ai_Math_im2col3d,         ZEND_ACC_PUBLIC)
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
