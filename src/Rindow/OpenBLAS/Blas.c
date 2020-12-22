#include <php.h>
#include <Zend/zend_interfaces.h>
#include <Zend/zend_exceptions.h>
#include <ext/spl/spl_iterators.h>
#include <ext/spl/spl_exceptions.h>
#include <cblas.h>
#include <Rindow/OpenBLAS/Buffer.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "php_rindow_openblas.h"


static zend_object_handlers rindow_openblas_blas_object_handlers;

// destractor
static void php_rindow_openblas_blas_free_object(zend_object* object)
{
    zend_object_std_dtor(object);
}

// constructor
static zend_object* php_rindow_openblas_blas_create_object(zend_class_entry* class_type) /* {{{ */
{
    zend_object* intern = NULL;

    intern = (zend_object*)ecalloc(1, sizeof(zend_object) + zend_object_properties_size(class_type));

    zend_object_std_init(intern, class_type);
    object_properties_init(intern, class_type);

    intern->handlers = &rindow_openblas_blas_object_handlers;

    return intern;
} /* }}} */

/* Method Rindow\OpenBLAS\Blas::
    public function getNumThreads() : int
 {{{ */
static PHP_METHOD(Blas, getNumThreads)
{
    int n;
    n = openblas_get_num_threads();
    RETURN_LONG(n);
}
/* }}} */

/* Method Rindow\OpenBLAS\Blas::
    public function getNumProcs() : int
 {{{ */
static PHP_METHOD(Blas, getNumProcs)
{
    int n;
    n = openblas_get_num_procs();
    RETURN_LONG(n);
}
/* }}} */

/* Method Rindow\OpenBLAS\Blas::
    public function getConfig() : string
 {{{ */
static PHP_METHOD(Blas, getConfig)
{
    char *s;
    s = openblas_get_config();
    RETURN_STRING(s);
}
/* }}} */

/* Method Rindow\OpenBLAS\Blas::
    public function getCorename() : string
 {{{ */
static PHP_METHOD(Blas, getCorename)
{
    char *s;
    s = openblas_get_corename();
    RETURN_STRING(s);
}
/* }}} */

/* Method Rindow\OpenBLAS\Blas::
    public function scal(
        int $n,
        float $alpha,
        Buffer $X, int $offsetX, int $incX) : void
 {{{ */
static PHP_METHOD(Blas, scal)
{
    php_rindow_openblas_buffer_t* buffer;
    zend_long n;
    double alpha;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;

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

    buffer = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, buffer,n,offsetX,incX)) {
        return;
    }
    switch (buffer->dtype) {
        case php_rindow_openblas_dtype_float32:
            cblas_sscal((blasint)n, (float)alpha, &(((float *)buffer->data)[offsetX]), (blasint)incX);
            break;
        case php_rindow_openblas_dtype_float64:
            cblas_dscal((blasint)n, (double)alpha, &(((double *)buffer->data)[offsetX]), (blasint)incX);
            break;
        default:
            zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
            return;
    }
}
/* }}} */

/* Method Rindow\OpenBLAS\Blas::
    public function axpy(
        int $n,
        float $alpha,
        Buffer $X, int $offsetX, int $incX,
        Buffer $Y, int $offsetY, int $incY ) : void
 {{{ */
static PHP_METHOD(Blas, axpy)
{
    php_rindow_openblas_buffer_t* bufferX;
    php_rindow_openblas_buffer_t* bufferY;
    zend_long n;
    double alpha;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zval* y=NULL;
    zend_long offsetY;
    zend_long incY;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "ldOllOll",
    //        &n,&alpha,
    //        &x,php_rindow_openblas_buffer_ce,&offsetX,&incX,
    //        &y,php_rindow_openblas_buffer_ce,&offsetY,&incY) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 8, 8)
        Z_PARAM_LONG(n)
        Z_PARAM_DOUBLE(alpha)
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
            cblas_saxpy((blasint)n, (float)alpha,
                &(((float *)bufferX->data)[offsetX]), (blasint)incX,
                &(((float *)bufferY->data)[offsetY]), (blasint)incY);
            break;
        case php_rindow_openblas_dtype_float64:
            cblas_daxpy((blasint)n, (double)alpha,
                &(((double *)bufferX->data)[offsetX]), (blasint)incX,
                &(((double *)bufferY->data)[offsetY]), (blasint)incY);
            break;
        default:
            zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
            return;
    }
}
/* }}} */

/* Method Rindow\OpenBLAS\Blas::
    public function dot(
        int $n,
        Buffer $X, int $offsetX, int $incX,
        Buffer $Y, int $offsetY, int $incY ) : float
 {{{ */
static PHP_METHOD(Blas, dot)
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
    double result;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "lOllOll",
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
            result = (double)cblas_sdot((blasint)n,
                &(((float *)bufferX->data)[offsetX]), (blasint)incX,
                &(((float *)bufferY->data)[offsetY]), (blasint)incY);
            break;
        case php_rindow_openblas_dtype_float64:
            result = (double)cblas_ddot((blasint)n,
                &(((double *)bufferX->data)[offsetX]), (blasint)incX,
                &(((double *)bufferY->data)[offsetY]), (blasint)incY);
            break;
        default:
            zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
            return;
    }
    RETURN_DOUBLE(result);
}
/* }}} */


/* Method Rindow\OpenBLAS\Blas::
    public function asum(
        int $n,
        Buffer $X, int $offsetX, int $incX ) : float
 {{{ */
static PHP_METHOD(Blas, asum)
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
        case php_rindow_openblas_dtype_float32:
            result = (double)cblas_sasum((blasint)n, &(((float *)buffer->data)[offsetX]), (blasint)incX);
            break;
        case php_rindow_openblas_dtype_float64:
            result = (double)cblas_dasum((blasint)n, &(((double *)buffer->data)[offsetX]), (blasint)incX);
            break;
        default:
            zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
            return;
    }
    RETURN_DOUBLE(result);
}
/* }}} */

/* Method Rindow\OpenBLAS\Blas::
    public function iamax(
        int $n,
        Buffer $X, int $offsetX, int $incX ) : int
 {{{ */
static PHP_METHOD(Blas, iamax)
{
    php_rindow_openblas_buffer_t* buffer;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
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
        case php_rindow_openblas_dtype_float32:
            resultIdx = (zend_long)cblas_isamax((blasint)n, &(((float *)buffer->data)[offsetX]), (blasint)incX);
            break;
        case php_rindow_openblas_dtype_float64:
            resultIdx = (zend_long)cblas_idamax((blasint)n, &(((double *)buffer->data)[offsetX]), (blasint)incX);
            break;
        default:
            zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
            return;
    }
    RETURN_LONG(resultIdx);
}
/* }}} */

#ifdef OPENBLAS_HAVE_IAMIN
/* Method Rindow\OpenBLAS\Blas::
    public function iamin(
        int $n,
        Buffer $X, int $offsetX, int $incX ) : int
 {{{ */
static PHP_METHOD(Blas, iamin)
{
    php_rindow_openblas_buffer_t* buffer;
    zend_long n;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
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
        case php_rindow_openblas_dtype_float32:
            resultIdx = (zend_long)cblas_isamin((blasint)n, &(((float *)buffer->data)[offsetX]), (blasint)incX);
            break;
        case php_rindow_openblas_dtype_float64:
            resultIdx = (zend_long)cblas_idamin((blasint)n, &(((double *)buffer->data)[offsetX]), (blasint)incX);
            break;
        default:
            zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
            return;
    }
    RETURN_LONG(resultIdx);
}
/* }}} */
#endif

/* Method Rindow\OpenBLAS\Blas::
    public function copy(
        int $n,
        Buffer $X, int $offsetX, int $incX,
        Buffer $Y, int $offsetY, int $incY ) : void
 {{{ */
static PHP_METHOD(Blas, copy)
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

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "lOllOll",
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
            cblas_scopy((blasint)n,
                &(((float *)bufferX->data)[offsetX]), (blasint)incX,
                &(((float *)bufferY->data)[offsetY]), (blasint)incY);
            break;
        case php_rindow_openblas_dtype_float64:
            cblas_dcopy((blasint)n,
                &(((double *)bufferX->data)[offsetX]), (blasint)incX,
                &(((double *)bufferY->data)[offsetY]), (blasint)incY);
            break;
        default:
            {
                zend_long i,idX,idY;
                int valueSize;
                uint8_t *x,*y;
                valueSize = php_rindow_openblas_dtype_to_valuesize(bufferX->dtype);
                x = php_rindow_openblas_get_address(bufferX,offsetX,valueSize);
                y = php_rindow_openblas_get_address(bufferY,offsetY,valueSize);
                if(incX==1 && incY==1) {
                    memcpy(y,x,valueSize*n);
                } else {
                    for(i=0,idX=0,idY=0; i<n; i++,idX+=incX,idY+=incY) {
                        memcpy(&y[valueSize*idY],&x[valueSize*idX],valueSize);
                    }
                }
            }
            break;
    }
}
/* }}} */

/* Method Rindow\OpenBLAS\Blas::
    public function nrm2(
        int $n,
        Buffer $X, int $offsetX, int $incX ) : float
 {{{ */
static PHP_METHOD(Blas, nrm2)
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
        case php_rindow_openblas_dtype_float32:
            result = (double)cblas_snrm2((blasint)n, &(((float *)buffer->data)[offsetX]), (blasint)incX);
            break;
        case php_rindow_openblas_dtype_float64:
            result = (double)cblas_dnrm2((blasint)n, &(((double *)buffer->data)[offsetX]), (blasint)incX);
            break;
        default:
            zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
            return;
    }
    RETURN_DOUBLE(result);
}
/* }}} */

/* Method Rindow\OpenBLAS\Blas::
    public function rotg(
        Buffer $A, int $offsetA,
        Buffer $B, int $offsetB,
        Buffer $C, int $offsetC,
        Buffer $S, int $offsetS,
    ) : void
 {{{ */
static PHP_METHOD(Blas, rotg)
{
    zval* a;
    zend_long offsetA;
    zval* b;
    zend_long offsetB;
    zval* c;
    zend_long offsetC;
    zval* s;
    zend_long offsetS;
    php_rindow_openblas_buffer_t* bufferA;
    php_rindow_openblas_buffer_t* bufferB;
    php_rindow_openblas_buffer_t* bufferC;
    php_rindow_openblas_buffer_t* bufferS;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 8, 8)
        Z_PARAM_OBJECT_OF_CLASS(a,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetA)
        Z_PARAM_OBJECT_OF_CLASS(b,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetB)
        Z_PARAM_OBJECT_OF_CLASS(c,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetC)
        Z_PARAM_OBJECT_OF_CLASS(s,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetS)
    ZEND_PARSE_PARAMETERS_END();

    // Check Buffer A
    bufferA = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, bufferA,1,offsetA,1)) {
        return;
    }

    // Check Buffer B
    bufferB = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(b);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_Y, bufferB,1,offsetB,1)) {
        return;
    }

    // Check Buffer C
    bufferC = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(c);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, bufferC,1,offsetC,1)) {
        return;
    }

    // Check Buffer S
    bufferS = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(s);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_Y, bufferS,1,offsetS,1)) {
        return;
    }


    // Check Buffer A and B and C and S
    if(bufferA->dtype!=bufferB->dtype || bufferB->dtype!=bufferC->dtype ||
        bufferC->dtype!=bufferS->dtype || bufferS->dtype!=bufferA->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for A,B,C and S", 0);
        return;
    }

    switch (bufferA->dtype) {
        case php_rindow_openblas_dtype_float32:
            cblas_srotg(
                &(((float *)bufferA->data)[offsetA]),
                &(((float *)bufferB->data)[offsetB]),
                &(((float *)bufferC->data)[offsetC]),
                &(((float *)bufferS->data)[offsetS])
            );
            break;
        case php_rindow_openblas_dtype_float64:
            cblas_drotg(
                &(((double *)bufferA->data)[offsetA]),
                &(((double *)bufferB->data)[offsetB]),
                &(((double *)bufferC->data)[offsetC]),
                &(((double *)bufferS->data)[offsetS])
            );
            break;
        default:
            zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
            return;
    }
}
/* }}} */

/* Method Rindow\OpenBLAS\Blas::
    public function rot(
        int $n,
        Buffer $X, int $offsetX, int $incX,
        Buffer $Y, int $offsetY, int $incY,
        Buffer $C, int $offsetC,
        Buffer $S, int $offsetS,
    ) : void
 {{{ */
static PHP_METHOD(Blas, rot)
{
    zend_long n;
    zval* x;
    zend_long offsetX;
    zend_long incX;
    zval* y;
    zend_long offsetY;
    zend_long incY;
    zval* c;
    zend_long offsetC;
    zval* s;
    zend_long offsetS;
    php_rindow_openblas_buffer_t* bufferX;
    php_rindow_openblas_buffer_t* bufferY;
    php_rindow_openblas_buffer_t* bufferC;
    php_rindow_openblas_buffer_t* bufferS;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 11, 11)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_OBJECT_OF_CLASS(y,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetY)
        Z_PARAM_LONG(incY)
        Z_PARAM_OBJECT_OF_CLASS(c,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetC)
        Z_PARAM_OBJECT_OF_CLASS(s,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetS)
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

    // Check Buffer C
    bufferC = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(c);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, bufferC,1,offsetC,1)) {
        return;
    }

    // Check Buffer S
    bufferS = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(s);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_Y, bufferS,1,offsetS,1)) {
        return;
    }

    // Check Buffer X and Y
    if(bufferX->dtype!=bufferY->dtype||bufferY->dtype!=bufferC->dtype||
        bufferC->dtype!=bufferS->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for X,Y,C and S", 0);
        return;
    }

    switch (bufferX->dtype) {
        case php_rindow_openblas_dtype_float32:
            cblas_srot((blasint)n,
                &(((float *)bufferX->data)[offsetX]), (blasint)incX,
                &(((float *)bufferY->data)[offsetY]), (blasint)incY,
                ((float *)bufferC->data)[offsetC],
                ((float *)bufferS->data)[offsetS]
            );
            break;
        case php_rindow_openblas_dtype_float64:
            cblas_drot((blasint)n,
                &(((double *)bufferX->data)[offsetX]), (blasint)incX,
                &(((double *)bufferY->data)[offsetY]), (blasint)incY,
                ((double *)bufferC->data)[offsetC],
                ((double *)bufferS->data)[offsetS]
            );
            break;
        default:
            zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
            return;
    }
}
/* }}} */

/* Method Rindow\OpenBLAS\Blas::
    public function swap(
        int $n,
        Buffer $X, int $offsetX, int $incX,
        Buffer $Y, int $offsetY, int $incY ) : void
 {{{ */
static PHP_METHOD(Blas, swap)
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
            cblas_sswap((blasint)n,
                &(((float *)bufferX->data)[offsetX]), (blasint)incX,
                &(((float *)bufferY->data)[offsetY]), (blasint)incY);
            break;
        case php_rindow_openblas_dtype_float64:
            cblas_dswap((blasint)n,
                &(((double *)bufferX->data)[offsetX]), (blasint)incX,
                &(((double *)bufferY->data)[offsetY]), (blasint)incY);
            break;
        default:
            zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
            return;
    }
}
/* }}} */

/* Method Rindow\OpenBLAS\Blas::
    public function gemv(
        int $order,
        int $trans,
        int $m,
        int $n,
        float $alpha,
        Buffer $A, int $offsetA, int $ldA,
        Buffer $X, int $offsetX, int $incX,
        float $beta,
        Buffer $Y, int $offsetY, int $incY ) : void
 {{{ */
static PHP_METHOD(Blas, gemv)
{
    php_rindow_openblas_buffer_t* bufferA;
    php_rindow_openblas_buffer_t* bufferX;
    php_rindow_openblas_buffer_t* bufferY;
    zend_long order;
    zend_long trans;
    zend_long m;
    zend_long n;
    double alpha;
    zval* a=NULL;
    zend_long offsetA;
    zend_long ldA;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    double beta;
    zval* y=NULL;
    zend_long offsetY;
    zend_long incY;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "lllldOllOlldOll",
    //        &order,&trans,&m,&n,&alpha,
    //        &a,php_rindow_openblas_buffer_ce,&offsetA,&ldA,
    //        &x,php_rindow_openblas_buffer_ce,&offsetX,&incX,
    //        &beta,
    //        &y,php_rindow_openblas_buffer_ce,&offsetY,&incY) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 15, 15)
        Z_PARAM_LONG(order)
        Z_PARAM_LONG(trans)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_DOUBLE(alpha)
        Z_PARAM_OBJECT_OF_CLASS(a,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_DOUBLE(beta)
        Z_PARAM_OBJECT_OF_CLASS(y,php_rindow_openblas_buffer_ce)
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
    bufferA = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(a);
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_A, bufferA,m,n,offsetA,ldA)) {
        return;
    }

    // Check Buffer X
    bufferX = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);

    // Check Buffer Y
    bufferY = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(y);

    // Check Buffer size X and Y
    {
        zend_long rows,cols;
        if(trans==CblasNoTrans || trans==CblasConjNoTrans ) {
            rows = m; cols = n;
        } else if(trans==CblasTrans || trans==CblasConjTrans) {
            rows = n; cols = m;
        } else {
            zend_throw_exception(spl_ce_RuntimeException, "unknown transpose mode for bufferA.", 0);
            return;
        }
        if(php_rindow_openblas_assert_vector_buffer_spec(
            PHP_RINDOW_OPENBLAS_ASSERT_X, bufferX,cols,offsetX,incX)) {
            return;
        }
        if(php_rindow_openblas_assert_vector_buffer_spec(
            PHP_RINDOW_OPENBLAS_ASSERT_Y, bufferY,rows,offsetY,incY)) {
            return;
        }
    }

    // Check Buffer A and X and Y
    if(bufferA->dtype!=bufferX->dtype || bufferX->dtype!=bufferY->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for A and X and Y", 0);
        return;
    }

    switch (bufferX->dtype) {
        case php_rindow_openblas_dtype_float32:
            cblas_sgemv(
                (OPENBLAS_CONST enum CBLAS_ORDER)order,
                (OPENBLAS_CONST enum CBLAS_TRANSPOSE)trans,
                (blasint)m,(blasint)n,
                (float)alpha,
                &(((float *)bufferA->data)[offsetA]), (blasint)ldA,
                &(((float *)bufferX->data)[offsetX]), (blasint)incX,
                (float)beta,
                &(((float *)bufferY->data)[offsetY]), (blasint)incY);
            break;
        case php_rindow_openblas_dtype_float64:
            cblas_dgemv(
                (OPENBLAS_CONST enum CBLAS_ORDER)order,
                (OPENBLAS_CONST enum CBLAS_TRANSPOSE)trans,
                (blasint)m,(blasint)n,
                (double)alpha,
                &(((double *)bufferA->data)[offsetA]), (blasint)ldA,
                &(((double *)bufferX->data)[offsetX]), (blasint)incX,
                (double)beta,
                &(((double *)bufferY->data)[offsetY]), (blasint)incY);
            break;
        default:
            zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
            return;
    }
}
/* }}} */

/* Method Rindow\OpenBLAS\Blas::
    public function gemm(
        int $order,
        int $transA,
        int $transB,
        int $m,
        int $n,
        int $k,
        float $alpha,
        Buffer $A, int $offsetA, int $ldA,
        Buffer $B, int $offsetB, int $ldB,
        float $beta,
        Buffer $C, int $offsetC, int $ldC ) : void
 {{{ */
static PHP_METHOD(Blas, gemm)
{
    php_rindow_openblas_buffer_t* bufferA;
    php_rindow_openblas_buffer_t* bufferB;
    php_rindow_openblas_buffer_t* bufferC;
    zend_long order;
    zend_long transA;
    zend_long transB;
    zend_long m;
    zend_long n;
    zend_long k;
    double alpha;
    zval* a=NULL;
    zend_long offsetA;
    zend_long ldA;
    zval* b=NULL;
    zend_long offsetB;
    zend_long ldB;
    double beta;
    zval* c=NULL;
    zend_long offsetC;
    zend_long ldC;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "lllllldOllOlldOll",
    //        &order,&transA,&transB,&m,&n,&k,&alpha,
    //        &a,php_rindow_openblas_buffer_ce,&offsetA,&ldA,
    //        &b,php_rindow_openblas_buffer_ce,&offsetB,&ldB,
    //        &beta,
    //        &c,php_rindow_openblas_buffer_ce,&offsetC,&ldC) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 17, 17)
        Z_PARAM_LONG(order)
        Z_PARAM_LONG(transA)
        Z_PARAM_LONG(transB)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_LONG(k)
        Z_PARAM_DOUBLE(alpha)
        Z_PARAM_OBJECT_OF_CLASS(a,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)
        Z_PARAM_OBJECT_OF_CLASS(b,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetB)
        Z_PARAM_LONG(ldB)
        Z_PARAM_DOUBLE(beta)
        Z_PARAM_OBJECT_OF_CLASS(c,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetC)
        Z_PARAM_LONG(ldC)
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
    bufferA = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(a);
    {
        zend_long rows,cols;
        if(transA==CblasNoTrans || transA==CblasConjNoTrans) {
            rows = m; cols = k;
        } else if(transA==CblasTrans || transA==CblasConjTrans) {
            rows = k; cols = m;
        } else {
            zend_throw_exception(spl_ce_RuntimeException, "unknown transpose mode for bufferA.", 0);
            return;
        }
        if(php_rindow_openblas_assert_matrix_buffer_spec(
            PHP_RINDOW_OPENBLAS_ASSERT_A, bufferA,rows,cols,offsetA,ldA)) {
            return;
        }
    }

    // Check Buffer B
    bufferB = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(b);
    {
        zend_long rows,cols;
        if(transB==CblasNoTrans || transB==CblasConjNoTrans) {
            rows = k; cols = n;
        } else if(transB==CblasTrans || transB==CblasConjTrans) {
            rows = n; cols = k;
        } else {
            zend_throw_exception(spl_ce_RuntimeException, "unknown transpose mode for bufferB.", 0);
            return;
        }
        if(php_rindow_openblas_assert_matrix_buffer_spec(
            PHP_RINDOW_OPENBLAS_ASSERT_B, bufferB,rows,cols,offsetB,ldB)) {
            return;
        }
    }

    // Check Buffer C
    bufferC = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(c);
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_C, bufferC,m,n,offsetC,ldC)) {
        return;
    }

    // Check Buffer A and B and C
    if(bufferA->dtype!=bufferB->dtype || bufferB->dtype!=bufferC->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for A and B and C", 0);
        return;
    }

    switch (bufferA->dtype) {
        case php_rindow_openblas_dtype_float32:
            cblas_sgemm(
                (OPENBLAS_CONST enum CBLAS_ORDER)order,
                (OPENBLAS_CONST enum CBLAS_TRANSPOSE)transA,
                (OPENBLAS_CONST enum CBLAS_TRANSPOSE)transB,
                (blasint)m,(blasint)n,(blasint)k,
                (float)alpha,
                &(((float *)bufferA->data)[offsetA]), (blasint)ldA,
                &(((float *)bufferB->data)[offsetB]), (blasint)ldB,
                (float)beta,
                &(((float *)bufferC->data)[offsetC]), (blasint)ldC);
            break;
        case php_rindow_openblas_dtype_float64:
            cblas_dgemm(
                (OPENBLAS_CONST enum CBLAS_ORDER)order,
                (OPENBLAS_CONST enum CBLAS_TRANSPOSE)transA,
                (OPENBLAS_CONST enum CBLAS_TRANSPOSE)transB,
                (blasint)m,(blasint)n,(blasint)k,
                (double)alpha,
                &(((double *)bufferA->data)[offsetA]), (blasint)ldA,
                &(((double *)bufferB->data)[offsetB]), (blasint)ldB,
                (double)beta,
                &(((double *)bufferC->data)[offsetC]), (blasint)ldC);
            break;
        default:
            zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
            return;
    }
}
/* }}} */



ZEND_BEGIN_ARG_INFO_EX(ai_Blas_scal, 0, 0, 5)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, alpha)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Blas_axpy, 0, 0, 8)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, alpha)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_OBJ_INFO(0, y, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, incY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Blas_dot, 0, 0, 7)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_OBJ_INFO(0, y, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, incY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Blas_asum, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Blas_iamax, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

#ifdef OPENBLAS_HAVE_IAMIN
ZEND_BEGIN_ARG_INFO_EX(ai_Blas_iamin, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()
#endif

ZEND_BEGIN_ARG_INFO_EX(ai_Blas_copy, 0, 0, 7)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_OBJ_INFO(0, y, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, incY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Blas_nrm2, 0, 0, 4)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Blas_rotg, 0, 0, 8)
    ZEND_ARG_OBJ_INFO(0, a, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_OBJ_INFO(0, b, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetB)
    ZEND_ARG_OBJ_INFO(0, c, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetC)
    ZEND_ARG_OBJ_INFO(0, s, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetS)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Blas_rot, 0, 0, 11)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_OBJ_INFO(0, y, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, incY)
    ZEND_ARG_OBJ_INFO(0, c, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetC)
    ZEND_ARG_OBJ_INFO(0, s, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetS)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Blas_swap, 0, 0, 7)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_OBJ_INFO(0, y, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, incY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Blas_gemv, 0, 0, 15)
    ZEND_ARG_INFO(0, order)
    ZEND_ARG_INFO(0, trans)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, alpha)
    ZEND_ARG_OBJ_INFO(0, a, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_INFO(0, ldA)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetX)
    ZEND_ARG_INFO(0, incX)
    ZEND_ARG_INFO(0, beta)
    ZEND_ARG_OBJ_INFO(0, x, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetY)
    ZEND_ARG_INFO(0, incY)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Blas_gemm, 0, 0, 17)
    ZEND_ARG_INFO(0, order)
    ZEND_ARG_INFO(0, transA)
    ZEND_ARG_INFO(0, transB)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)
    ZEND_ARG_INFO(0, k)
    ZEND_ARG_INFO(0, alpha)
    ZEND_ARG_OBJ_INFO(0, a, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_INFO(0, ldA)
    ZEND_ARG_OBJ_INFO(0, b, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetB)
    ZEND_ARG_INFO(0, ldB)
    ZEND_ARG_INFO(0, beta)
    ZEND_ARG_OBJ_INFO(0, c, Rindow\\OpenBLAS\\Buffer, 0)
    ZEND_ARG_INFO(0, offsetC)
    ZEND_ARG_INFO(0, ldC)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Blas_void, 0, 0, 0)
ZEND_END_ARG_INFO()

/* {{{ Rindow\OpenBLAS\Blas function entries */
static zend_function_entry php_rindow_openblas_blas_me[] = {
    /* clang-format off */
    PHP_ME(Blas, getNumThreads, ai_Blas_void, ZEND_ACC_PUBLIC)
    PHP_ME(Blas, getNumProcs,   ai_Blas_void, ZEND_ACC_PUBLIC)
    PHP_ME(Blas, getConfig,     ai_Blas_void, ZEND_ACC_PUBLIC)
    PHP_ME(Blas, getCorename,   ai_Blas_void, ZEND_ACC_PUBLIC)
    PHP_ME(Blas, scal,  ai_Blas_scal,  ZEND_ACC_PUBLIC)
    PHP_ME(Blas, axpy,  ai_Blas_axpy,  ZEND_ACC_PUBLIC)
    PHP_ME(Blas, dot,   ai_Blas_dot,   ZEND_ACC_PUBLIC)
    PHP_ME(Blas, asum,  ai_Blas_asum,  ZEND_ACC_PUBLIC)
    PHP_ME(Blas, iamax, ai_Blas_iamax, ZEND_ACC_PUBLIC)
#ifdef OPENBLAS_HAVE_IAMIN
    PHP_ME(Blas, iamin, ai_Blas_iamin, ZEND_ACC_PUBLIC)
#endif
    PHP_ME(Blas, copy,  ai_Blas_copy,  ZEND_ACC_PUBLIC)
    PHP_ME(Blas, nrm2,  ai_Blas_nrm2,  ZEND_ACC_PUBLIC)
    PHP_ME(Blas, rotg,  ai_Blas_rotg,  ZEND_ACC_PUBLIC)
    PHP_ME(Blas, rot,   ai_Blas_rot,   ZEND_ACC_PUBLIC)
    PHP_ME(Blas, swap,  ai_Blas_swap,  ZEND_ACC_PUBLIC)
    PHP_ME(Blas, gemv,  ai_Blas_gemv,  ZEND_ACC_PUBLIC)
    PHP_ME(Blas, gemm,  ai_Blas_gemm,  ZEND_ACC_PUBLIC)
    PHP_FE_END
    /* clang-format on */
};
/* }}} */

/* Class Rindow\OpenBLAS\Blas {{{ */
static zend_class_entry* rindow_openblas_blas_ce;

void php_rindow_openblas_blas_init_ce(INIT_FUNC_ARGS)
{
    zend_class_entry ce;

    INIT_NS_CLASS_ENTRY(ce, "Rindow\\OpenBLAS", "Blas", php_rindow_openblas_blas_me);
    rindow_openblas_blas_ce = zend_register_internal_class(&ce);
    rindow_openblas_blas_ce->create_object = php_rindow_openblas_blas_create_object;

    memcpy(&rindow_openblas_blas_object_handlers, zend_get_std_object_handlers(), sizeof(zend_object_handlers));
    rindow_openblas_blas_object_handlers.offset    = 0;
    rindow_openblas_blas_object_handlers.free_obj  = php_rindow_openblas_blas_free_object;
    rindow_openblas_blas_object_handlers.clone_obj = NULL;

    //zend_class_implements(rindow_openblas_blas_ce, 2, spl_ce_ArrayAccess, spl_ce_Countable);
}
/* }}} */
