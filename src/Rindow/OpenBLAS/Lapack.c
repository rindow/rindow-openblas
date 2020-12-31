#include <php.h>
#include <Zend/zend_interfaces.h>
#include <Zend/zend_exceptions.h>
#include <ext/spl/spl_iterators.h>
#include <ext/spl/spl_exceptions.h>
#if _MSC_VER
#include <complex.h>
#define lapack_complex_float _Fcomplex
#define lapack_complex_double _Dcomplex
#endif
#include <lapacke.h>
#include <Interop/Polite/Math/Matrix.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "php_rindow_openblas.h"

#ifndef MAX
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif

static zend_object_handlers rindow_openblas_lapack_object_handlers;

// destractor
static void php_rindow_openblas_lapack_free_object(zend_object* object)
{
    zend_object_std_dtor(object);
}

// constructor
static zend_object* php_rindow_openblas_lapack_create_object(zend_class_entry* class_type) /* {{{ */
{
    zend_object* intern = NULL;

    intern = (zend_object*)ecalloc(1, sizeof(zend_object) + zend_object_properties_size(class_type));

    zend_object_std_init(intern, class_type);
    object_properties_init(intern, class_type);

    intern->handlers = &rindow_openblas_lapack_object_handlers;

    return intern;
} /* }}} */


/* Method Rindow\OpenBLAS\Lapack::
    public function gesvd(
        int $matrix_layout,
        string $jobu,
        string $jobvt,
        int $m,
        int $n,
        Buffer $A,  int $offsetA,  int $ldA,
        Buffer $S,  int $offsetS,
        Buffer $U,  int $offsetU,  int $ldU,
        Buffer $VT, int $offsetVT, int $ldVT,
        Buffer $SuperB,  int $offsetSuperB
     ) : void
 {{{ */
static PHP_METHOD(Lapack, gesvd)
{
/*
lapack_int LAPACKE_dgesvd( int matrix_layout, char jobu, char jobvt,
                           lapack_int m, lapack_int n, double* a,
                           lapack_int lda, double* s, double* u, lapack_int ldu,
                           double* vt, lapack_int ldvt, double* superb );
*/
    php_interop_polite_math_matrix_linear_buffer_t* bufferA;
    php_interop_polite_math_matrix_linear_buffer_t* bufferS;
    php_interop_polite_math_matrix_linear_buffer_t* bufferU;
    php_interop_polite_math_matrix_linear_buffer_t* bufferVT;
    php_interop_polite_math_matrix_linear_buffer_t* bufferSuperB;
    zend_long matrix_layout;
    zend_long jobu;
    zend_long jobvt;
    zend_long m;
    zend_long n;
    zval* objA=NULL;
    zend_long offsetA;
    zend_long ldA;
    zval* objS=NULL;
    zend_long offsetS;
    zval* objU=NULL;
    zend_long offsetU;
    zend_long ldU;
    zval* objVT=NULL;
    zend_long offsetVT;
    zend_long ldVT;
    zval* objSuperB=NULL;
    zend_long offsetSuperB;
    int info;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 18, 18)
        Z_PARAM_LONG(matrix_layout)
        Z_PARAM_LONG(jobu)
        Z_PARAM_LONG(jobvt)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)

        Z_PARAM_OBJECT(objA) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)
        Z_PARAM_OBJECT(objS) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetS)

        Z_PARAM_OBJECT(objU) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetU)
        Z_PARAM_LONG(ldU)
        Z_PARAM_OBJECT(objVT) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetVT)

        Z_PARAM_LONG(ldVT)
        Z_PARAM_OBJECT(objSuperB) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(offsetSuperB)
    ZEND_PARSE_PARAMETERS_END();

    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_M, m)) {
        return;
    }
    if(php_rindow_openblas_assert_shape_parameter(
        PHP_RINDOW_OPENBLAS_ASSERT_N, n)) {
        return;
    }
    if( offsetS < 0 ) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "offsetS must be greater than zero or equal", 0);
        return;
    }
    if( offsetU < 0 ) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "offsetU must be greater than zero or equal", 0);
        return;
    }
    if( ldU <= 0 ) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "ldU must be greater than zero", 0);
        return;
    }
    if( offsetVT < 0 ) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "offsetVT must be greater than zero or equal", 0);
        return;
    }
    if( ldVT <= 0 ) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "ldVT must be greater than zero", 0);
        return;
    }
    if( offsetSuperB < 0 ) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "offsetVT must be greater than zero or equal", 0);
        return;
    }
    // Check Buffer A
    bufferA = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(objA);
    if(php_rindow_openblas_assert_buffer_type(bufferA,"a")) {
        return;
    }
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_A, bufferA,m,n,offsetA,ldA)) {
        return;
    }

    // Check Buffer S
    bufferS = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(objS);
    if(php_rindow_openblas_assert_buffer_type(bufferS,"s")) {
        return;
    }
    if( offsetS+MIN(m,n) > bufferS->size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "BufferS size is too small", 0);
        return;
    }

    // Check Buffer U
    bufferU = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(objU);
    if(php_rindow_openblas_assert_buffer_type(bufferU,"u")) {
        return;
    }
    if( offsetU+m*ldU > bufferU->size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "BufferU size is too small", 0);
        return;
    }

    // Check Buffer VT
    bufferVT = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(objVT);
    if(php_rindow_openblas_assert_buffer_type(bufferVT,"vt")) {
        return;
    }
    if( offsetVT+ldVT*n > bufferVT->size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "BufferVT size is too small", 0);
        return;
    }

    // Check Buffer SuperB
    bufferSuperB = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(objSuperB);
    if(php_rindow_openblas_assert_buffer_type(bufferSuperB,"b")) {
        return;
    }
    if( offsetSuperB+MIN(m,n)-1 > bufferSuperB->size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "bufferSuperB size is too small", 0);
        return;
    }

    // Check Buffer A and B and C
    if(bufferA->dtype!=bufferS->dtype ||
       bufferA->dtype!=bufferU->dtype ||
       bufferA->dtype!=bufferVT->dtype ||
       bufferA->dtype!=bufferSuperB->dtype
    ) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type", 0);
        return;
    }

    switch (bufferA->dtype) {
        case php_interop_polite_math_matrix_dtype_float32:
            info = LAPACKE_sgesvd(
                (int)matrix_layout,
                (char)jobu,
                (char)jobvt,
                (lapack_int)m,(lapack_int)n,
                &(((float *)bufferA->data)[offsetA]), (lapack_int)ldA,
                &(((float *)bufferS->data)[offsetS]),
                &(((float *)bufferU->data)[offsetU]), (lapack_int)ldU,
                &(((float *)bufferVT->data)[offsetVT]), (lapack_int)ldVT,
                &(((float *)bufferSuperB->data)[offsetSuperB])
            );
            break;
        case php_interop_polite_math_matrix_dtype_float64:
            info = LAPACKE_dgesvd(
                (int)matrix_layout,
                (char)jobu,
                (char)jobvt,
                (lapack_int)m,(lapack_int)n,
                &(((double *)bufferA->data)[offsetA]), (lapack_int)ldA,
                &(((double *)bufferS->data)[offsetS]),
                &(((double *)bufferU->data)[offsetU]), (lapack_int)ldU,
                &(((double *)bufferVT->data)[offsetVT]), (lapack_int)ldVT,
                &(((double *)bufferSuperB->data)[offsetSuperB])
            );
            break;
        default:
            zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
            return;
    }
    if( info == LAPACK_WORK_MEMORY_ERROR ) {
        zend_throw_exception(spl_ce_RuntimeException, "Not enough memory to allocate work array.", info);
        return;
    } else if( info == LAPACK_TRANSPOSE_MEMORY_ERROR ) {
        zend_throw_exception(spl_ce_RuntimeException, "Not enough memory to transpose matrix.", info);
        return;
    } else if( info < 0 ) {
        zend_throw_exception(spl_ce_RuntimeException, "Wrong parameter.", info);
        return;
    }
}
/* }}} */


ZEND_BEGIN_ARG_INFO_EX(ai_Lapack_gesvd, 0, 0, 18)
    ZEND_ARG_INFO(0, matrix_layout)
    ZEND_ARG_INFO(0, jobu)
    ZEND_ARG_INFO(0, jobvt)
    ZEND_ARG_INFO(0, m)
    ZEND_ARG_INFO(0, n)

    ZEND_ARG_OBJ_INFO(0, objA, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetA)
    ZEND_ARG_INFO(0, ldA)
    ZEND_ARG_OBJ_INFO(0, objS, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetS)

    ZEND_ARG_OBJ_INFO(0, objU, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetU)
    ZEND_ARG_INFO(0, ldU)
    ZEND_ARG_OBJ_INFO(0, objVT, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetVT)

    ZEND_ARG_INFO(0, ldVT)
    ZEND_ARG_OBJ_INFO(0, objSuperB, Interop\\Polite\\Math\\Matrix\\LinearBuffer, 0)
    ZEND_ARG_INFO(0, offsetSuperB)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Lapack_void, 0, 0, 0)
ZEND_END_ARG_INFO()

/* {{{ Rindow\OpenBLAS\Lapack function entries */
static zend_function_entry php_rindow_openblas_lapack_me[] = {
    /* clang-format off */
    PHP_ME(Lapack, gesvd,  ai_Lapack_gesvd,  ZEND_ACC_PUBLIC)
    PHP_FE_END
    /* clang-format on */
};
/* }}} */

/* Class Rindow\OpenBLAS\Lapack {{{ */
static zend_class_entry* rindow_openblas_lapack_ce;

void php_rindow_openblas_lapack_init_ce(INIT_FUNC_ARGS)
{
    zend_class_entry ce;

    INIT_NS_CLASS_ENTRY(ce, "Rindow\\OpenBLAS", "Lapack", php_rindow_openblas_lapack_me);
    rindow_openblas_lapack_ce = zend_register_internal_class(&ce);
    rindow_openblas_lapack_ce->create_object = php_rindow_openblas_lapack_create_object;

    memcpy(&rindow_openblas_lapack_object_handlers, zend_get_std_object_handlers(), sizeof(zend_object_handlers));
    rindow_openblas_lapack_object_handlers.offset    = 0;
    rindow_openblas_lapack_object_handlers.free_obj  = php_rindow_openblas_lapack_free_object;
    rindow_openblas_lapack_object_handlers.clone_obj = NULL;

    //zend_class_implements(rindow_openblas_lapack_ce, 2, spl_ce_ArrayAccess, spl_ce_Countable);
}
/* }}} */
