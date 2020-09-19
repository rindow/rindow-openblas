/* rindow_openblas extension for PHP */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <php.h>
#include <Zend/zend_interfaces.h>
#include <Zend/zend_exceptions.h>
#include <ext/spl/spl_iterators.h>
#include <ext/spl/spl_exceptions.h>
#include "ext/standard/info.h"
#include <Rindow/OpenBLAS/Buffer.h>
#include "php_rindow_openblas.h"

/* For compatibility with older PHP versions */
#ifndef ZEND_PARSE_PARAMETERS_NONE
#define ZEND_PARSE_PARAMETERS_NONE() \
    ZEND_PARSE_PARAMETERS_START(0, 0) \
    ZEND_PARSE_PARAMETERS_END()
#endif

int php_rindow_openblas_dtype_to_valuesize(zend_long dtype)
{
    switch (dtype) {
        case php_rindow_openblas_dtype_bool:
        case php_rindow_openblas_dtype_int8:
        case php_rindow_openblas_dtype_uint8:
        case php_rindow_openblas_dtype_float8:
            return 1;
        case php_rindow_openblas_dtype_int16:
        case php_rindow_openblas_dtype_uint16:
        case php_rindow_openblas_dtype_float16:
            return 2;
        case php_rindow_openblas_dtype_int32:
        case php_rindow_openblas_dtype_uint32:
        case php_rindow_openblas_dtype_float32:
            return 4;
        case php_rindow_openblas_dtype_int64:
        case php_rindow_openblas_dtype_uint64:
        case php_rindow_openblas_dtype_float64:
            return 8;
    }
    return 0;
}


int php_rindow_openblas_dtype_is_int(zend_long dtype)
{
    switch (dtype) {
        case php_rindow_openblas_dtype_int8:
        case php_rindow_openblas_dtype_uint8:
        case php_rindow_openblas_dtype_int16:
        case php_rindow_openblas_dtype_uint16:
        case php_rindow_openblas_dtype_int32:
        case php_rindow_openblas_dtype_uint32:
        case php_rindow_openblas_dtype_int64:
        case php_rindow_openblas_dtype_uint64:
            return 1;
    }
    return 0;
}

int php_rindow_openblas_dtype_is_float(zend_long dtype)
{
    switch (dtype) {
        case php_rindow_openblas_dtype_float8:
        case php_rindow_openblas_dtype_float16:
        case php_rindow_openblas_dtype_float32:
        case php_rindow_openblas_dtype_float64:
            return 1;
    }
    return 0;
}

int php_rindow_openblas_dtype_is_bool(zend_long dtype)
{
    switch (dtype) {
        case php_rindow_openblas_dtype_bool:
            return 1;
    }
    return 0;
}

int php_rindow_openblas_assert_shape_parameter(
    int name, zend_long n)
{
    static const char *message[3] = {
        "Argument m must be greater than 0.",
        "Argument n must be greater than 0.",
        "Argument k must be greater than 0.",
    };
    if(n<1) {
        zend_throw_exception(spl_ce_RuntimeException, message[name], 0);
        return -1;
    }
    return 0;
}

int php_rindow_openblas_assert_vector_buffer_spec(
    int name,php_rindow_openblas_buffer_t *buffer,
    zend_long n, zend_long offset, zend_long inc)
{
    static const char *message[2][4] = {
        {
            "uninitialized array",
            "Argument offsetX must be greater than equals 0.",
            "Argument incX must be greater than 0.",
            "Vector specification too large for bufferX."
        },
        {
            "uninitialized array",
            "Argument offsetY must be greater than equals 0.",
            "Argument incY must be greater than 0.",
            "Vector specification too large for bufferY."
        },
    };
    if(buffer->data==NULL) {
        zend_throw_exception(spl_ce_DomainException, message[name][0], 0);
        return -1;
    }
    if(offset<0) {
        zend_throw_exception(spl_ce_RuntimeException, message[name][1], 0);
        return -1;
    }
    if(inc<1) {
        zend_throw_exception(spl_ce_RuntimeException, message[name][2], 0);
        return -1;
    }
    if(offset+(n-1)*inc >= buffer->size) {
        zend_throw_exception(spl_ce_RuntimeException, message[name][3], 0);
        return -1;
    }

    return 0;
}

int php_rindow_openblas_assert_matrix_buffer_spec(
    int name, php_rindow_openblas_buffer_t *buffer,
    zend_long m,zend_long n, zend_long offset, zend_long ld)
{
    static const char *message[4][4] = {
        {
            "uninitialized array",
            "Argument offsetA must be greater than equals 0.",
            "Argument ldA must be greater than 0.",
            "Matrix specification too large for bufferA."
        },
        {
            "uninitialized array",
            "Argument offsetY must be greater than equals 0.",
            "Argument ldY must be greater than 0.",
            "Matrix specification too large for bufferY."
        },
        {
            "uninitialized array",
            "Argument offsetB must be greater than equals 0.",
            "Argument ldB must be greater than 0.",
            "Matrix specification too large for bufferB."
        },
        {
            "uninitialized array",
            "Argument offsetC must be greater than equals 0.",
            "Argument ldC must be greater than 0.",
            "Matrix specification too large for bufferC."
        },
    };
    if(buffer->data==NULL) {
        zend_throw_exception(spl_ce_DomainException, message[name][0], 0);
        return -1;
    }
    if(offset<0) {
        zend_throw_exception(spl_ce_RuntimeException, message[name][1], 0);
        return -1;
    }
    if(ld<1) {
        zend_throw_exception(spl_ce_RuntimeException, message[name][2], 0);
        return -1;
    }
    if(offset+(m-1)*ld+(n-1) >= buffer->size) {
        zend_throw_exception(spl_ce_RuntimeException, message[name][3], 0);
        return -1;
    }

    return 0;
}

int php_rindow_openblas_assert_buffer_size(
    php_rindow_openblas_buffer_t *buffer,
    zend_long offset,zend_long size,
    char* message)
{
    if(size<1 || offset<0 ||
        buffer->size < offset+size) {
        zend_throw_exception(spl_ce_RuntimeException, message, 0);
        return -1;
    }
    return 0;
}

/* {{{ PHP_RINIT_FUNCTION
 */
PHP_RINIT_FUNCTION(rindow_openblas)
{
#if defined(ZTS) && defined(COMPILE_DL_RINDOW_OPENBLAS)
    ZEND_TSRMLS_CACHE_UPDATE();
#endif

    return SUCCESS;
}
/* }}} */

/* {{{ PHP_MINFO_FUNCTION
 */
PHP_MINFO_FUNCTION(rindow_openblas)
{
    php_info_print_table_start();
    php_info_print_table_header(2, "rindow_openblas support", "enabled");
    php_info_print_table_end();
}
/* }}} */

/* {{{ rindow_openblas_functions[]
 */
//static const zend_function_entry rindow_openblas_functions[] = {
//};
/* }}} */

PHP_MINIT_FUNCTION(rindow_openblas)
{
    php_rindow_openblas_buffer_init_ce(INIT_FUNC_ARGS_PASSTHRU);
    php_rindow_openblas_blas_init_ce(INIT_FUNC_ARGS_PASSTHRU);
    php_rindow_openblas_lapack_init_ce(INIT_FUNC_ARGS_PASSTHRU);
    php_rindow_openblas_math_init_ce(INIT_FUNC_ARGS_PASSTHRU);
    return SUCCESS;
}

/* {{{ php_rindow_openblas_module_entry
 */
zend_module_entry rindow_openblas_module_entry = {
    STANDARD_MODULE_HEADER,
    "rindow_openblas",					/* Extension name */
    NULL,			                    /* zend_function_entry */
    PHP_MINIT(rindow_openblas),			/* PHP_MINIT - Module initialization */
    NULL,							    /* PHP_MSHUTDOWN - Module shutdown */
    PHP_RINIT(rindow_openblas),			/* PHP_RINIT - Request initialization */
    NULL,							    /* PHP_RSHUTDOWN - Request shutdown */
    PHP_MINFO(rindow_openblas),			/* PHP_MINFO - Module info */
    PHP_RINDOW_OPENBLAS_VERSION,		/* Version */
    STANDARD_MODULE_PROPERTIES
};
/* }}} */

#ifdef COMPILE_DL_RINDOW_OPENBLAS
# ifdef ZTS
ZEND_TSRMLS_CACHE_DEFINE()
# endif
ZEND_GET_MODULE(rindow_openblas)
#endif
