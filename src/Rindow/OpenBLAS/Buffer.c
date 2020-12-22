#include <php.h>
#include <Zend/zend_interfaces.h>
#include <Zend/zend_exceptions.h>
#include <ext/spl/spl_iterators.h>
#include <ext/spl/spl_exceptions.h>
#include <stdint.h>
#include <Rindow/OpenBLAS/Buffer.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "php_rindow_openblas.h"


static zend_object_handlers rindow_openblas_buffer_object_handlers;

// destractor
static void php_rindow_openblas_buffer_free_object(zend_object* object)
{
    php_rindow_openblas_buffer_t* obj = php_rindow_openblas_buffer_fetch_object(object);
    if (obj->data) {
        efree(obj->data);
    }
    zend_object_std_dtor(&obj->std);
}

// constructor
static zend_object* php_rindow_openblas_buffer_create_object(zend_class_entry* class_type) /* {{{ */
{
    php_rindow_openblas_buffer_t* intern = NULL;

    intern = (php_rindow_openblas_buffer_t*)ecalloc(1, sizeof(php_rindow_openblas_buffer_t) + zend_object_properties_size(class_type));

    zend_object_std_init(&intern->std, class_type);
    object_properties_init(&intern->std, class_type);

    intern->std.handlers = &rindow_openblas_buffer_object_handlers;

    return &intern->std;
} /* }}} */

static void php_rindow_openblas_buffer_do_set(
    php_rindow_openblas_buffer_t* intern,
    zend_long offset,
    zend_long intvalue,
    double floatvalue)
{
    switch(intern->dtype) {
        case php_rindow_openblas_dtype_bool:
        case php_rindow_openblas_dtype_int8:
        case php_rindow_openblas_dtype_uint8:
            {
                uint8_t* bufint8=(uint8_t*)intern->data;
                bufint8[offset]=(uint8_t)(intvalue & 0xff);
            }
            break;
        case php_rindow_openblas_dtype_int16:
        case php_rindow_openblas_dtype_uint16:
            {
                int16_t* bufint16=(int16_t*)intern->data;
                bufint16[offset]=(int16_t)(intvalue & 0xffff);
            }
            break;
        case php_rindow_openblas_dtype_int32:
        case php_rindow_openblas_dtype_uint32:
            {
                int32_t* bufint32=(int32_t*)intern->data;
                bufint32[offset]=(int32_t)(intvalue);
            }
            break;
        case php_rindow_openblas_dtype_int64:
        case php_rindow_openblas_dtype_uint64:
            {
                int64_t* bufint64=(int64_t*)intern->data;
                bufint64[offset]=(int64_t)(intvalue);
            }
            break;
        case php_rindow_openblas_dtype_float32:
            {
                float* buffloat32=(float*)intern->data;
                buffloat32[offset]=(float)floatvalue;
            }
            break;
        case php_rindow_openblas_dtype_float64:
            {
                double* buffloat64=(double*)intern->data;
                buffloat64[offset]=(double)floatvalue;
            }
            break;
        default:
            zend_throw_exception(spl_ce_DomainException, "invalid dtype", 0);
    }
}

/* Method Rindow\OpenBLAS\Buffer::__construct($size,$dtype) {{{ */
static PHP_METHOD(Buffer, __construct)
{
    php_rindow_openblas_buffer_t* intern;
    zend_long size = 0;
    zend_long dtype = 0;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "ll",
    //        &size, &dtype) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "The size and the dtype are required", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 2, 2)
        Z_PARAM_LONG(size)
        Z_PARAM_LONG(dtype)
    ZEND_PARSE_PARAMETERS_END();

    switch(dtype) {
        case php_rindow_openblas_dtype_bool:
        case php_rindow_openblas_dtype_int8:
        case php_rindow_openblas_dtype_uint8:
        case php_rindow_openblas_dtype_int16:
        case php_rindow_openblas_dtype_uint16:
        case php_rindow_openblas_dtype_int32:
        case php_rindow_openblas_dtype_uint32:
        case php_rindow_openblas_dtype_int64:
        case php_rindow_openblas_dtype_uint64:
        case php_rindow_openblas_dtype_float32:
        case php_rindow_openblas_dtype_float64:
            break;
        default:
            zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
            return;
    }

    intern = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(getThis());
    if(size<=0 || dtype<=0) {
        intern->size = 0;
        intern->dtype = 0;
        intern->data = NULL;
        zend_throw_exception(spl_ce_InvalidArgumentException, "The size must be at least 1 and The dtype must not be 0.", 0);
        return;
    }
    intern->size = size;
    intern->dtype = dtype;
    intern->value_size = php_rindow_openblas_dtype_to_valuesize(dtype);
    if(intern->value_size==0) {
        intern->data = NULL;
        return;
    }
    intern->data = ecalloc(size, intern->value_size);
}
/* }}} */

/* Method Rindow\OpenBLAS\Buffer::offsetExists($offset) {{{ */
static PHP_METHOD(Buffer, offsetExists)
{
    php_rindow_openblas_buffer_t* intern;
    zend_long offset;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "l",
    //        &offset) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "1 parameter must be integer", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 1, 1)
        Z_PARAM_LONG(offset)
    ZEND_PARSE_PARAMETERS_END();

    intern = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(getThis());
    if(intern->data==NULL) {
        zend_throw_exception(spl_ce_DomainException, "uninitialized array", 0);
        return;
    }
    if(offset<0 || offset>=intern->size) {
        RETURN_FALSE;
    }
    RETURN_TRUE;
}
/* }}} */

/* Method Rindow\OpenBLAS\Buffer::offsetGet($offset) {{{ */
static PHP_METHOD(Buffer, offsetGet)
{
    php_rindow_openblas_buffer_t* intern;
    zend_long offset;
    zend_long intvalue;
    double floatvalue;
    zend_bool boolvalue;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "l",
    //        &offset) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "1 parameter must be integer", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 1, 1)
        Z_PARAM_LONG(offset)
    ZEND_PARSE_PARAMETERS_END();

    intern = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(getThis());
    if(intern->data==NULL) {
        zend_throw_exception(spl_ce_DomainException, "uninitialized array", 0);
        return;
    }
    if(offset<0 || offset>=intern->size) {
        zend_throw_exception(spl_ce_RuntimeException, "Index invalid or out of range", 0);
        return;
    }
    switch(intern->dtype) {
        case php_rindow_openblas_dtype_bool:
            {
                uint8_t* bufbool=(uint8_t*)intern->data;
                boolvalue =(zend_bool)bufbool[offset];
            }
            RETURN_BOOL(boolvalue);
        case php_rindow_openblas_dtype_int8:
            {
                int8_t* bufint8=(int8_t*)intern->data;
                intvalue =(zend_long)bufint8[offset];
            }
            RETURN_LONG(intvalue);
        case php_rindow_openblas_dtype_uint8:
            {
                uint8_t* bufuint8=(uint8_t*)intern->data;
                intvalue =(zend_long)bufuint8[offset];
            }
            RETURN_LONG(intvalue);
        case php_rindow_openblas_dtype_int16:
            {
                int16_t* bufint16=(int16_t*)intern->data;
                intvalue = (zend_long)bufint16[offset];
            }
            RETURN_LONG(intvalue);
        case php_rindow_openblas_dtype_uint16:
            {
                uint16_t* bufuint16=(uint16_t*)intern->data;
                intvalue = (zend_long)bufuint16[offset];
            }
            RETURN_LONG(intvalue);
        case php_rindow_openblas_dtype_int32:
            {
                int32_t* bufint32=(int32_t*)intern->data;
                intvalue = (zend_long)bufint32[offset];
            }
            RETURN_LONG(intvalue);
        case php_rindow_openblas_dtype_uint32:
            {
                uint32_t* bufuint32=(uint32_t*)intern->data;
                intvalue = (zend_long)bufuint32[offset];
            }
            RETURN_LONG(intvalue);
        case php_rindow_openblas_dtype_int64:
            {
                int64_t* bufint64=(int64_t*)intern->data;
                intvalue = (zend_long)bufint64[offset];
            }
            RETURN_LONG(intvalue);
        case php_rindow_openblas_dtype_uint64:
            {
                uint64_t* bufuint64=(uint64_t*)intern->data;
                intvalue = (zend_long)bufuint64[offset];
            }
            RETURN_LONG(intvalue);
        case php_rindow_openblas_dtype_float32:
            {
                float* buffloat32=(float*)intern->data;
                floatvalue = (double)buffloat32[offset];
            }
            RETURN_DOUBLE(floatvalue);
        case php_rindow_openblas_dtype_float64:
            {
                double* buffloat64=(double*)intern->data;
                floatvalue = (double)buffloat64[offset];
            }
            RETURN_DOUBLE(floatvalue);
        default:
            zend_throw_exception(spl_ce_DomainException, "invalid dtype", 0);
    }
}
/* }}} */

/* Method Rindow\OpenBLAS\Buffer::offsetSet($offset,$value) {{{ */
static PHP_METHOD(Buffer, offsetSet)
{
    php_rindow_openblas_buffer_t* intern;
    zend_long offset;
    zend_long intvalue;
    double floatvalue;
    zend_bool boolvalue;

    intern = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(getThis());
    if(intern->data==NULL) {
        zend_throw_exception(spl_ce_DomainException, "uninitialized array", 0);
        return;
    }
    if(php_rindow_openblas_dtype_is_int(intern->dtype)) {
        //if (zend_parse_parameters(ZEND_NUM_ARGS(), "ll",
        //        &offset,&intvalue) == FAILURE) {
        //    zend_throw_exception(spl_ce_InvalidArgumentException, "The offset and the value are required", 0);
        //    return;
        //}
        ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 2, 2)
            Z_PARAM_LONG(offset)
            Z_PARAM_LONG(intvalue)
        ZEND_PARSE_PARAMETERS_END();

    } else if(php_rindow_openblas_dtype_is_float(intern->dtype)) {
        //if (zend_parse_parameters(ZEND_NUM_ARGS(), "ld",
        //        &offset,&floatvalue) == FAILURE) {
        //    zend_throw_exception(spl_ce_InvalidArgumentException, "The offset and the value are required", 0);
        //    return;
        //}
        ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 2, 2)
            Z_PARAM_LONG(offset)
            Z_PARAM_DOUBLE(floatvalue)
        ZEND_PARSE_PARAMETERS_END();

    } else if(php_rindow_openblas_dtype_is_bool(intern->dtype)) {
        //if (zend_parse_parameters(ZEND_NUM_ARGS(), "lb",
        //        &offset,&boolvalue) == FAILURE) {
        //    zend_throw_exception(spl_ce_InvalidArgumentException, "The offset and the value are required", 0);
        //    return;
        //}
        ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 2, 2)
            Z_PARAM_LONG(offset)
            Z_PARAM_BOOL(boolvalue)
        ZEND_PARSE_PARAMETERS_END();

        if(boolvalue) {
            intvalue = 1;
        } else {
            intvalue = 0;
        }
    } else {
        zend_throw_exception(spl_ce_DomainException, "invalid dtype", 0);
        return;
    }

    if(offset<0 || offset>=intern->size) {
        zend_throw_exception(spl_ce_RuntimeException, "Index invalid or out of range", 0);
        return;
    }
    php_rindow_openblas_buffer_do_set(intern,offset,intvalue,floatvalue);
}

/* }}} */

/* Method Rindow\OpenBLAS\Buffer::offsetUnset($offset) {{{ */
static PHP_METHOD(Buffer, offsetUnset)
{
    php_rindow_openblas_buffer_t* intern;
    zend_long offset;

    intern = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(getThis());
    if(intern->data==NULL) {
        zend_throw_exception(spl_ce_DomainException, "uninitialized array", 0);
        return;
    }
    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "l",
    //        &offset) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "1 parameter must be integer", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 1, 1)
        Z_PARAM_LONG(offset)
    ZEND_PARSE_PARAMETERS_END();

    if(offset<0 || offset>=intern->size) {
        return;
    }
    php_rindow_openblas_buffer_do_set(intern,offset,0,0.0);
}
/* }}} */

/* Method Rindow\OpenBLAS\Buffer::count() {{{ */
static PHP_METHOD(Buffer, count)
{
    php_rindow_openblas_buffer_t* intern;

    intern = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(getThis());
    if(intern->data==NULL) {
        zend_throw_exception(spl_ce_DomainException, "uninitialized array", 0);
        return;
    }
    RETURN_LONG(intern->size);
}
/* }}} */

/* Method Rindow\OpenBLAS\Buffer::dtype() {{{ */
static PHP_METHOD(Buffer, dtype)
{
    php_rindow_openblas_buffer_t* intern;

    intern = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(getThis());
    if(intern->data==NULL) {
        zend_throw_exception(spl_ce_DomainException, "uninitialized array", 0);
        return;
    }
    RETURN_LONG(intern->dtype);
}
/* }}} */

/* Method Rindow\OpenBLAS\Buffer::value_size() {{{ */
static PHP_METHOD(Buffer, value_size)
{
    php_rindow_openblas_buffer_t* intern;

    intern = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(getThis());
    if(intern->data==NULL) {
        zend_throw_exception(spl_ce_DomainException, "uninitialized array", 0);
        return;
    }
    RETURN_LONG(intern->value_size);
}
/* }}} */

/* Method Rindow\OpenBLAS\Buffer::dump() : string {{{ */
static PHP_METHOD(Buffer, dump)
{
    php_rindow_openblas_buffer_t* intern;
    size_t len=0;

    intern = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(getThis());
    if(intern->data==NULL) {
        zend_throw_exception(spl_ce_DomainException, "uninitialized array", 0);
        return;
    }

    len = intern->size * php_rindow_openblas_dtype_to_valuesize(intern->dtype);
    RETURN_STRINGL(intern->data,len);
}
/* }}} */

/* Method Rindow\OpenBLAS\Buffer::load(string $data) {{{ */
static PHP_METHOD(Buffer, load)
{
    php_rindow_openblas_buffer_t* intern;
    size_t len=0;
    char* str=NULL;

    intern = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(getThis());
    if(intern->data==NULL) {
        zend_throw_exception(spl_ce_DomainException, "uninitialized array", 0);
        return;
    }
    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "s",
    //        &str,&len) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "1 parameter must be string", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 1, 1)
        Z_PARAM_STRING(str,len)
    ZEND_PARSE_PARAMETERS_END();

    if(len != intern->size * php_rindow_openblas_dtype_to_valuesize(intern->dtype)) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "unmatch data size", 0);
        return;
    }
    memcpy(intern->data,str,len);
}
/* }}} */

ZEND_BEGIN_ARG_INFO_EX(ai_Buffer___construct, 0, 0, 2)
    ZEND_ARG_INFO(0, size)
    ZEND_ARG_INFO(0, dtype)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Buffer_offsetExists, 0, 0, 1)
    ZEND_ARG_INFO(0, offset)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Buffer_offsetGet, 0, 0, 1)
    ZEND_ARG_INFO(0, offset)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Buffer_offsetSet, 0, 0, 2)
    ZEND_ARG_INFO(0, offset)
    ZEND_ARG_INFO(0, value)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Buffer_offsetUnset, 0, 0, 1)
    ZEND_ARG_INFO(0, offset)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Buffer_load, 0, 0, 1)
    ZEND_ARG_INFO(0, data)
ZEND_END_ARG_INFO()

ZEND_BEGIN_ARG_INFO_EX(ai_Buffer_void, 0, 0, 0)
ZEND_END_ARG_INFO()

/* {{{ Rindow\OpenBLAS\Buffer function entries */
static zend_function_entry php_rindow_openblas_buffer_me[] = {
    /* clang-format off */
    PHP_ME(Buffer, __construct, ai_Buffer___construct, ZEND_ACC_PUBLIC)
    PHP_ME(Buffer, offsetExists, ai_Buffer_offsetExists, ZEND_ACC_PUBLIC)
    PHP_ME(Buffer, offsetGet, ai_Buffer_offsetGet, ZEND_ACC_PUBLIC)
    PHP_ME(Buffer, offsetSet, ai_Buffer_offsetSet, ZEND_ACC_PUBLIC)
    PHP_ME(Buffer, offsetUnset, ai_Buffer_offsetUnset, ZEND_ACC_PUBLIC)
    PHP_ME(Buffer, count, ai_Buffer_void, ZEND_ACC_PUBLIC)
    PHP_ME(Buffer, value_size, ai_Buffer_void, ZEND_ACC_PUBLIC)
    PHP_ME(Buffer, dtype, ai_Buffer_void, ZEND_ACC_PUBLIC)
    PHP_ME(Buffer, dump, ai_Buffer_void, ZEND_ACC_PUBLIC)
    PHP_ME(Buffer, load, ai_Buffer_void, ZEND_ACC_PUBLIC)
    PHP_FE_END
    /* clang-format on */
};
/* }}} */

/* Class Rindow\OpenBLAS\Buffer {{{ */
zend_class_entry* php_rindow_openblas_buffer_ce;

void php_rindow_openblas_buffer_init_ce(INIT_FUNC_ARGS)
{
    zend_class_entry ce;

    INIT_NS_CLASS_ENTRY(ce, "Rindow\\OpenBLAS", "Buffer", php_rindow_openblas_buffer_me);
    php_rindow_openblas_buffer_ce = zend_register_internal_class(&ce);
    php_rindow_openblas_buffer_ce->create_object = php_rindow_openblas_buffer_create_object;

    memcpy(&rindow_openblas_buffer_object_handlers, zend_get_std_object_handlers(), sizeof(zend_object_handlers));
    rindow_openblas_buffer_object_handlers.offset    = XtOffsetOf(php_rindow_openblas_buffer_t, std);
    rindow_openblas_buffer_object_handlers.free_obj  = php_rindow_openblas_buffer_free_object;
    rindow_openblas_buffer_object_handlers.clone_obj = NULL;

    zend_class_implements(php_rindow_openblas_buffer_ce, 2, spl_ce_ArrayAccess, spl_ce_Countable);
}
/* }}} */
