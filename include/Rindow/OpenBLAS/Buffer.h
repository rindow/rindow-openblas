#ifndef PHP_RINDOW_OPENBLAS_BUFFER_H
# define PHP_RINDOW_OPENBLAS_BUFFER_H

#define PHP_RINDOW_OPENBLAS_BUFFER_CLASSNAME "Rindow\\OpenBLAS\\Buffer"

typedef struct {
    zend_long size;
    zend_long dtype;
    zend_long value_size;
    void* data;
    zend_object std;
} php_rindow_openblas_buffer_t;
static inline php_rindow_openblas_buffer_t* php_rindow_openblas_buffer_fetch_object(zend_object* obj)
{
	return (php_rindow_openblas_buffer_t*) ((char*) obj - XtOffsetOf(php_rindow_openblas_buffer_t, std));
}
#define Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(zv) (php_rindow_openblas_buffer_fetch_object(Z_OBJ_P(zv)))

enum php_rindow_openblas_dtype {
    php_rindow_openblas_dtype_unknown = 0,
    php_rindow_openblas_dtype_bool    = 1,
    php_rindow_openblas_dtype_int8    = 2,
    php_rindow_openblas_dtype_int16   = 3,
    php_rindow_openblas_dtype_int32   = 4,
    php_rindow_openblas_dtype_int64   = 5,
    php_rindow_openblas_dtype_uint8   = 6,
    php_rindow_openblas_dtype_uint16  = 7,
    php_rindow_openblas_dtype_uint32  = 8,
    php_rindow_openblas_dtype_uint64  = 9,
    php_rindow_openblas_dtype_float8  = 10,
    php_rindow_openblas_dtype_float16 = 11,
    php_rindow_openblas_dtype_float32 = 12,
    php_rindow_openblas_dtype_float64 = 13,
};

#endif	/* PHP_RINDOW_OPENBLAS_BUFFER_H */
