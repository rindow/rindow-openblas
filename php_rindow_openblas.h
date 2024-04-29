/* rindow_openblas extension for PHP */

#ifndef PHP_RINDOW_OPENBLAS_H
# define PHP_RINDOW_OPENBLAS_H

# define phpext_rindow_openblas_ptr &rindow_openblas_module_entry

# define PHP_RINDOW_OPENBLAS_VERSION "0.4.1"

# if defined(ZTS) && defined(COMPILE_DL_RINDOW_OPENBLAS)
ZEND_TSRMLS_CACHE_EXTERN()
# endif


// Rindow\OpenBLAS\Blas object structures
typedef struct {
    zend_object std;
} php_rindow_openblas_blas_t;
static inline php_rindow_openblas_blas_t* php_rindow_openblas_blas_fetch_object(zend_object* obj)
{
	return (php_rindow_openblas_blas_t*) ((char*) obj - XtOffsetOf(php_rindow_openblas_blas_t, std));
}
#define Z_RINDOW_OPENBLAS_BLAS_OBJ_P(zv) (php_rindow_openblas_blas_fetch_object(Z_OBJ_P(zv)))

// Rindow\OpenBLAS\Lapack object structures
typedef struct {
    zend_object std;
} php_rindow_openblas_lapack_t;
static inline php_rindow_openblas_lapack_t* php_rindow_openblas_lapack_fetch_object(zend_object* obj)
{
	return (php_rindow_openblas_lapack_t*) ((char*) obj - XtOffsetOf(php_rindow_openblas_lapack_t, std));
}
#define Z_RINDOW_OPENBLAS_LAPACK_OBJ_P(zv) (php_rindow_openblas_lapack_fetch_object(Z_OBJ_P(zv)))

// Rindow\OpenBLAS\Math object structures
typedef struct {
    zend_object std;
} php_rindow_openblas_math_t;
static inline php_rindow_openblas_math_t* php_rindow_openblas_math_fetch_object(zend_object* obj)
{
	return (php_rindow_openblas_math_t*) ((char*) obj - XtOffsetOf(php_rindow_openblas_math_t, std));
}
#define Z_RINDOW_OPENBLAS_MATH_OBJ_P(zv) (php_rindow_openblas_math_fetch_object(Z_OBJ_P(zv)))

static inline void *php_rindow_openblas_get_address(
    php_interop_polite_math_matrix_linear_buffer_t* buffer,zend_long offset,int valueSize)
{
    return (uint8_t *)(buffer->data)+(offset*valueSize);
}


extern int php_rindow_openblas_common_dtype_to_valuesize(zend_long dtype);
extern int php_rindow_openblas_common_dtype_is_int(zend_long dtype);
extern int php_rindow_openblas_common_dtype_is_float(zend_long dtype);
extern int php_rindow_openblas_common_dtype_is_bool(zend_long dtype);

extern void php_rindow_openblas_buffer_init_ce(INIT_FUNC_ARGS);
extern void php_rindow_openblas_blas_init_ce(INIT_FUNC_ARGS);
extern void php_rindow_openblas_lapack_init_ce(INIT_FUNC_ARGS);
extern void php_rindow_openblas_math_init_ce(INIT_FUNC_ARGS);
extern zend_class_entry* php_rindow_openblas_buffer_ce;
extern zend_module_entry rindow_openblas_module_entry;

extern int php_rindow_openblas_assert_shape_parameter(
    char* name, zend_long n);
extern int php_rindow_openblas_assert_vector_buffer_spec(
    char* name,php_interop_polite_math_matrix_linear_buffer_t *buffer,
    zend_long n, zend_long offset, zend_long inc);
extern int php_rindow_openblas_assert_matrix_buffer_spec(
    char *name,php_interop_polite_math_matrix_linear_buffer_t *buffer,
    zend_long m, zend_long n, zend_long offset, zend_long ld);
extern int php_rindow_openblas_assert_buffer_size(
    php_interop_polite_math_matrix_linear_buffer_t *buffer,
    zend_long offset,zend_long size,
    char* message);
extern int php_rindow_openblas_assert_buffer_type(
    php_interop_polite_math_matrix_linear_buffer_t *buffer,
    char* name);

#endif	/* PHP_RINDOW_OPENBLAS_H */
