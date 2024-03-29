#include <stdio.h>
#include <Windows.h>
#include <cblas.h>
#if _MSC_VER
#include <complex.h>
#define lapack_complex_float _Fcomplex
#define lapack_complex_double _Dcomplex
#endif
#include <lapacke.h>

#define LOADFUNC(funcname) \
_g_##funcname = (PFN##funcname)GetProcAddress( _h_openblas, #funcname ); \
if(_g_##funcname==NULL) { \
    printf("load error: %s\n",  #funcname); \
    return -1; \
} \

static HMODULE _h_openblas = NULL;
typedef void (CALLBACK* PFNopenblas_set_num_threads)( /* openblas_set_num_threads */
    int            /* num_threads */
);
static PFNopenblas_set_num_threads _g_openblas_set_num_threads = NULL;
void openblas_set_num_threads(
    int            num_threads
)
{
    if(_h_openblas==NULL || _g_openblas_set_num_threads==NULL) {
        return;
    }
    _g_openblas_set_num_threads(
        num_threads    
    );
}
typedef void (CALLBACK* PFNgoto_set_num_threads)( /* goto_set_num_threads */
    int            /* num_threads */
);
static PFNgoto_set_num_threads _g_goto_set_num_threads = NULL;
void goto_set_num_threads(
    int            num_threads
)
{
    if(_h_openblas==NULL || _g_goto_set_num_threads==NULL) {
        return;
    }
    _g_goto_set_num_threads(
        num_threads    
    );
}
typedef int (CALLBACK* PFNopenblas_get_num_threads)( /* openblas_get_num_threads */
    void            /*  */
);
static PFNopenblas_get_num_threads _g_openblas_get_num_threads = NULL;
int openblas_get_num_threads(
    void            
)
{
    if(_h_openblas==NULL || _g_openblas_get_num_threads==NULL) {
        return 0;
    }
    return _g_openblas_get_num_threads(
    
    );
}
typedef int (CALLBACK* PFNopenblas_get_num_procs)( /* openblas_get_num_procs */
    void            /*  */
);
static PFNopenblas_get_num_procs _g_openblas_get_num_procs = NULL;
int openblas_get_num_procs(
    void            
)
{
    if(_h_openblas==NULL || _g_openblas_get_num_procs==NULL) {
        return 0;
    }
    return _g_openblas_get_num_procs(
    
    );
}
typedef char* (CALLBACK* PFNopenblas_get_config)( /* openblas_get_config */
    void            /*  */
);
static PFNopenblas_get_config _g_openblas_get_config = NULL;
char* openblas_get_config(
    void            
)
{
    if(_h_openblas==NULL || _g_openblas_get_config==NULL) {
        return 0;
    }
    return _g_openblas_get_config(
    
    );
}
typedef char* (CALLBACK* PFNopenblas_get_corename)( /* openblas_get_corename */
    void            /*  */
);
static PFNopenblas_get_corename _g_openblas_get_corename = NULL;
char* openblas_get_corename(
    void            
)
{
    if(_h_openblas==NULL || _g_openblas_get_corename==NULL) {
        return 0;
    }
    return _g_openblas_get_corename(
    
    );
}
typedef int (CALLBACK* PFNopenblas_get_parallel)( /* openblas_get_parallel */
    void            /*  */
);
static PFNopenblas_get_parallel _g_openblas_get_parallel = NULL;
int openblas_get_parallel(
    void            
)
{
    if(_h_openblas==NULL || _g_openblas_get_parallel==NULL) {
        return 0;
    }
    return _g_openblas_get_parallel(
    
    );
}
typedef float (CALLBACK* PFNcblas_sdsdot)( /* cblas_sdsdot */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST float            /* alpha */,
    OPENBLAS_CONST float *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    OPENBLAS_CONST float *            /* y */,
    OPENBLAS_CONST blasint            /* incy */
);
static PFNcblas_sdsdot _g_cblas_sdsdot = NULL;
float cblas_sdsdot(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST float            alpha,
    OPENBLAS_CONST float *            x,
    OPENBLAS_CONST blasint            incx,
    OPENBLAS_CONST float *            y,
    OPENBLAS_CONST blasint            incy
)
{
    if(_h_openblas==NULL || _g_cblas_sdsdot==NULL) {
        return 0;
    }
    return _g_cblas_sdsdot(
        n,
        alpha,
        x,
        incx,
        y,
        incy    
    );
}
typedef double (CALLBACK* PFNcblas_dsdot)( /* cblas_dsdot */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST float *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    OPENBLAS_CONST float *            /* y */,
    OPENBLAS_CONST blasint            /* incy */
);
static PFNcblas_dsdot _g_cblas_dsdot = NULL;
double cblas_dsdot(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST float *            x,
    OPENBLAS_CONST blasint            incx,
    OPENBLAS_CONST float *            y,
    OPENBLAS_CONST blasint            incy
)
{
    if(_h_openblas==NULL || _g_cblas_dsdot==NULL) {
        return 0;
    }
    return _g_cblas_dsdot(
        n,
        x,
        incx,
        y,
        incy    
    );
}
typedef float (CALLBACK* PFNcblas_sdot)( /* cblas_sdot */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST float *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    OPENBLAS_CONST float *            /* y */,
    OPENBLAS_CONST blasint            /* incy */
);
static PFNcblas_sdot _g_cblas_sdot = NULL;
float cblas_sdot(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST float *            x,
    OPENBLAS_CONST blasint            incx,
    OPENBLAS_CONST float *            y,
    OPENBLAS_CONST blasint            incy
)
{
    if(_h_openblas==NULL || _g_cblas_sdot==NULL) {
        return 0;
    }
    return _g_cblas_sdot(
        n,
        x,
        incx,
        y,
        incy    
    );
}
typedef double (CALLBACK* PFNcblas_ddot)( /* cblas_ddot */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST double *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    OPENBLAS_CONST double *            /* y */,
    OPENBLAS_CONST blasint            /* incy */
);
static PFNcblas_ddot _g_cblas_ddot = NULL;
double cblas_ddot(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST double *            x,
    OPENBLAS_CONST blasint            incx,
    OPENBLAS_CONST double *            y,
    OPENBLAS_CONST blasint            incy
)
{
    if(_h_openblas==NULL || _g_cblas_ddot==NULL) {
        return 0;
    }
    return _g_cblas_ddot(
        n,
        x,
        incx,
        y,
        incy    
    );
}
typedef void (CALLBACK* PFNcblas_cdotu_sub)( /* cblas_cdotu_sub */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    OPENBLAS_CONST void *            /* y */,
    OPENBLAS_CONST blasint            /* incy */,
    void *            /* ret */
);
static PFNcblas_cdotu_sub _g_cblas_cdotu_sub = NULL;
void cblas_cdotu_sub(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST void *            x,
    OPENBLAS_CONST blasint            incx,
    OPENBLAS_CONST void *            y,
    OPENBLAS_CONST blasint            incy,
    void *            ret
)
{
    if(_h_openblas==NULL || _g_cblas_cdotu_sub==NULL) {
        return;
    }
    _g_cblas_cdotu_sub(
        n,
        x,
        incx,
        y,
        incy,
        ret    
    );
}
typedef void (CALLBACK* PFNcblas_cdotc_sub)( /* cblas_cdotc_sub */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    OPENBLAS_CONST void *            /* y */,
    OPENBLAS_CONST blasint            /* incy */,
    void *            /* ret */
);
static PFNcblas_cdotc_sub _g_cblas_cdotc_sub = NULL;
void cblas_cdotc_sub(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST void *            x,
    OPENBLAS_CONST blasint            incx,
    OPENBLAS_CONST void *            y,
    OPENBLAS_CONST blasint            incy,
    void *            ret
)
{
    if(_h_openblas==NULL || _g_cblas_cdotc_sub==NULL) {
        return;
    }
    _g_cblas_cdotc_sub(
        n,
        x,
        incx,
        y,
        incy,
        ret    
    );
}
typedef void (CALLBACK* PFNcblas_zdotu_sub)( /* cblas_zdotu_sub */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    OPENBLAS_CONST void *            /* y */,
    OPENBLAS_CONST blasint            /* incy */,
    void *            /* ret */
);
static PFNcblas_zdotu_sub _g_cblas_zdotu_sub = NULL;
void cblas_zdotu_sub(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST void *            x,
    OPENBLAS_CONST blasint            incx,
    OPENBLAS_CONST void *            y,
    OPENBLAS_CONST blasint            incy,
    void *            ret
)
{
    if(_h_openblas==NULL || _g_cblas_zdotu_sub==NULL) {
        return;
    }
    _g_cblas_zdotu_sub(
        n,
        x,
        incx,
        y,
        incy,
        ret    
    );
}
typedef void (CALLBACK* PFNcblas_zdotc_sub)( /* cblas_zdotc_sub */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    OPENBLAS_CONST void *            /* y */,
    OPENBLAS_CONST blasint            /* incy */,
    void *            /* ret */
);
static PFNcblas_zdotc_sub _g_cblas_zdotc_sub = NULL;
void cblas_zdotc_sub(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST void *            x,
    OPENBLAS_CONST blasint            incx,
    OPENBLAS_CONST void *            y,
    OPENBLAS_CONST blasint            incy,
    void *            ret
)
{
    if(_h_openblas==NULL || _g_cblas_zdotc_sub==NULL) {
        return;
    }
    _g_cblas_zdotc_sub(
        n,
        x,
        incx,
        y,
        incy,
        ret    
    );
}
typedef float (CALLBACK* PFNcblas_sasum)( /* cblas_sasum */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST float *            /* x */,
    OPENBLAS_CONST blasint            /* incx */
);
static PFNcblas_sasum _g_cblas_sasum = NULL;
float cblas_sasum(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST float *            x,
    OPENBLAS_CONST blasint            incx
)
{
    if(_h_openblas==NULL || _g_cblas_sasum==NULL) {
        return 0;
    }
    return _g_cblas_sasum(
        n,
        x,
        incx    
    );
}
typedef double (CALLBACK* PFNcblas_dasum)( /* cblas_dasum */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST double *            /* x */,
    OPENBLAS_CONST blasint            /* incx */
);
static PFNcblas_dasum _g_cblas_dasum = NULL;
double cblas_dasum(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST double *            x,
    OPENBLAS_CONST blasint            incx
)
{
    if(_h_openblas==NULL || _g_cblas_dasum==NULL) {
        return 0;
    }
    return _g_cblas_dasum(
        n,
        x,
        incx    
    );
}
typedef float (CALLBACK* PFNcblas_scasum)( /* cblas_scasum */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */
);
static PFNcblas_scasum _g_cblas_scasum = NULL;
float cblas_scasum(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST void *            x,
    OPENBLAS_CONST blasint            incx
)
{
    if(_h_openblas==NULL || _g_cblas_scasum==NULL) {
        return 0;
    }
    return _g_cblas_scasum(
        n,
        x,
        incx    
    );
}
typedef double (CALLBACK* PFNcblas_dzasum)( /* cblas_dzasum */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */
);
static PFNcblas_dzasum _g_cblas_dzasum = NULL;
double cblas_dzasum(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST void *            x,
    OPENBLAS_CONST blasint            incx
)
{
    if(_h_openblas==NULL || _g_cblas_dzasum==NULL) {
        return 0;
    }
    return _g_cblas_dzasum(
        n,
        x,
        incx    
    );
}
typedef float (CALLBACK* PFNcblas_ssum)( /* cblas_ssum */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST float *            /* x */,
    OPENBLAS_CONST blasint            /* incx */
);
static PFNcblas_ssum _g_cblas_ssum = NULL;
float cblas_ssum(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST float *            x,
    OPENBLAS_CONST blasint            incx
)
{
    if(_h_openblas==NULL || _g_cblas_ssum==NULL) {
        return 0;
    }
    return _g_cblas_ssum(
        n,
        x,
        incx    
    );
}
typedef double (CALLBACK* PFNcblas_dsum)( /* cblas_dsum */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST double *            /* x */,
    OPENBLAS_CONST blasint            /* incx */
);
static PFNcblas_dsum _g_cblas_dsum = NULL;
double cblas_dsum(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST double *            x,
    OPENBLAS_CONST blasint            incx
)
{
    if(_h_openblas==NULL || _g_cblas_dsum==NULL) {
        return 0;
    }
    return _g_cblas_dsum(
        n,
        x,
        incx    
    );
}
typedef float (CALLBACK* PFNcblas_scsum)( /* cblas_scsum */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */
);
static PFNcblas_scsum _g_cblas_scsum = NULL;
float cblas_scsum(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST void *            x,
    OPENBLAS_CONST blasint            incx
)
{
    if(_h_openblas==NULL || _g_cblas_scsum==NULL) {
        return 0;
    }
    return _g_cblas_scsum(
        n,
        x,
        incx    
    );
}
typedef double (CALLBACK* PFNcblas_dzsum)( /* cblas_dzsum */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */
);
static PFNcblas_dzsum _g_cblas_dzsum = NULL;
double cblas_dzsum(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST void *            x,
    OPENBLAS_CONST blasint            incx
)
{
    if(_h_openblas==NULL || _g_cblas_dzsum==NULL) {
        return 0;
    }
    return _g_cblas_dzsum(
        n,
        x,
        incx    
    );
}
typedef float (CALLBACK* PFNcblas_snrm2)( /* cblas_snrm2 */
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST float *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_snrm2 _g_cblas_snrm2 = NULL;
float cblas_snrm2(
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST float *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_snrm2==NULL) {
        return 0;
    }
    return _g_cblas_snrm2(
        N,
        X,
        incX    
    );
}
typedef double (CALLBACK* PFNcblas_dnrm2)( /* cblas_dnrm2 */
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST double *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_dnrm2 _g_cblas_dnrm2 = NULL;
double cblas_dnrm2(
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST double *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_dnrm2==NULL) {
        return 0;
    }
    return _g_cblas_dnrm2(
        N,
        X,
        incX    
    );
}
typedef float (CALLBACK* PFNcblas_scnrm2)( /* cblas_scnrm2 */
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_scnrm2 _g_cblas_scnrm2 = NULL;
float cblas_scnrm2(
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_scnrm2==NULL) {
        return 0;
    }
    return _g_cblas_scnrm2(
        N,
        X,
        incX    
    );
}
typedef double (CALLBACK* PFNcblas_dznrm2)( /* cblas_dznrm2 */
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_dznrm2 _g_cblas_dznrm2 = NULL;
double cblas_dznrm2(
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_dznrm2==NULL) {
        return 0;
    }
    return _g_cblas_dznrm2(
        N,
        X,
        incX    
    );
}
typedef CBLAS_INDEX (CALLBACK* PFNcblas_isamax)( /* cblas_isamax */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST float *            /* x */,
    OPENBLAS_CONST blasint            /* incx */
);
static PFNcblas_isamax _g_cblas_isamax = NULL;
CBLAS_INDEX cblas_isamax(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST float *            x,
    OPENBLAS_CONST blasint            incx
)
{
    if(_h_openblas==NULL || _g_cblas_isamax==NULL) {
        return 0;
    }
    return _g_cblas_isamax(
        n,
        x,
        incx    
    );
}
typedef CBLAS_INDEX (CALLBACK* PFNcblas_idamax)( /* cblas_idamax */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST double *            /* x */,
    OPENBLAS_CONST blasint            /* incx */
);
static PFNcblas_idamax _g_cblas_idamax = NULL;
CBLAS_INDEX cblas_idamax(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST double *            x,
    OPENBLAS_CONST blasint            incx
)
{
    if(_h_openblas==NULL || _g_cblas_idamax==NULL) {
        return 0;
    }
    return _g_cblas_idamax(
        n,
        x,
        incx    
    );
}
typedef CBLAS_INDEX (CALLBACK* PFNcblas_icamax)( /* cblas_icamax */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */
);
static PFNcblas_icamax _g_cblas_icamax = NULL;
CBLAS_INDEX cblas_icamax(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST void *            x,
    OPENBLAS_CONST blasint            incx
)
{
    if(_h_openblas==NULL || _g_cblas_icamax==NULL) {
        return 0;
    }
    return _g_cblas_icamax(
        n,
        x,
        incx    
    );
}
typedef CBLAS_INDEX (CALLBACK* PFNcblas_izamax)( /* cblas_izamax */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */
);
static PFNcblas_izamax _g_cblas_izamax = NULL;
CBLAS_INDEX cblas_izamax(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST void *            x,
    OPENBLAS_CONST blasint            incx
)
{
    if(_h_openblas==NULL || _g_cblas_izamax==NULL) {
        return 0;
    }
    return _g_cblas_izamax(
        n,
        x,
        incx    
    );
}
typedef CBLAS_INDEX (CALLBACK* PFNcblas_isamin)( /* cblas_isamin */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST float *            /* x */,
    OPENBLAS_CONST blasint            /* incx */
);
static PFNcblas_isamin _g_cblas_isamin = NULL;
CBLAS_INDEX cblas_isamin(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST float *            x,
    OPENBLAS_CONST blasint            incx
)
{
    if(_h_openblas==NULL || _g_cblas_isamin==NULL) {
        return 0;
    }
    return _g_cblas_isamin(
        n,
        x,
        incx    
    );
}
typedef CBLAS_INDEX (CALLBACK* PFNcblas_idamin)( /* cblas_idamin */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST double *            /* x */,
    OPENBLAS_CONST blasint            /* incx */
);
static PFNcblas_idamin _g_cblas_idamin = NULL;
CBLAS_INDEX cblas_idamin(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST double *            x,
    OPENBLAS_CONST blasint            incx
)
{
    if(_h_openblas==NULL || _g_cblas_idamin==NULL) {
        return 0;
    }
    return _g_cblas_idamin(
        n,
        x,
        incx    
    );
}
typedef CBLAS_INDEX (CALLBACK* PFNcblas_icamin)( /* cblas_icamin */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */
);
static PFNcblas_icamin _g_cblas_icamin = NULL;
CBLAS_INDEX cblas_icamin(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST void *            x,
    OPENBLAS_CONST blasint            incx
)
{
    if(_h_openblas==NULL || _g_cblas_icamin==NULL) {
        return 0;
    }
    return _g_cblas_icamin(
        n,
        x,
        incx    
    );
}
typedef CBLAS_INDEX (CALLBACK* PFNcblas_izamin)( /* cblas_izamin */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */
);
static PFNcblas_izamin _g_cblas_izamin = NULL;
CBLAS_INDEX cblas_izamin(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST void *            x,
    OPENBLAS_CONST blasint            incx
)
{
    if(_h_openblas==NULL || _g_cblas_izamin==NULL) {
        return 0;
    }
    return _g_cblas_izamin(
        n,
        x,
        incx    
    );
}
typedef CBLAS_INDEX (CALLBACK* PFNcblas_ismax)( /* cblas_ismax */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST float *            /* x */,
    OPENBLAS_CONST blasint            /* incx */
);
static PFNcblas_ismax _g_cblas_ismax = NULL;
CBLAS_INDEX cblas_ismax(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST float *            x,
    OPENBLAS_CONST blasint            incx
)
{
    if(_h_openblas==NULL || _g_cblas_ismax==NULL) {
        return 0;
    }
    return _g_cblas_ismax(
        n,
        x,
        incx    
    );
}
typedef CBLAS_INDEX (CALLBACK* PFNcblas_idmax)( /* cblas_idmax */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST double *            /* x */,
    OPENBLAS_CONST blasint            /* incx */
);
static PFNcblas_idmax _g_cblas_idmax = NULL;
CBLAS_INDEX cblas_idmax(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST double *            x,
    OPENBLAS_CONST blasint            incx
)
{
    if(_h_openblas==NULL || _g_cblas_idmax==NULL) {
        return 0;
    }
    return _g_cblas_idmax(
        n,
        x,
        incx    
    );
}
typedef CBLAS_INDEX (CALLBACK* PFNcblas_icmax)( /* cblas_icmax */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */
);
static PFNcblas_icmax _g_cblas_icmax = NULL;
CBLAS_INDEX cblas_icmax(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST void *            x,
    OPENBLAS_CONST blasint            incx
)
{
    if(_h_openblas==NULL || _g_cblas_icmax==NULL) {
        return 0;
    }
    return _g_cblas_icmax(
        n,
        x,
        incx    
    );
}
typedef CBLAS_INDEX (CALLBACK* PFNcblas_izmax)( /* cblas_izmax */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */
);
static PFNcblas_izmax _g_cblas_izmax = NULL;
CBLAS_INDEX cblas_izmax(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST void *            x,
    OPENBLAS_CONST blasint            incx
)
{
    if(_h_openblas==NULL || _g_cblas_izmax==NULL) {
        return 0;
    }
    return _g_cblas_izmax(
        n,
        x,
        incx    
    );
}
typedef CBLAS_INDEX (CALLBACK* PFNcblas_ismin)( /* cblas_ismin */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST float *            /* x */,
    OPENBLAS_CONST blasint            /* incx */
);
static PFNcblas_ismin _g_cblas_ismin = NULL;
CBLAS_INDEX cblas_ismin(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST float *            x,
    OPENBLAS_CONST blasint            incx
)
{
    if(_h_openblas==NULL || _g_cblas_ismin==NULL) {
        return 0;
    }
    return _g_cblas_ismin(
        n,
        x,
        incx    
    );
}
typedef CBLAS_INDEX (CALLBACK* PFNcblas_idmin)( /* cblas_idmin */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST double *            /* x */,
    OPENBLAS_CONST blasint            /* incx */
);
static PFNcblas_idmin _g_cblas_idmin = NULL;
CBLAS_INDEX cblas_idmin(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST double *            x,
    OPENBLAS_CONST blasint            incx
)
{
    if(_h_openblas==NULL || _g_cblas_idmin==NULL) {
        return 0;
    }
    return _g_cblas_idmin(
        n,
        x,
        incx    
    );
}
typedef CBLAS_INDEX (CALLBACK* PFNcblas_icmin)( /* cblas_icmin */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */
);
static PFNcblas_icmin _g_cblas_icmin = NULL;
CBLAS_INDEX cblas_icmin(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST void *            x,
    OPENBLAS_CONST blasint            incx
)
{
    if(_h_openblas==NULL || _g_cblas_icmin==NULL) {
        return 0;
    }
    return _g_cblas_icmin(
        n,
        x,
        incx    
    );
}
typedef CBLAS_INDEX (CALLBACK* PFNcblas_izmin)( /* cblas_izmin */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */
);
static PFNcblas_izmin _g_cblas_izmin = NULL;
CBLAS_INDEX cblas_izmin(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST void *            x,
    OPENBLAS_CONST blasint            incx
)
{
    if(_h_openblas==NULL || _g_cblas_izmin==NULL) {
        return 0;
    }
    return _g_cblas_izmin(
        n,
        x,
        incx    
    );
}
typedef void (CALLBACK* PFNcblas_saxpy)( /* cblas_saxpy */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST float            /* alpha */,
    OPENBLAS_CONST float *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    float *            /* y */,
    OPENBLAS_CONST blasint            /* incy */
);
static PFNcblas_saxpy _g_cblas_saxpy = NULL;
void cblas_saxpy(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST float            alpha,
    OPENBLAS_CONST float *            x,
    OPENBLAS_CONST blasint            incx,
    float *            y,
    OPENBLAS_CONST blasint            incy
)
{
    if(_h_openblas==NULL || _g_cblas_saxpy==NULL) {
        return;
    }
    _g_cblas_saxpy(
        n,
        alpha,
        x,
        incx,
        y,
        incy    
    );
}
typedef void (CALLBACK* PFNcblas_daxpy)( /* cblas_daxpy */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST double            /* alpha */,
    OPENBLAS_CONST double *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    double *            /* y */,
    OPENBLAS_CONST blasint            /* incy */
);
static PFNcblas_daxpy _g_cblas_daxpy = NULL;
void cblas_daxpy(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST double            alpha,
    OPENBLAS_CONST double *            x,
    OPENBLAS_CONST blasint            incx,
    double *            y,
    OPENBLAS_CONST blasint            incy
)
{
    if(_h_openblas==NULL || _g_cblas_daxpy==NULL) {
        return;
    }
    _g_cblas_daxpy(
        n,
        alpha,
        x,
        incx,
        y,
        incy    
    );
}
typedef void (CALLBACK* PFNcblas_caxpy)( /* cblas_caxpy */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    void *            /* y */,
    OPENBLAS_CONST blasint            /* incy */
);
static PFNcblas_caxpy _g_cblas_caxpy = NULL;
void cblas_caxpy(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            x,
    OPENBLAS_CONST blasint            incx,
    void *            y,
    OPENBLAS_CONST blasint            incy
)
{
    if(_h_openblas==NULL || _g_cblas_caxpy==NULL) {
        return;
    }
    _g_cblas_caxpy(
        n,
        alpha,
        x,
        incx,
        y,
        incy    
    );
}
typedef void (CALLBACK* PFNcblas_zaxpy)( /* cblas_zaxpy */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    void *            /* y */,
    OPENBLAS_CONST blasint            /* incy */
);
static PFNcblas_zaxpy _g_cblas_zaxpy = NULL;
void cblas_zaxpy(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            x,
    OPENBLAS_CONST blasint            incx,
    void *            y,
    OPENBLAS_CONST blasint            incy
)
{
    if(_h_openblas==NULL || _g_cblas_zaxpy==NULL) {
        return;
    }
    _g_cblas_zaxpy(
        n,
        alpha,
        x,
        incx,
        y,
        incy    
    );
}
typedef void (CALLBACK* PFNcblas_scopy)( /* cblas_scopy */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST float *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    float *            /* y */,
    OPENBLAS_CONST blasint            /* incy */
);
static PFNcblas_scopy _g_cblas_scopy = NULL;
void cblas_scopy(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST float *            x,
    OPENBLAS_CONST blasint            incx,
    float *            y,
    OPENBLAS_CONST blasint            incy
)
{
    if(_h_openblas==NULL || _g_cblas_scopy==NULL) {
        return;
    }
    _g_cblas_scopy(
        n,
        x,
        incx,
        y,
        incy    
    );
}
typedef void (CALLBACK* PFNcblas_dcopy)( /* cblas_dcopy */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST double *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    double *            /* y */,
    OPENBLAS_CONST blasint            /* incy */
);
static PFNcblas_dcopy _g_cblas_dcopy = NULL;
void cblas_dcopy(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST double *            x,
    OPENBLAS_CONST blasint            incx,
    double *            y,
    OPENBLAS_CONST blasint            incy
)
{
    if(_h_openblas==NULL || _g_cblas_dcopy==NULL) {
        return;
    }
    _g_cblas_dcopy(
        n,
        x,
        incx,
        y,
        incy    
    );
}
typedef void (CALLBACK* PFNcblas_ccopy)( /* cblas_ccopy */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    void *            /* y */,
    OPENBLAS_CONST blasint            /* incy */
);
static PFNcblas_ccopy _g_cblas_ccopy = NULL;
void cblas_ccopy(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST void *            x,
    OPENBLAS_CONST blasint            incx,
    void *            y,
    OPENBLAS_CONST blasint            incy
)
{
    if(_h_openblas==NULL || _g_cblas_ccopy==NULL) {
        return;
    }
    _g_cblas_ccopy(
        n,
        x,
        incx,
        y,
        incy    
    );
}
typedef void (CALLBACK* PFNcblas_zcopy)( /* cblas_zcopy */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    void *            /* y */,
    OPENBLAS_CONST blasint            /* incy */
);
static PFNcblas_zcopy _g_cblas_zcopy = NULL;
void cblas_zcopy(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST void *            x,
    OPENBLAS_CONST blasint            incx,
    void *            y,
    OPENBLAS_CONST blasint            incy
)
{
    if(_h_openblas==NULL || _g_cblas_zcopy==NULL) {
        return;
    }
    _g_cblas_zcopy(
        n,
        x,
        incx,
        y,
        incy    
    );
}
typedef void (CALLBACK* PFNcblas_sswap)( /* cblas_sswap */
    OPENBLAS_CONST blasint            /* n */,
    float *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    float *            /* y */,
    OPENBLAS_CONST blasint            /* incy */
);
static PFNcblas_sswap _g_cblas_sswap = NULL;
void cblas_sswap(
    OPENBLAS_CONST blasint            n,
    float *            x,
    OPENBLAS_CONST blasint            incx,
    float *            y,
    OPENBLAS_CONST blasint            incy
)
{
    if(_h_openblas==NULL || _g_cblas_sswap==NULL) {
        return;
    }
    _g_cblas_sswap(
        n,
        x,
        incx,
        y,
        incy    
    );
}
typedef void (CALLBACK* PFNcblas_dswap)( /* cblas_dswap */
    OPENBLAS_CONST blasint            /* n */,
    double *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    double *            /* y */,
    OPENBLAS_CONST blasint            /* incy */
);
static PFNcblas_dswap _g_cblas_dswap = NULL;
void cblas_dswap(
    OPENBLAS_CONST blasint            n,
    double *            x,
    OPENBLAS_CONST blasint            incx,
    double *            y,
    OPENBLAS_CONST blasint            incy
)
{
    if(_h_openblas==NULL || _g_cblas_dswap==NULL) {
        return;
    }
    _g_cblas_dswap(
        n,
        x,
        incx,
        y,
        incy    
    );
}
typedef void (CALLBACK* PFNcblas_cswap)( /* cblas_cswap */
    OPENBLAS_CONST blasint            /* n */,
    void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    void *            /* y */,
    OPENBLAS_CONST blasint            /* incy */
);
static PFNcblas_cswap _g_cblas_cswap = NULL;
void cblas_cswap(
    OPENBLAS_CONST blasint            n,
    void *            x,
    OPENBLAS_CONST blasint            incx,
    void *            y,
    OPENBLAS_CONST blasint            incy
)
{
    if(_h_openblas==NULL || _g_cblas_cswap==NULL) {
        return;
    }
    _g_cblas_cswap(
        n,
        x,
        incx,
        y,
        incy    
    );
}
typedef void (CALLBACK* PFNcblas_zswap)( /* cblas_zswap */
    OPENBLAS_CONST blasint            /* n */,
    void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    void *            /* y */,
    OPENBLAS_CONST blasint            /* incy */
);
static PFNcblas_zswap _g_cblas_zswap = NULL;
void cblas_zswap(
    OPENBLAS_CONST blasint            n,
    void *            x,
    OPENBLAS_CONST blasint            incx,
    void *            y,
    OPENBLAS_CONST blasint            incy
)
{
    if(_h_openblas==NULL || _g_cblas_zswap==NULL) {
        return;
    }
    _g_cblas_zswap(
        n,
        x,
        incx,
        y,
        incy    
    );
}
typedef void (CALLBACK* PFNcblas_srot)( /* cblas_srot */
    OPENBLAS_CONST blasint            /* N */,
    float *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    float *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */,
    OPENBLAS_CONST float            /* c */,
    OPENBLAS_CONST float            /* s */
);
static PFNcblas_srot _g_cblas_srot = NULL;
void cblas_srot(
    OPENBLAS_CONST blasint            N,
    float *            X,
    OPENBLAS_CONST blasint            incX,
    float *            Y,
    OPENBLAS_CONST blasint            incY,
    OPENBLAS_CONST float            c,
    OPENBLAS_CONST float            s
)
{
    if(_h_openblas==NULL || _g_cblas_srot==NULL) {
        return;
    }
    _g_cblas_srot(
        N,
        X,
        incX,
        Y,
        incY,
        c,
        s    
    );
}
typedef void (CALLBACK* PFNcblas_drot)( /* cblas_drot */
    OPENBLAS_CONST blasint            /* N */,
    double *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    double *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */,
    OPENBLAS_CONST double            /* c */,
    OPENBLAS_CONST double            /* s */
);
static PFNcblas_drot _g_cblas_drot = NULL;
void cblas_drot(
    OPENBLAS_CONST blasint            N,
    double *            X,
    OPENBLAS_CONST blasint            incX,
    double *            Y,
    OPENBLAS_CONST blasint            incY,
    OPENBLAS_CONST double            c,
    OPENBLAS_CONST double            s
)
{
    if(_h_openblas==NULL || _g_cblas_drot==NULL) {
        return;
    }
    _g_cblas_drot(
        N,
        X,
        incX,
        Y,
        incY,
        c,
        s    
    );
}
typedef void (CALLBACK* PFNcblas_srotg)( /* cblas_srotg */
    float *            /* a */,
    float *            /* b */,
    float *            /* c */,
    float *            /* s */
);
static PFNcblas_srotg _g_cblas_srotg = NULL;
void cblas_srotg(
    float *            a,
    float *            b,
    float *            c,
    float *            s
)
{
    if(_h_openblas==NULL || _g_cblas_srotg==NULL) {
        return;
    }
    _g_cblas_srotg(
        a,
        b,
        c,
        s    
    );
}
typedef void (CALLBACK* PFNcblas_drotg)( /* cblas_drotg */
    double *            /* a */,
    double *            /* b */,
    double *            /* c */,
    double *            /* s */
);
static PFNcblas_drotg _g_cblas_drotg = NULL;
void cblas_drotg(
    double *            a,
    double *            b,
    double *            c,
    double *            s
)
{
    if(_h_openblas==NULL || _g_cblas_drotg==NULL) {
        return;
    }
    _g_cblas_drotg(
        a,
        b,
        c,
        s    
    );
}
typedef void (CALLBACK* PFNcblas_srotm)( /* cblas_srotm */
    OPENBLAS_CONST blasint            /* N */,
    float *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    float *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */,
    OPENBLAS_CONST float *            /* P */
);
static PFNcblas_srotm _g_cblas_srotm = NULL;
void cblas_srotm(
    OPENBLAS_CONST blasint            N,
    float *            X,
    OPENBLAS_CONST blasint            incX,
    float *            Y,
    OPENBLAS_CONST blasint            incY,
    OPENBLAS_CONST float *            P
)
{
    if(_h_openblas==NULL || _g_cblas_srotm==NULL) {
        return;
    }
    _g_cblas_srotm(
        N,
        X,
        incX,
        Y,
        incY,
        P    
    );
}
typedef void (CALLBACK* PFNcblas_drotm)( /* cblas_drotm */
    OPENBLAS_CONST blasint            /* N */,
    double *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    double *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */,
    OPENBLAS_CONST double *            /* P */
);
static PFNcblas_drotm _g_cblas_drotm = NULL;
void cblas_drotm(
    OPENBLAS_CONST blasint            N,
    double *            X,
    OPENBLAS_CONST blasint            incX,
    double *            Y,
    OPENBLAS_CONST blasint            incY,
    OPENBLAS_CONST double *            P
)
{
    if(_h_openblas==NULL || _g_cblas_drotm==NULL) {
        return;
    }
    _g_cblas_drotm(
        N,
        X,
        incX,
        Y,
        incY,
        P    
    );
}
typedef void (CALLBACK* PFNcblas_srotmg)( /* cblas_srotmg */
    float *            /* d1 */,
    float *            /* d2 */,
    float *            /* b1 */,
    OPENBLAS_CONST float            /* b2 */,
    float *            /* P */
);
static PFNcblas_srotmg _g_cblas_srotmg = NULL;
void cblas_srotmg(
    float *            d1,
    float *            d2,
    float *            b1,
    OPENBLAS_CONST float            b2,
    float *            P
)
{
    if(_h_openblas==NULL || _g_cblas_srotmg==NULL) {
        return;
    }
    _g_cblas_srotmg(
        d1,
        d2,
        b1,
        b2,
        P    
    );
}
typedef void (CALLBACK* PFNcblas_drotmg)( /* cblas_drotmg */
    double *            /* d1 */,
    double *            /* d2 */,
    double *            /* b1 */,
    OPENBLAS_CONST double            /* b2 */,
    double *            /* P */
);
static PFNcblas_drotmg _g_cblas_drotmg = NULL;
void cblas_drotmg(
    double *            d1,
    double *            d2,
    double *            b1,
    OPENBLAS_CONST double            b2,
    double *            P
)
{
    if(_h_openblas==NULL || _g_cblas_drotmg==NULL) {
        return;
    }
    _g_cblas_drotmg(
        d1,
        d2,
        b1,
        b2,
        P    
    );
}
typedef void (CALLBACK* PFNcblas_sscal)( /* cblas_sscal */
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST float            /* alpha */,
    float *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_sscal _g_cblas_sscal = NULL;
void cblas_sscal(
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST float            alpha,
    float *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_sscal==NULL) {
        return;
    }
    _g_cblas_sscal(
        N,
        alpha,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_dscal)( /* cblas_dscal */
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST double            /* alpha */,
    double *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_dscal _g_cblas_dscal = NULL;
void cblas_dscal(
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST double            alpha,
    double *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_dscal==NULL) {
        return;
    }
    _g_cblas_dscal(
        N,
        alpha,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_cscal)( /* cblas_cscal */
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* alpha */,
    void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_cscal _g_cblas_cscal = NULL;
void cblas_cscal(
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            alpha,
    void *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_cscal==NULL) {
        return;
    }
    _g_cblas_cscal(
        N,
        alpha,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_zscal)( /* cblas_zscal */
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* alpha */,
    void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_zscal _g_cblas_zscal = NULL;
void cblas_zscal(
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            alpha,
    void *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_zscal==NULL) {
        return;
    }
    _g_cblas_zscal(
        N,
        alpha,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_csscal)( /* cblas_csscal */
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST float            /* alpha */,
    void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_csscal _g_cblas_csscal = NULL;
void cblas_csscal(
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST float            alpha,
    void *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_csscal==NULL) {
        return;
    }
    _g_cblas_csscal(
        N,
        alpha,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_zdscal)( /* cblas_zdscal */
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST double            /* alpha */,
    void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_zdscal _g_cblas_zdscal = NULL;
void cblas_zdscal(
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST double            alpha,
    void *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_zdscal==NULL) {
        return;
    }
    _g_cblas_zdscal(
        N,
        alpha,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_sgemv)( /* cblas_sgemv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* trans */,
    OPENBLAS_CONST blasint            /* m */,
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST float            /* alpha */,
    OPENBLAS_CONST float *            /* a */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST float *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    OPENBLAS_CONST float            /* beta */,
    float *            /* y */,
    OPENBLAS_CONST blasint            /* incy */
);
static PFNcblas_sgemv _g_cblas_sgemv = NULL;
void cblas_sgemv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            trans,
    OPENBLAS_CONST blasint            m,
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST float            alpha,
    OPENBLAS_CONST float *            a,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST float *            x,
    OPENBLAS_CONST blasint            incx,
    OPENBLAS_CONST float            beta,
    float *            y,
    OPENBLAS_CONST blasint            incy
)
{
    if(_h_openblas==NULL || _g_cblas_sgemv==NULL) {
        return;
    }
    _g_cblas_sgemv(
        order,
        trans,
        m,
        n,
        alpha,
        a,
        lda,
        x,
        incx,
        beta,
        y,
        incy    
    );
}
typedef void (CALLBACK* PFNcblas_dgemv)( /* cblas_dgemv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* trans */,
    OPENBLAS_CONST blasint            /* m */,
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST double            /* alpha */,
    OPENBLAS_CONST double *            /* a */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST double *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    OPENBLAS_CONST double            /* beta */,
    double *            /* y */,
    OPENBLAS_CONST blasint            /* incy */
);
static PFNcblas_dgemv _g_cblas_dgemv = NULL;
void cblas_dgemv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            trans,
    OPENBLAS_CONST blasint            m,
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST double            alpha,
    OPENBLAS_CONST double *            a,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST double *            x,
    OPENBLAS_CONST blasint            incx,
    OPENBLAS_CONST double            beta,
    double *            y,
    OPENBLAS_CONST blasint            incy
)
{
    if(_h_openblas==NULL || _g_cblas_dgemv==NULL) {
        return;
    }
    _g_cblas_dgemv(
        order,
        trans,
        m,
        n,
        alpha,
        a,
        lda,
        x,
        incx,
        beta,
        y,
        incy    
    );
}
typedef void (CALLBACK* PFNcblas_cgemv)( /* cblas_cgemv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* trans */,
    OPENBLAS_CONST blasint            /* m */,
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* a */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    OPENBLAS_CONST void *            /* beta */,
    void *            /* y */,
    OPENBLAS_CONST blasint            /* incy */
);
static PFNcblas_cgemv _g_cblas_cgemv = NULL;
void cblas_cgemv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            trans,
    OPENBLAS_CONST blasint            m,
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            a,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST void *            x,
    OPENBLAS_CONST blasint            incx,
    OPENBLAS_CONST void *            beta,
    void *            y,
    OPENBLAS_CONST blasint            incy
)
{
    if(_h_openblas==NULL || _g_cblas_cgemv==NULL) {
        return;
    }
    _g_cblas_cgemv(
        order,
        trans,
        m,
        n,
        alpha,
        a,
        lda,
        x,
        incx,
        beta,
        y,
        incy    
    );
}
typedef void (CALLBACK* PFNcblas_zgemv)( /* cblas_zgemv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* trans */,
    OPENBLAS_CONST blasint            /* m */,
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* a */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    OPENBLAS_CONST void *            /* beta */,
    void *            /* y */,
    OPENBLAS_CONST blasint            /* incy */
);
static PFNcblas_zgemv _g_cblas_zgemv = NULL;
void cblas_zgemv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            trans,
    OPENBLAS_CONST blasint            m,
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            a,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST void *            x,
    OPENBLAS_CONST blasint            incx,
    OPENBLAS_CONST void *            beta,
    void *            y,
    OPENBLAS_CONST blasint            incy
)
{
    if(_h_openblas==NULL || _g_cblas_zgemv==NULL) {
        return;
    }
    _g_cblas_zgemv(
        order,
        trans,
        m,
        n,
        alpha,
        a,
        lda,
        x,
        incx,
        beta,
        y,
        incy    
    );
}
typedef void (CALLBACK* PFNcblas_sger)( /* cblas_sger */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST float            /* alpha */,
    OPENBLAS_CONST float *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST float *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */,
    float *            /* A */,
    OPENBLAS_CONST blasint            /* lda */
);
static PFNcblas_sger _g_cblas_sger = NULL;
void cblas_sger(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST float            alpha,
    OPENBLAS_CONST float *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST float *            Y,
    OPENBLAS_CONST blasint            incY,
    float *            A,
    OPENBLAS_CONST blasint            lda
)
{
    if(_h_openblas==NULL || _g_cblas_sger==NULL) {
        return;
    }
    _g_cblas_sger(
        order,
        M,
        N,
        alpha,
        X,
        incX,
        Y,
        incY,
        A,
        lda    
    );
}
typedef void (CALLBACK* PFNcblas_dger)( /* cblas_dger */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST double            /* alpha */,
    OPENBLAS_CONST double *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST double *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */,
    double *            /* A */,
    OPENBLAS_CONST blasint            /* lda */
);
static PFNcblas_dger _g_cblas_dger = NULL;
void cblas_dger(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST double            alpha,
    OPENBLAS_CONST double *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST double *            Y,
    OPENBLAS_CONST blasint            incY,
    double *            A,
    OPENBLAS_CONST blasint            lda
)
{
    if(_h_openblas==NULL || _g_cblas_dger==NULL) {
        return;
    }
    _g_cblas_dger(
        order,
        M,
        N,
        alpha,
        X,
        incX,
        Y,
        incY,
        A,
        lda    
    );
}
typedef void (CALLBACK* PFNcblas_cgeru)( /* cblas_cgeru */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST void *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */,
    void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */
);
static PFNcblas_cgeru _g_cblas_cgeru = NULL;
void cblas_cgeru(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST void *            Y,
    OPENBLAS_CONST blasint            incY,
    void *            A,
    OPENBLAS_CONST blasint            lda
)
{
    if(_h_openblas==NULL || _g_cblas_cgeru==NULL) {
        return;
    }
    _g_cblas_cgeru(
        order,
        M,
        N,
        alpha,
        X,
        incX,
        Y,
        incY,
        A,
        lda    
    );
}
typedef void (CALLBACK* PFNcblas_cgerc)( /* cblas_cgerc */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST void *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */,
    void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */
);
static PFNcblas_cgerc _g_cblas_cgerc = NULL;
void cblas_cgerc(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST void *            Y,
    OPENBLAS_CONST blasint            incY,
    void *            A,
    OPENBLAS_CONST blasint            lda
)
{
    if(_h_openblas==NULL || _g_cblas_cgerc==NULL) {
        return;
    }
    _g_cblas_cgerc(
        order,
        M,
        N,
        alpha,
        X,
        incX,
        Y,
        incY,
        A,
        lda    
    );
}
typedef void (CALLBACK* PFNcblas_zgeru)( /* cblas_zgeru */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST void *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */,
    void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */
);
static PFNcblas_zgeru _g_cblas_zgeru = NULL;
void cblas_zgeru(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST void *            Y,
    OPENBLAS_CONST blasint            incY,
    void *            A,
    OPENBLAS_CONST blasint            lda
)
{
    if(_h_openblas==NULL || _g_cblas_zgeru==NULL) {
        return;
    }
    _g_cblas_zgeru(
        order,
        M,
        N,
        alpha,
        X,
        incX,
        Y,
        incY,
        A,
        lda    
    );
}
typedef void (CALLBACK* PFNcblas_zgerc)( /* cblas_zgerc */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST void *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */,
    void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */
);
static PFNcblas_zgerc _g_cblas_zgerc = NULL;
void cblas_zgerc(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST void *            Y,
    OPENBLAS_CONST blasint            incY,
    void *            A,
    OPENBLAS_CONST blasint            lda
)
{
    if(_h_openblas==NULL || _g_cblas_zgerc==NULL) {
        return;
    }
    _g_cblas_zgerc(
        order,
        M,
        N,
        alpha,
        X,
        incX,
        Y,
        incY,
        A,
        lda    
    );
}
typedef void (CALLBACK* PFNcblas_strsv)( /* cblas_strsv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST float *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    float *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_strsv _g_cblas_strsv = NULL;
void cblas_strsv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST float *            A,
    OPENBLAS_CONST blasint            lda,
    float *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_strsv==NULL) {
        return;
    }
    _g_cblas_strsv(
        order,
        Uplo,
        TransA,
        Diag,
        N,
        A,
        lda,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_dtrsv)( /* cblas_dtrsv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST double *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    double *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_dtrsv _g_cblas_dtrsv = NULL;
void cblas_dtrsv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST double *            A,
    OPENBLAS_CONST blasint            lda,
    double *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_dtrsv==NULL) {
        return;
    }
    _g_cblas_dtrsv(
        order,
        Uplo,
        TransA,
        Diag,
        N,
        A,
        lda,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_ctrsv)( /* cblas_ctrsv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_ctrsv _g_cblas_ctrsv = NULL;
void cblas_ctrsv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    void *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_ctrsv==NULL) {
        return;
    }
    _g_cblas_ctrsv(
        order,
        Uplo,
        TransA,
        Diag,
        N,
        A,
        lda,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_ztrsv)( /* cblas_ztrsv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_ztrsv _g_cblas_ztrsv = NULL;
void cblas_ztrsv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    void *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_ztrsv==NULL) {
        return;
    }
    _g_cblas_ztrsv(
        order,
        Uplo,
        TransA,
        Diag,
        N,
        A,
        lda,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_strmv)( /* cblas_strmv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST float *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    float *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_strmv _g_cblas_strmv = NULL;
void cblas_strmv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST float *            A,
    OPENBLAS_CONST blasint            lda,
    float *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_strmv==NULL) {
        return;
    }
    _g_cblas_strmv(
        order,
        Uplo,
        TransA,
        Diag,
        N,
        A,
        lda,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_dtrmv)( /* cblas_dtrmv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST double *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    double *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_dtrmv _g_cblas_dtrmv = NULL;
void cblas_dtrmv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST double *            A,
    OPENBLAS_CONST blasint            lda,
    double *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_dtrmv==NULL) {
        return;
    }
    _g_cblas_dtrmv(
        order,
        Uplo,
        TransA,
        Diag,
        N,
        A,
        lda,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_ctrmv)( /* cblas_ctrmv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_ctrmv _g_cblas_ctrmv = NULL;
void cblas_ctrmv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    void *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_ctrmv==NULL) {
        return;
    }
    _g_cblas_ctrmv(
        order,
        Uplo,
        TransA,
        Diag,
        N,
        A,
        lda,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_ztrmv)( /* cblas_ztrmv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_ztrmv _g_cblas_ztrmv = NULL;
void cblas_ztrmv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    void *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_ztrmv==NULL) {
        return;
    }
    _g_cblas_ztrmv(
        order,
        Uplo,
        TransA,
        Diag,
        N,
        A,
        lda,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_ssyr)( /* cblas_ssyr */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST float            /* alpha */,
    OPENBLAS_CONST float *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    float *            /* A */,
    OPENBLAS_CONST blasint            /* lda */
);
static PFNcblas_ssyr _g_cblas_ssyr = NULL;
void cblas_ssyr(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST float            alpha,
    OPENBLAS_CONST float *            X,
    OPENBLAS_CONST blasint            incX,
    float *            A,
    OPENBLAS_CONST blasint            lda
)
{
    if(_h_openblas==NULL || _g_cblas_ssyr==NULL) {
        return;
    }
    _g_cblas_ssyr(
        order,
        Uplo,
        N,
        alpha,
        X,
        incX,
        A,
        lda    
    );
}
typedef void (CALLBACK* PFNcblas_dsyr)( /* cblas_dsyr */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST double            /* alpha */,
    OPENBLAS_CONST double *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    double *            /* A */,
    OPENBLAS_CONST blasint            /* lda */
);
static PFNcblas_dsyr _g_cblas_dsyr = NULL;
void cblas_dsyr(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST double            alpha,
    OPENBLAS_CONST double *            X,
    OPENBLAS_CONST blasint            incX,
    double *            A,
    OPENBLAS_CONST blasint            lda
)
{
    if(_h_openblas==NULL || _g_cblas_dsyr==NULL) {
        return;
    }
    _g_cblas_dsyr(
        order,
        Uplo,
        N,
        alpha,
        X,
        incX,
        A,
        lda    
    );
}
typedef void (CALLBACK* PFNcblas_cher)( /* cblas_cher */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST float            /* alpha */,
    OPENBLAS_CONST void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */
);
static PFNcblas_cher _g_cblas_cher = NULL;
void cblas_cher(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST float            alpha,
    OPENBLAS_CONST void *            X,
    OPENBLAS_CONST blasint            incX,
    void *            A,
    OPENBLAS_CONST blasint            lda
)
{
    if(_h_openblas==NULL || _g_cblas_cher==NULL) {
        return;
    }
    _g_cblas_cher(
        order,
        Uplo,
        N,
        alpha,
        X,
        incX,
        A,
        lda    
    );
}
typedef void (CALLBACK* PFNcblas_zher)( /* cblas_zher */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST double            /* alpha */,
    OPENBLAS_CONST void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */
);
static PFNcblas_zher _g_cblas_zher = NULL;
void cblas_zher(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST double            alpha,
    OPENBLAS_CONST void *            X,
    OPENBLAS_CONST blasint            incX,
    void *            A,
    OPENBLAS_CONST blasint            lda
)
{
    if(_h_openblas==NULL || _g_cblas_zher==NULL) {
        return;
    }
    _g_cblas_zher(
        order,
        Uplo,
        N,
        alpha,
        X,
        incX,
        A,
        lda    
    );
}
typedef void (CALLBACK* PFNcblas_ssyr2)( /* cblas_ssyr2 */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST float            /* alpha */,
    OPENBLAS_CONST float *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST float *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */,
    float *            /* A */,
    OPENBLAS_CONST blasint            /* lda */
);
static PFNcblas_ssyr2 _g_cblas_ssyr2 = NULL;
void cblas_ssyr2(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST float            alpha,
    OPENBLAS_CONST float *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST float *            Y,
    OPENBLAS_CONST blasint            incY,
    float *            A,
    OPENBLAS_CONST blasint            lda
)
{
    if(_h_openblas==NULL || _g_cblas_ssyr2==NULL) {
        return;
    }
    _g_cblas_ssyr2(
        order,
        Uplo,
        N,
        alpha,
        X,
        incX,
        Y,
        incY,
        A,
        lda    
    );
}
typedef void (CALLBACK* PFNcblas_dsyr2)( /* cblas_dsyr2 */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST double            /* alpha */,
    OPENBLAS_CONST double *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST double *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */,
    double *            /* A */,
    OPENBLAS_CONST blasint            /* lda */
);
static PFNcblas_dsyr2 _g_cblas_dsyr2 = NULL;
void cblas_dsyr2(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST double            alpha,
    OPENBLAS_CONST double *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST double *            Y,
    OPENBLAS_CONST blasint            incY,
    double *            A,
    OPENBLAS_CONST blasint            lda
)
{
    if(_h_openblas==NULL || _g_cblas_dsyr2==NULL) {
        return;
    }
    _g_cblas_dsyr2(
        order,
        Uplo,
        N,
        alpha,
        X,
        incX,
        Y,
        incY,
        A,
        lda    
    );
}
typedef void (CALLBACK* PFNcblas_cher2)( /* cblas_cher2 */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST void *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */,
    void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */
);
static PFNcblas_cher2 _g_cblas_cher2 = NULL;
void cblas_cher2(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST void *            Y,
    OPENBLAS_CONST blasint            incY,
    void *            A,
    OPENBLAS_CONST blasint            lda
)
{
    if(_h_openblas==NULL || _g_cblas_cher2==NULL) {
        return;
    }
    _g_cblas_cher2(
        order,
        Uplo,
        N,
        alpha,
        X,
        incX,
        Y,
        incY,
        A,
        lda    
    );
}
typedef void (CALLBACK* PFNcblas_zher2)( /* cblas_zher2 */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST void *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */,
    void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */
);
static PFNcblas_zher2 _g_cblas_zher2 = NULL;
void cblas_zher2(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST void *            Y,
    OPENBLAS_CONST blasint            incY,
    void *            A,
    OPENBLAS_CONST blasint            lda
)
{
    if(_h_openblas==NULL || _g_cblas_zher2==NULL) {
        return;
    }
    _g_cblas_zher2(
        order,
        Uplo,
        N,
        alpha,
        X,
        incX,
        Y,
        incY,
        A,
        lda    
    );
}
typedef void (CALLBACK* PFNcblas_sgbmv)( /* cblas_sgbmv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* KL */,
    OPENBLAS_CONST blasint            /* KU */,
    OPENBLAS_CONST float            /* alpha */,
    OPENBLAS_CONST float *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST float *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST float            /* beta */,
    float *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */
);
static PFNcblas_sgbmv _g_cblas_sgbmv = NULL;
void cblas_sgbmv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            KL,
    OPENBLAS_CONST blasint            KU,
    OPENBLAS_CONST float            alpha,
    OPENBLAS_CONST float *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST float *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST float            beta,
    float *            Y,
    OPENBLAS_CONST blasint            incY
)
{
    if(_h_openblas==NULL || _g_cblas_sgbmv==NULL) {
        return;
    }
    _g_cblas_sgbmv(
        order,
        TransA,
        M,
        N,
        KL,
        KU,
        alpha,
        A,
        lda,
        X,
        incX,
        beta,
        Y,
        incY    
    );
}
typedef void (CALLBACK* PFNcblas_dgbmv)( /* cblas_dgbmv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* KL */,
    OPENBLAS_CONST blasint            /* KU */,
    OPENBLAS_CONST double            /* alpha */,
    OPENBLAS_CONST double *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST double *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST double            /* beta */,
    double *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */
);
static PFNcblas_dgbmv _g_cblas_dgbmv = NULL;
void cblas_dgbmv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            KL,
    OPENBLAS_CONST blasint            KU,
    OPENBLAS_CONST double            alpha,
    OPENBLAS_CONST double *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST double *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST double            beta,
    double *            Y,
    OPENBLAS_CONST blasint            incY
)
{
    if(_h_openblas==NULL || _g_cblas_dgbmv==NULL) {
        return;
    }
    _g_cblas_dgbmv(
        order,
        TransA,
        M,
        N,
        KL,
        KU,
        alpha,
        A,
        lda,
        X,
        incX,
        beta,
        Y,
        incY    
    );
}
typedef void (CALLBACK* PFNcblas_cgbmv)( /* cblas_cgbmv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* KL */,
    OPENBLAS_CONST blasint            /* KU */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST void *            /* beta */,
    void *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */
);
static PFNcblas_cgbmv _g_cblas_cgbmv = NULL;
void cblas_cgbmv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            KL,
    OPENBLAS_CONST blasint            KU,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST void *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST void *            beta,
    void *            Y,
    OPENBLAS_CONST blasint            incY
)
{
    if(_h_openblas==NULL || _g_cblas_cgbmv==NULL) {
        return;
    }
    _g_cblas_cgbmv(
        order,
        TransA,
        M,
        N,
        KL,
        KU,
        alpha,
        A,
        lda,
        X,
        incX,
        beta,
        Y,
        incY    
    );
}
typedef void (CALLBACK* PFNcblas_zgbmv)( /* cblas_zgbmv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* KL */,
    OPENBLAS_CONST blasint            /* KU */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST void *            /* beta */,
    void *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */
);
static PFNcblas_zgbmv _g_cblas_zgbmv = NULL;
void cblas_zgbmv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            KL,
    OPENBLAS_CONST blasint            KU,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST void *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST void *            beta,
    void *            Y,
    OPENBLAS_CONST blasint            incY
)
{
    if(_h_openblas==NULL || _g_cblas_zgbmv==NULL) {
        return;
    }
    _g_cblas_zgbmv(
        order,
        TransA,
        M,
        N,
        KL,
        KU,
        alpha,
        A,
        lda,
        X,
        incX,
        beta,
        Y,
        incY    
    );
}
typedef void (CALLBACK* PFNcblas_ssbmv)( /* cblas_ssbmv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST float            /* alpha */,
    OPENBLAS_CONST float *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST float *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST float            /* beta */,
    float *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */
);
static PFNcblas_ssbmv _g_cblas_ssbmv = NULL;
void cblas_ssbmv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST float            alpha,
    OPENBLAS_CONST float *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST float *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST float            beta,
    float *            Y,
    OPENBLAS_CONST blasint            incY
)
{
    if(_h_openblas==NULL || _g_cblas_ssbmv==NULL) {
        return;
    }
    _g_cblas_ssbmv(
        order,
        Uplo,
        N,
        K,
        alpha,
        A,
        lda,
        X,
        incX,
        beta,
        Y,
        incY    
    );
}
typedef void (CALLBACK* PFNcblas_dsbmv)( /* cblas_dsbmv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST double            /* alpha */,
    OPENBLAS_CONST double *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST double *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST double            /* beta */,
    double *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */
);
static PFNcblas_dsbmv _g_cblas_dsbmv = NULL;
void cblas_dsbmv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST double            alpha,
    OPENBLAS_CONST double *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST double *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST double            beta,
    double *            Y,
    OPENBLAS_CONST blasint            incY
)
{
    if(_h_openblas==NULL || _g_cblas_dsbmv==NULL) {
        return;
    }
    _g_cblas_dsbmv(
        order,
        Uplo,
        N,
        K,
        alpha,
        A,
        lda,
        X,
        incX,
        beta,
        Y,
        incY    
    );
}
typedef void (CALLBACK* PFNcblas_stbmv)( /* cblas_stbmv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST float *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    float *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_stbmv _g_cblas_stbmv = NULL;
void cblas_stbmv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST float *            A,
    OPENBLAS_CONST blasint            lda,
    float *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_stbmv==NULL) {
        return;
    }
    _g_cblas_stbmv(
        order,
        Uplo,
        TransA,
        Diag,
        N,
        K,
        A,
        lda,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_dtbmv)( /* cblas_dtbmv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST double *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    double *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_dtbmv _g_cblas_dtbmv = NULL;
void cblas_dtbmv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST double *            A,
    OPENBLAS_CONST blasint            lda,
    double *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_dtbmv==NULL) {
        return;
    }
    _g_cblas_dtbmv(
        order,
        Uplo,
        TransA,
        Diag,
        N,
        K,
        A,
        lda,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_ctbmv)( /* cblas_ctbmv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_ctbmv _g_cblas_ctbmv = NULL;
void cblas_ctbmv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    void *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_ctbmv==NULL) {
        return;
    }
    _g_cblas_ctbmv(
        order,
        Uplo,
        TransA,
        Diag,
        N,
        K,
        A,
        lda,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_ztbmv)( /* cblas_ztbmv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_ztbmv _g_cblas_ztbmv = NULL;
void cblas_ztbmv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    void *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_ztbmv==NULL) {
        return;
    }
    _g_cblas_ztbmv(
        order,
        Uplo,
        TransA,
        Diag,
        N,
        K,
        A,
        lda,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_stbsv)( /* cblas_stbsv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST float *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    float *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_stbsv _g_cblas_stbsv = NULL;
void cblas_stbsv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST float *            A,
    OPENBLAS_CONST blasint            lda,
    float *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_stbsv==NULL) {
        return;
    }
    _g_cblas_stbsv(
        order,
        Uplo,
        TransA,
        Diag,
        N,
        K,
        A,
        lda,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_dtbsv)( /* cblas_dtbsv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST double *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    double *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_dtbsv _g_cblas_dtbsv = NULL;
void cblas_dtbsv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST double *            A,
    OPENBLAS_CONST blasint            lda,
    double *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_dtbsv==NULL) {
        return;
    }
    _g_cblas_dtbsv(
        order,
        Uplo,
        TransA,
        Diag,
        N,
        K,
        A,
        lda,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_ctbsv)( /* cblas_ctbsv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_ctbsv _g_cblas_ctbsv = NULL;
void cblas_ctbsv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    void *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_ctbsv==NULL) {
        return;
    }
    _g_cblas_ctbsv(
        order,
        Uplo,
        TransA,
        Diag,
        N,
        K,
        A,
        lda,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_ztbsv)( /* cblas_ztbsv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_ztbsv _g_cblas_ztbsv = NULL;
void cblas_ztbsv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    void *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_ztbsv==NULL) {
        return;
    }
    _g_cblas_ztbsv(
        order,
        Uplo,
        TransA,
        Diag,
        N,
        K,
        A,
        lda,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_stpmv)( /* cblas_stpmv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST float *            /* Ap */,
    float *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_stpmv _g_cblas_stpmv = NULL;
void cblas_stpmv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST float *            Ap,
    float *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_stpmv==NULL) {
        return;
    }
    _g_cblas_stpmv(
        order,
        Uplo,
        TransA,
        Diag,
        N,
        Ap,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_dtpmv)( /* cblas_dtpmv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST double *            /* Ap */,
    double *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_dtpmv _g_cblas_dtpmv = NULL;
void cblas_dtpmv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST double *            Ap,
    double *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_dtpmv==NULL) {
        return;
    }
    _g_cblas_dtpmv(
        order,
        Uplo,
        TransA,
        Diag,
        N,
        Ap,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_ctpmv)( /* cblas_ctpmv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* Ap */,
    void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_ctpmv _g_cblas_ctpmv = NULL;
void cblas_ctpmv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            Ap,
    void *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_ctpmv==NULL) {
        return;
    }
    _g_cblas_ctpmv(
        order,
        Uplo,
        TransA,
        Diag,
        N,
        Ap,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_ztpmv)( /* cblas_ztpmv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* Ap */,
    void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_ztpmv _g_cblas_ztpmv = NULL;
void cblas_ztpmv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            Ap,
    void *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_ztpmv==NULL) {
        return;
    }
    _g_cblas_ztpmv(
        order,
        Uplo,
        TransA,
        Diag,
        N,
        Ap,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_stpsv)( /* cblas_stpsv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST float *            /* Ap */,
    float *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_stpsv _g_cblas_stpsv = NULL;
void cblas_stpsv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST float *            Ap,
    float *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_stpsv==NULL) {
        return;
    }
    _g_cblas_stpsv(
        order,
        Uplo,
        TransA,
        Diag,
        N,
        Ap,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_dtpsv)( /* cblas_dtpsv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST double *            /* Ap */,
    double *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_dtpsv _g_cblas_dtpsv = NULL;
void cblas_dtpsv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST double *            Ap,
    double *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_dtpsv==NULL) {
        return;
    }
    _g_cblas_dtpsv(
        order,
        Uplo,
        TransA,
        Diag,
        N,
        Ap,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_ctpsv)( /* cblas_ctpsv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* Ap */,
    void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_ctpsv _g_cblas_ctpsv = NULL;
void cblas_ctpsv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            Ap,
    void *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_ctpsv==NULL) {
        return;
    }
    _g_cblas_ctpsv(
        order,
        Uplo,
        TransA,
        Diag,
        N,
        Ap,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_ztpsv)( /* cblas_ztpsv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* Ap */,
    void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */
);
static PFNcblas_ztpsv _g_cblas_ztpsv = NULL;
void cblas_ztpsv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            Ap,
    void *            X,
    OPENBLAS_CONST blasint            incX
)
{
    if(_h_openblas==NULL || _g_cblas_ztpsv==NULL) {
        return;
    }
    _g_cblas_ztpsv(
        order,
        Uplo,
        TransA,
        Diag,
        N,
        Ap,
        X,
        incX    
    );
}
typedef void (CALLBACK* PFNcblas_ssymv)( /* cblas_ssymv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST float            /* alpha */,
    OPENBLAS_CONST float *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST float *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST float            /* beta */,
    float *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */
);
static PFNcblas_ssymv _g_cblas_ssymv = NULL;
void cblas_ssymv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST float            alpha,
    OPENBLAS_CONST float *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST float *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST float            beta,
    float *            Y,
    OPENBLAS_CONST blasint            incY
)
{
    if(_h_openblas==NULL || _g_cblas_ssymv==NULL) {
        return;
    }
    _g_cblas_ssymv(
        order,
        Uplo,
        N,
        alpha,
        A,
        lda,
        X,
        incX,
        beta,
        Y,
        incY    
    );
}
typedef void (CALLBACK* PFNcblas_dsymv)( /* cblas_dsymv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST double            /* alpha */,
    OPENBLAS_CONST double *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST double *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST double            /* beta */,
    double *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */
);
static PFNcblas_dsymv _g_cblas_dsymv = NULL;
void cblas_dsymv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST double            alpha,
    OPENBLAS_CONST double *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST double *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST double            beta,
    double *            Y,
    OPENBLAS_CONST blasint            incY
)
{
    if(_h_openblas==NULL || _g_cblas_dsymv==NULL) {
        return;
    }
    _g_cblas_dsymv(
        order,
        Uplo,
        N,
        alpha,
        A,
        lda,
        X,
        incX,
        beta,
        Y,
        incY    
    );
}
typedef void (CALLBACK* PFNcblas_chemv)( /* cblas_chemv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST void *            /* beta */,
    void *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */
);
static PFNcblas_chemv _g_cblas_chemv = NULL;
void cblas_chemv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST void *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST void *            beta,
    void *            Y,
    OPENBLAS_CONST blasint            incY
)
{
    if(_h_openblas==NULL || _g_cblas_chemv==NULL) {
        return;
    }
    _g_cblas_chemv(
        order,
        Uplo,
        N,
        alpha,
        A,
        lda,
        X,
        incX,
        beta,
        Y,
        incY    
    );
}
typedef void (CALLBACK* PFNcblas_zhemv)( /* cblas_zhemv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST void *            /* beta */,
    void *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */
);
static PFNcblas_zhemv _g_cblas_zhemv = NULL;
void cblas_zhemv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST void *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST void *            beta,
    void *            Y,
    OPENBLAS_CONST blasint            incY
)
{
    if(_h_openblas==NULL || _g_cblas_zhemv==NULL) {
        return;
    }
    _g_cblas_zhemv(
        order,
        Uplo,
        N,
        alpha,
        A,
        lda,
        X,
        incX,
        beta,
        Y,
        incY    
    );
}
typedef void (CALLBACK* PFNcblas_sspmv)( /* cblas_sspmv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST float            /* alpha */,
    OPENBLAS_CONST float *            /* Ap */,
    OPENBLAS_CONST float *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST float            /* beta */,
    float *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */
);
static PFNcblas_sspmv _g_cblas_sspmv = NULL;
void cblas_sspmv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST float            alpha,
    OPENBLAS_CONST float *            Ap,
    OPENBLAS_CONST float *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST float            beta,
    float *            Y,
    OPENBLAS_CONST blasint            incY
)
{
    if(_h_openblas==NULL || _g_cblas_sspmv==NULL) {
        return;
    }
    _g_cblas_sspmv(
        order,
        Uplo,
        N,
        alpha,
        Ap,
        X,
        incX,
        beta,
        Y,
        incY    
    );
}
typedef void (CALLBACK* PFNcblas_dspmv)( /* cblas_dspmv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST double            /* alpha */,
    OPENBLAS_CONST double *            /* Ap */,
    OPENBLAS_CONST double *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST double            /* beta */,
    double *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */
);
static PFNcblas_dspmv _g_cblas_dspmv = NULL;
void cblas_dspmv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST double            alpha,
    OPENBLAS_CONST double *            Ap,
    OPENBLAS_CONST double *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST double            beta,
    double *            Y,
    OPENBLAS_CONST blasint            incY
)
{
    if(_h_openblas==NULL || _g_cblas_dspmv==NULL) {
        return;
    }
    _g_cblas_dspmv(
        order,
        Uplo,
        N,
        alpha,
        Ap,
        X,
        incX,
        beta,
        Y,
        incY    
    );
}
typedef void (CALLBACK* PFNcblas_sspr)( /* cblas_sspr */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST float            /* alpha */,
    OPENBLAS_CONST float *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    float *            /* Ap */
);
static PFNcblas_sspr _g_cblas_sspr = NULL;
void cblas_sspr(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST float            alpha,
    OPENBLAS_CONST float *            X,
    OPENBLAS_CONST blasint            incX,
    float *            Ap
)
{
    if(_h_openblas==NULL || _g_cblas_sspr==NULL) {
        return;
    }
    _g_cblas_sspr(
        order,
        Uplo,
        N,
        alpha,
        X,
        incX,
        Ap    
    );
}
typedef void (CALLBACK* PFNcblas_dspr)( /* cblas_dspr */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST double            /* alpha */,
    OPENBLAS_CONST double *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    double *            /* Ap */
);
static PFNcblas_dspr _g_cblas_dspr = NULL;
void cblas_dspr(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST double            alpha,
    OPENBLAS_CONST double *            X,
    OPENBLAS_CONST blasint            incX,
    double *            Ap
)
{
    if(_h_openblas==NULL || _g_cblas_dspr==NULL) {
        return;
    }
    _g_cblas_dspr(
        order,
        Uplo,
        N,
        alpha,
        X,
        incX,
        Ap    
    );
}
typedef void (CALLBACK* PFNcblas_chpr)( /* cblas_chpr */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST float            /* alpha */,
    OPENBLAS_CONST void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    void *            /* A */
);
static PFNcblas_chpr _g_cblas_chpr = NULL;
void cblas_chpr(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST float            alpha,
    OPENBLAS_CONST void *            X,
    OPENBLAS_CONST blasint            incX,
    void *            A
)
{
    if(_h_openblas==NULL || _g_cblas_chpr==NULL) {
        return;
    }
    _g_cblas_chpr(
        order,
        Uplo,
        N,
        alpha,
        X,
        incX,
        A    
    );
}
typedef void (CALLBACK* PFNcblas_zhpr)( /* cblas_zhpr */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST double            /* alpha */,
    OPENBLAS_CONST void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    void *            /* A */
);
static PFNcblas_zhpr _g_cblas_zhpr = NULL;
void cblas_zhpr(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST double            alpha,
    OPENBLAS_CONST void *            X,
    OPENBLAS_CONST blasint            incX,
    void *            A
)
{
    if(_h_openblas==NULL || _g_cblas_zhpr==NULL) {
        return;
    }
    _g_cblas_zhpr(
        order,
        Uplo,
        N,
        alpha,
        X,
        incX,
        A    
    );
}
typedef void (CALLBACK* PFNcblas_sspr2)( /* cblas_sspr2 */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST float            /* alpha */,
    OPENBLAS_CONST float *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST float *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */,
    float *            /* A */
);
static PFNcblas_sspr2 _g_cblas_sspr2 = NULL;
void cblas_sspr2(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST float            alpha,
    OPENBLAS_CONST float *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST float *            Y,
    OPENBLAS_CONST blasint            incY,
    float *            A
)
{
    if(_h_openblas==NULL || _g_cblas_sspr2==NULL) {
        return;
    }
    _g_cblas_sspr2(
        order,
        Uplo,
        N,
        alpha,
        X,
        incX,
        Y,
        incY,
        A    
    );
}
typedef void (CALLBACK* PFNcblas_dspr2)( /* cblas_dspr2 */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST double            /* alpha */,
    OPENBLAS_CONST double *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST double *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */,
    double *            /* A */
);
static PFNcblas_dspr2 _g_cblas_dspr2 = NULL;
void cblas_dspr2(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST double            alpha,
    OPENBLAS_CONST double *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST double *            Y,
    OPENBLAS_CONST blasint            incY,
    double *            A
)
{
    if(_h_openblas==NULL || _g_cblas_dspr2==NULL) {
        return;
    }
    _g_cblas_dspr2(
        order,
        Uplo,
        N,
        alpha,
        X,
        incX,
        Y,
        incY,
        A    
    );
}
typedef void (CALLBACK* PFNcblas_chpr2)( /* cblas_chpr2 */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST void *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */,
    void *            /* Ap */
);
static PFNcblas_chpr2 _g_cblas_chpr2 = NULL;
void cblas_chpr2(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST void *            Y,
    OPENBLAS_CONST blasint            incY,
    void *            Ap
)
{
    if(_h_openblas==NULL || _g_cblas_chpr2==NULL) {
        return;
    }
    _g_cblas_chpr2(
        order,
        Uplo,
        N,
        alpha,
        X,
        incX,
        Y,
        incY,
        Ap    
    );
}
typedef void (CALLBACK* PFNcblas_zhpr2)( /* cblas_zhpr2 */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST void *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */,
    void *            /* Ap */
);
static PFNcblas_zhpr2 _g_cblas_zhpr2 = NULL;
void cblas_zhpr2(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST void *            Y,
    OPENBLAS_CONST blasint            incY,
    void *            Ap
)
{
    if(_h_openblas==NULL || _g_cblas_zhpr2==NULL) {
        return;
    }
    _g_cblas_zhpr2(
        order,
        Uplo,
        N,
        alpha,
        X,
        incX,
        Y,
        incY,
        Ap    
    );
}
typedef void (CALLBACK* PFNcblas_chbmv)( /* cblas_chbmv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST void *            /* beta */,
    void *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */
);
static PFNcblas_chbmv _g_cblas_chbmv = NULL;
void cblas_chbmv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST void *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST void *            beta,
    void *            Y,
    OPENBLAS_CONST blasint            incY
)
{
    if(_h_openblas==NULL || _g_cblas_chbmv==NULL) {
        return;
    }
    _g_cblas_chbmv(
        order,
        Uplo,
        N,
        K,
        alpha,
        A,
        lda,
        X,
        incX,
        beta,
        Y,
        incY    
    );
}
typedef void (CALLBACK* PFNcblas_zhbmv)( /* cblas_zhbmv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST void *            /* beta */,
    void *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */
);
static PFNcblas_zhbmv _g_cblas_zhbmv = NULL;
void cblas_zhbmv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST void *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST void *            beta,
    void *            Y,
    OPENBLAS_CONST blasint            incY
)
{
    if(_h_openblas==NULL || _g_cblas_zhbmv==NULL) {
        return;
    }
    _g_cblas_zhbmv(
        order,
        Uplo,
        N,
        K,
        alpha,
        A,
        lda,
        X,
        incX,
        beta,
        Y,
        incY    
    );
}
typedef void (CALLBACK* PFNcblas_chpmv)( /* cblas_chpmv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* Ap */,
    OPENBLAS_CONST void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST void *            /* beta */,
    void *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */
);
static PFNcblas_chpmv _g_cblas_chpmv = NULL;
void cblas_chpmv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            Ap,
    OPENBLAS_CONST void *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST void *            beta,
    void *            Y,
    OPENBLAS_CONST blasint            incY
)
{
    if(_h_openblas==NULL || _g_cblas_chpmv==NULL) {
        return;
    }
    _g_cblas_chpmv(
        order,
        Uplo,
        N,
        alpha,
        Ap,
        X,
        incX,
        beta,
        Y,
        incY    
    );
}
typedef void (CALLBACK* PFNcblas_zhpmv)( /* cblas_zhpmv */
    OPENBLAS_CONST enum CBLAS_ORDER            /* order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* Ap */,
    OPENBLAS_CONST void *            /* X */,
    OPENBLAS_CONST blasint            /* incX */,
    OPENBLAS_CONST void *            /* beta */,
    void *            /* Y */,
    OPENBLAS_CONST blasint            /* incY */
);
static PFNcblas_zhpmv _g_cblas_zhpmv = NULL;
void cblas_zhpmv(
    OPENBLAS_CONST enum CBLAS_ORDER            order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            Ap,
    OPENBLAS_CONST void *            X,
    OPENBLAS_CONST blasint            incX,
    OPENBLAS_CONST void *            beta,
    void *            Y,
    OPENBLAS_CONST blasint            incY
)
{
    if(_h_openblas==NULL || _g_cblas_zhpmv==NULL) {
        return;
    }
    _g_cblas_zhpmv(
        order,
        Uplo,
        N,
        alpha,
        Ap,
        X,
        incX,
        beta,
        Y,
        incY    
    );
}
typedef void (CALLBACK* PFNcblas_sgemm)( /* cblas_sgemm */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransB */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST float            /* alpha */,
    OPENBLAS_CONST float *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST float *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */,
    OPENBLAS_CONST float            /* beta */,
    float *            /* C */,
    OPENBLAS_CONST blasint            /* ldc */
);
static PFNcblas_sgemm _g_cblas_sgemm = NULL;
void cblas_sgemm(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransB,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST float            alpha,
    OPENBLAS_CONST float *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST float *            B,
    OPENBLAS_CONST blasint            ldb,
    OPENBLAS_CONST float            beta,
    float *            C,
    OPENBLAS_CONST blasint            ldc
)
{
    if(_h_openblas==NULL || _g_cblas_sgemm==NULL) {
        return;
    }
    _g_cblas_sgemm(
        Order,
        TransA,
        TransB,
        M,
        N,
        K,
        alpha,
        A,
        lda,
        B,
        ldb,
        beta,
        C,
        ldc    
    );
}
typedef void (CALLBACK* PFNcblas_dgemm)( /* cblas_dgemm */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransB */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST double            /* alpha */,
    OPENBLAS_CONST double *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST double *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */,
    OPENBLAS_CONST double            /* beta */,
    double *            /* C */,
    OPENBLAS_CONST blasint            /* ldc */
);
static PFNcblas_dgemm _g_cblas_dgemm = NULL;
void cblas_dgemm(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransB,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST double            alpha,
    OPENBLAS_CONST double *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST double *            B,
    OPENBLAS_CONST blasint            ldb,
    OPENBLAS_CONST double            beta,
    double *            C,
    OPENBLAS_CONST blasint            ldc
)
{
    if(_h_openblas==NULL || _g_cblas_dgemm==NULL) {
        return;
    }
    _g_cblas_dgemm(
        Order,
        TransA,
        TransB,
        M,
        N,
        K,
        alpha,
        A,
        lda,
        B,
        ldb,
        beta,
        C,
        ldc    
    );
}
typedef void (CALLBACK* PFNcblas_cgemm)( /* cblas_cgemm */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransB */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST void *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */,
    OPENBLAS_CONST void *            /* beta */,
    void *            /* C */,
    OPENBLAS_CONST blasint            /* ldc */
);
static PFNcblas_cgemm _g_cblas_cgemm = NULL;
void cblas_cgemm(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransB,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST void *            B,
    OPENBLAS_CONST blasint            ldb,
    OPENBLAS_CONST void *            beta,
    void *            C,
    OPENBLAS_CONST blasint            ldc
)
{
    if(_h_openblas==NULL || _g_cblas_cgemm==NULL) {
        return;
    }
    _g_cblas_cgemm(
        Order,
        TransA,
        TransB,
        M,
        N,
        K,
        alpha,
        A,
        lda,
        B,
        ldb,
        beta,
        C,
        ldc    
    );
}
typedef void (CALLBACK* PFNcblas_cgemm3m)( /* cblas_cgemm3m */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransB */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST void *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */,
    OPENBLAS_CONST void *            /* beta */,
    void *            /* C */,
    OPENBLAS_CONST blasint            /* ldc */
);
static PFNcblas_cgemm3m _g_cblas_cgemm3m = NULL;
void cblas_cgemm3m(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransB,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST void *            B,
    OPENBLAS_CONST blasint            ldb,
    OPENBLAS_CONST void *            beta,
    void *            C,
    OPENBLAS_CONST blasint            ldc
)
{
    if(_h_openblas==NULL || _g_cblas_cgemm3m==NULL) {
        return;
    }
    _g_cblas_cgemm3m(
        Order,
        TransA,
        TransB,
        M,
        N,
        K,
        alpha,
        A,
        lda,
        B,
        ldb,
        beta,
        C,
        ldc    
    );
}
typedef void (CALLBACK* PFNcblas_zgemm)( /* cblas_zgemm */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransB */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST void *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */,
    OPENBLAS_CONST void *            /* beta */,
    void *            /* C */,
    OPENBLAS_CONST blasint            /* ldc */
);
static PFNcblas_zgemm _g_cblas_zgemm = NULL;
void cblas_zgemm(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransB,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST void *            B,
    OPENBLAS_CONST blasint            ldb,
    OPENBLAS_CONST void *            beta,
    void *            C,
    OPENBLAS_CONST blasint            ldc
)
{
    if(_h_openblas==NULL || _g_cblas_zgemm==NULL) {
        return;
    }
    _g_cblas_zgemm(
        Order,
        TransA,
        TransB,
        M,
        N,
        K,
        alpha,
        A,
        lda,
        B,
        ldb,
        beta,
        C,
        ldc    
    );
}
typedef void (CALLBACK* PFNcblas_zgemm3m)( /* cblas_zgemm3m */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransB */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST void *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */,
    OPENBLAS_CONST void *            /* beta */,
    void *            /* C */,
    OPENBLAS_CONST blasint            /* ldc */
);
static PFNcblas_zgemm3m _g_cblas_zgemm3m = NULL;
void cblas_zgemm3m(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransB,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST void *            B,
    OPENBLAS_CONST blasint            ldb,
    OPENBLAS_CONST void *            beta,
    void *            C,
    OPENBLAS_CONST blasint            ldc
)
{
    if(_h_openblas==NULL || _g_cblas_zgemm3m==NULL) {
        return;
    }
    _g_cblas_zgemm3m(
        Order,
        TransA,
        TransB,
        M,
        N,
        K,
        alpha,
        A,
        lda,
        B,
        ldb,
        beta,
        C,
        ldc    
    );
}
typedef void (CALLBACK* PFNcblas_ssymm)( /* cblas_ssymm */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_SIDE            /* Side */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST float            /* alpha */,
    OPENBLAS_CONST float *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST float *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */,
    OPENBLAS_CONST float            /* beta */,
    float *            /* C */,
    OPENBLAS_CONST blasint            /* ldc */
);
static PFNcblas_ssymm _g_cblas_ssymm = NULL;
void cblas_ssymm(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_SIDE            Side,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST float            alpha,
    OPENBLAS_CONST float *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST float *            B,
    OPENBLAS_CONST blasint            ldb,
    OPENBLAS_CONST float            beta,
    float *            C,
    OPENBLAS_CONST blasint            ldc
)
{
    if(_h_openblas==NULL || _g_cblas_ssymm==NULL) {
        return;
    }
    _g_cblas_ssymm(
        Order,
        Side,
        Uplo,
        M,
        N,
        alpha,
        A,
        lda,
        B,
        ldb,
        beta,
        C,
        ldc    
    );
}
typedef void (CALLBACK* PFNcblas_dsymm)( /* cblas_dsymm */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_SIDE            /* Side */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST double            /* alpha */,
    OPENBLAS_CONST double *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST double *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */,
    OPENBLAS_CONST double            /* beta */,
    double *            /* C */,
    OPENBLAS_CONST blasint            /* ldc */
);
static PFNcblas_dsymm _g_cblas_dsymm = NULL;
void cblas_dsymm(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_SIDE            Side,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST double            alpha,
    OPENBLAS_CONST double *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST double *            B,
    OPENBLAS_CONST blasint            ldb,
    OPENBLAS_CONST double            beta,
    double *            C,
    OPENBLAS_CONST blasint            ldc
)
{
    if(_h_openblas==NULL || _g_cblas_dsymm==NULL) {
        return;
    }
    _g_cblas_dsymm(
        Order,
        Side,
        Uplo,
        M,
        N,
        alpha,
        A,
        lda,
        B,
        ldb,
        beta,
        C,
        ldc    
    );
}
typedef void (CALLBACK* PFNcblas_csymm)( /* cblas_csymm */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_SIDE            /* Side */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST void *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */,
    OPENBLAS_CONST void *            /* beta */,
    void *            /* C */,
    OPENBLAS_CONST blasint            /* ldc */
);
static PFNcblas_csymm _g_cblas_csymm = NULL;
void cblas_csymm(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_SIDE            Side,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST void *            B,
    OPENBLAS_CONST blasint            ldb,
    OPENBLAS_CONST void *            beta,
    void *            C,
    OPENBLAS_CONST blasint            ldc
)
{
    if(_h_openblas==NULL || _g_cblas_csymm==NULL) {
        return;
    }
    _g_cblas_csymm(
        Order,
        Side,
        Uplo,
        M,
        N,
        alpha,
        A,
        lda,
        B,
        ldb,
        beta,
        C,
        ldc    
    );
}
typedef void (CALLBACK* PFNcblas_zsymm)( /* cblas_zsymm */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_SIDE            /* Side */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST void *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */,
    OPENBLAS_CONST void *            /* beta */,
    void *            /* C */,
    OPENBLAS_CONST blasint            /* ldc */
);
static PFNcblas_zsymm _g_cblas_zsymm = NULL;
void cblas_zsymm(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_SIDE            Side,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST void *            B,
    OPENBLAS_CONST blasint            ldb,
    OPENBLAS_CONST void *            beta,
    void *            C,
    OPENBLAS_CONST blasint            ldc
)
{
    if(_h_openblas==NULL || _g_cblas_zsymm==NULL) {
        return;
    }
    _g_cblas_zsymm(
        Order,
        Side,
        Uplo,
        M,
        N,
        alpha,
        A,
        lda,
        B,
        ldb,
        beta,
        C,
        ldc    
    );
}
typedef void (CALLBACK* PFNcblas_ssyrk)( /* cblas_ssyrk */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* Trans */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST float            /* alpha */,
    OPENBLAS_CONST float *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST float            /* beta */,
    float *            /* C */,
    OPENBLAS_CONST blasint            /* ldc */
);
static PFNcblas_ssyrk _g_cblas_ssyrk = NULL;
void cblas_ssyrk(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            Trans,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST float            alpha,
    OPENBLAS_CONST float *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST float            beta,
    float *            C,
    OPENBLAS_CONST blasint            ldc
)
{
    if(_h_openblas==NULL || _g_cblas_ssyrk==NULL) {
        return;
    }
    _g_cblas_ssyrk(
        Order,
        Uplo,
        Trans,
        N,
        K,
        alpha,
        A,
        lda,
        beta,
        C,
        ldc    
    );
}
typedef void (CALLBACK* PFNcblas_dsyrk)( /* cblas_dsyrk */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* Trans */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST double            /* alpha */,
    OPENBLAS_CONST double *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST double            /* beta */,
    double *            /* C */,
    OPENBLAS_CONST blasint            /* ldc */
);
static PFNcblas_dsyrk _g_cblas_dsyrk = NULL;
void cblas_dsyrk(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            Trans,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST double            alpha,
    OPENBLAS_CONST double *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST double            beta,
    double *            C,
    OPENBLAS_CONST blasint            ldc
)
{
    if(_h_openblas==NULL || _g_cblas_dsyrk==NULL) {
        return;
    }
    _g_cblas_dsyrk(
        Order,
        Uplo,
        Trans,
        N,
        K,
        alpha,
        A,
        lda,
        beta,
        C,
        ldc    
    );
}
typedef void (CALLBACK* PFNcblas_csyrk)( /* cblas_csyrk */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* Trans */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST void *            /* beta */,
    void *            /* C */,
    OPENBLAS_CONST blasint            /* ldc */
);
static PFNcblas_csyrk _g_cblas_csyrk = NULL;
void cblas_csyrk(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            Trans,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST void *            beta,
    void *            C,
    OPENBLAS_CONST blasint            ldc
)
{
    if(_h_openblas==NULL || _g_cblas_csyrk==NULL) {
        return;
    }
    _g_cblas_csyrk(
        Order,
        Uplo,
        Trans,
        N,
        K,
        alpha,
        A,
        lda,
        beta,
        C,
        ldc    
    );
}
typedef void (CALLBACK* PFNcblas_zsyrk)( /* cblas_zsyrk */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* Trans */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST void *            /* beta */,
    void *            /* C */,
    OPENBLAS_CONST blasint            /* ldc */
);
static PFNcblas_zsyrk _g_cblas_zsyrk = NULL;
void cblas_zsyrk(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            Trans,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST void *            beta,
    void *            C,
    OPENBLAS_CONST blasint            ldc
)
{
    if(_h_openblas==NULL || _g_cblas_zsyrk==NULL) {
        return;
    }
    _g_cblas_zsyrk(
        Order,
        Uplo,
        Trans,
        N,
        K,
        alpha,
        A,
        lda,
        beta,
        C,
        ldc    
    );
}
typedef void (CALLBACK* PFNcblas_ssyr2k)( /* cblas_ssyr2k */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* Trans */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST float            /* alpha */,
    OPENBLAS_CONST float *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST float *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */,
    OPENBLAS_CONST float            /* beta */,
    float *            /* C */,
    OPENBLAS_CONST blasint            /* ldc */
);
static PFNcblas_ssyr2k _g_cblas_ssyr2k = NULL;
void cblas_ssyr2k(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            Trans,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST float            alpha,
    OPENBLAS_CONST float *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST float *            B,
    OPENBLAS_CONST blasint            ldb,
    OPENBLAS_CONST float            beta,
    float *            C,
    OPENBLAS_CONST blasint            ldc
)
{
    if(_h_openblas==NULL || _g_cblas_ssyr2k==NULL) {
        return;
    }
    _g_cblas_ssyr2k(
        Order,
        Uplo,
        Trans,
        N,
        K,
        alpha,
        A,
        lda,
        B,
        ldb,
        beta,
        C,
        ldc    
    );
}
typedef void (CALLBACK* PFNcblas_dsyr2k)( /* cblas_dsyr2k */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* Trans */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST double            /* alpha */,
    OPENBLAS_CONST double *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST double *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */,
    OPENBLAS_CONST double            /* beta */,
    double *            /* C */,
    OPENBLAS_CONST blasint            /* ldc */
);
static PFNcblas_dsyr2k _g_cblas_dsyr2k = NULL;
void cblas_dsyr2k(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            Trans,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST double            alpha,
    OPENBLAS_CONST double *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST double *            B,
    OPENBLAS_CONST blasint            ldb,
    OPENBLAS_CONST double            beta,
    double *            C,
    OPENBLAS_CONST blasint            ldc
)
{
    if(_h_openblas==NULL || _g_cblas_dsyr2k==NULL) {
        return;
    }
    _g_cblas_dsyr2k(
        Order,
        Uplo,
        Trans,
        N,
        K,
        alpha,
        A,
        lda,
        B,
        ldb,
        beta,
        C,
        ldc    
    );
}
typedef void (CALLBACK* PFNcblas_csyr2k)( /* cblas_csyr2k */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* Trans */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST void *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */,
    OPENBLAS_CONST void *            /* beta */,
    void *            /* C */,
    OPENBLAS_CONST blasint            /* ldc */
);
static PFNcblas_csyr2k _g_cblas_csyr2k = NULL;
void cblas_csyr2k(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            Trans,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST void *            B,
    OPENBLAS_CONST blasint            ldb,
    OPENBLAS_CONST void *            beta,
    void *            C,
    OPENBLAS_CONST blasint            ldc
)
{
    if(_h_openblas==NULL || _g_cblas_csyr2k==NULL) {
        return;
    }
    _g_cblas_csyr2k(
        Order,
        Uplo,
        Trans,
        N,
        K,
        alpha,
        A,
        lda,
        B,
        ldb,
        beta,
        C,
        ldc    
    );
}
typedef void (CALLBACK* PFNcblas_zsyr2k)( /* cblas_zsyr2k */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* Trans */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST void *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */,
    OPENBLAS_CONST void *            /* beta */,
    void *            /* C */,
    OPENBLAS_CONST blasint            /* ldc */
);
static PFNcblas_zsyr2k _g_cblas_zsyr2k = NULL;
void cblas_zsyr2k(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            Trans,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST void *            B,
    OPENBLAS_CONST blasint            ldb,
    OPENBLAS_CONST void *            beta,
    void *            C,
    OPENBLAS_CONST blasint            ldc
)
{
    if(_h_openblas==NULL || _g_cblas_zsyr2k==NULL) {
        return;
    }
    _g_cblas_zsyr2k(
        Order,
        Uplo,
        Trans,
        N,
        K,
        alpha,
        A,
        lda,
        B,
        ldb,
        beta,
        C,
        ldc    
    );
}
typedef void (CALLBACK* PFNcblas_strmm)( /* cblas_strmm */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_SIDE            /* Side */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST float            /* alpha */,
    OPENBLAS_CONST float *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    float *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */
);
static PFNcblas_strmm _g_cblas_strmm = NULL;
void cblas_strmm(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_SIDE            Side,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST float            alpha,
    OPENBLAS_CONST float *            A,
    OPENBLAS_CONST blasint            lda,
    float *            B,
    OPENBLAS_CONST blasint            ldb
)
{
    if(_h_openblas==NULL || _g_cblas_strmm==NULL) {
        return;
    }
    _g_cblas_strmm(
        Order,
        Side,
        Uplo,
        TransA,
        Diag,
        M,
        N,
        alpha,
        A,
        lda,
        B,
        ldb    
    );
}
typedef void (CALLBACK* PFNcblas_dtrmm)( /* cblas_dtrmm */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_SIDE            /* Side */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST double            /* alpha */,
    OPENBLAS_CONST double *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    double *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */
);
static PFNcblas_dtrmm _g_cblas_dtrmm = NULL;
void cblas_dtrmm(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_SIDE            Side,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST double            alpha,
    OPENBLAS_CONST double *            A,
    OPENBLAS_CONST blasint            lda,
    double *            B,
    OPENBLAS_CONST blasint            ldb
)
{
    if(_h_openblas==NULL || _g_cblas_dtrmm==NULL) {
        return;
    }
    _g_cblas_dtrmm(
        Order,
        Side,
        Uplo,
        TransA,
        Diag,
        M,
        N,
        alpha,
        A,
        lda,
        B,
        ldb    
    );
}
typedef void (CALLBACK* PFNcblas_ctrmm)( /* cblas_ctrmm */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_SIDE            /* Side */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    void *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */
);
static PFNcblas_ctrmm _g_cblas_ctrmm = NULL;
void cblas_ctrmm(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_SIDE            Side,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    void *            B,
    OPENBLAS_CONST blasint            ldb
)
{
    if(_h_openblas==NULL || _g_cblas_ctrmm==NULL) {
        return;
    }
    _g_cblas_ctrmm(
        Order,
        Side,
        Uplo,
        TransA,
        Diag,
        M,
        N,
        alpha,
        A,
        lda,
        B,
        ldb    
    );
}
typedef void (CALLBACK* PFNcblas_ztrmm)( /* cblas_ztrmm */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_SIDE            /* Side */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    void *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */
);
static PFNcblas_ztrmm _g_cblas_ztrmm = NULL;
void cblas_ztrmm(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_SIDE            Side,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    void *            B,
    OPENBLAS_CONST blasint            ldb
)
{
    if(_h_openblas==NULL || _g_cblas_ztrmm==NULL) {
        return;
    }
    _g_cblas_ztrmm(
        Order,
        Side,
        Uplo,
        TransA,
        Diag,
        M,
        N,
        alpha,
        A,
        lda,
        B,
        ldb    
    );
}
typedef void (CALLBACK* PFNcblas_strsm)( /* cblas_strsm */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_SIDE            /* Side */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST float            /* alpha */,
    OPENBLAS_CONST float *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    float *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */
);
static PFNcblas_strsm _g_cblas_strsm = NULL;
void cblas_strsm(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_SIDE            Side,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST float            alpha,
    OPENBLAS_CONST float *            A,
    OPENBLAS_CONST blasint            lda,
    float *            B,
    OPENBLAS_CONST blasint            ldb
)
{
    if(_h_openblas==NULL || _g_cblas_strsm==NULL) {
        return;
    }
    _g_cblas_strsm(
        Order,
        Side,
        Uplo,
        TransA,
        Diag,
        M,
        N,
        alpha,
        A,
        lda,
        B,
        ldb    
    );
}
typedef void (CALLBACK* PFNcblas_dtrsm)( /* cblas_dtrsm */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_SIDE            /* Side */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST double            /* alpha */,
    OPENBLAS_CONST double *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    double *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */
);
static PFNcblas_dtrsm _g_cblas_dtrsm = NULL;
void cblas_dtrsm(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_SIDE            Side,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST double            alpha,
    OPENBLAS_CONST double *            A,
    OPENBLAS_CONST blasint            lda,
    double *            B,
    OPENBLAS_CONST blasint            ldb
)
{
    if(_h_openblas==NULL || _g_cblas_dtrsm==NULL) {
        return;
    }
    _g_cblas_dtrsm(
        Order,
        Side,
        Uplo,
        TransA,
        Diag,
        M,
        N,
        alpha,
        A,
        lda,
        B,
        ldb    
    );
}
typedef void (CALLBACK* PFNcblas_ctrsm)( /* cblas_ctrsm */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_SIDE            /* Side */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    void *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */
);
static PFNcblas_ctrsm _g_cblas_ctrsm = NULL;
void cblas_ctrsm(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_SIDE            Side,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    void *            B,
    OPENBLAS_CONST blasint            ldb
)
{
    if(_h_openblas==NULL || _g_cblas_ctrsm==NULL) {
        return;
    }
    _g_cblas_ctrsm(
        Order,
        Side,
        Uplo,
        TransA,
        Diag,
        M,
        N,
        alpha,
        A,
        lda,
        B,
        ldb    
    );
}
typedef void (CALLBACK* PFNcblas_ztrsm)( /* cblas_ztrsm */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_SIDE            /* Side */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* TransA */,
    OPENBLAS_CONST enum CBLAS_DIAG            /* Diag */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    void *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */
);
static PFNcblas_ztrsm _g_cblas_ztrsm = NULL;
void cblas_ztrsm(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_SIDE            Side,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            TransA,
    OPENBLAS_CONST enum CBLAS_DIAG            Diag,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    void *            B,
    OPENBLAS_CONST blasint            ldb
)
{
    if(_h_openblas==NULL || _g_cblas_ztrsm==NULL) {
        return;
    }
    _g_cblas_ztrsm(
        Order,
        Side,
        Uplo,
        TransA,
        Diag,
        M,
        N,
        alpha,
        A,
        lda,
        B,
        ldb    
    );
}
typedef void (CALLBACK* PFNcblas_chemm)( /* cblas_chemm */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_SIDE            /* Side */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST void *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */,
    OPENBLAS_CONST void *            /* beta */,
    void *            /* C */,
    OPENBLAS_CONST blasint            /* ldc */
);
static PFNcblas_chemm _g_cblas_chemm = NULL;
void cblas_chemm(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_SIDE            Side,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST void *            B,
    OPENBLAS_CONST blasint            ldb,
    OPENBLAS_CONST void *            beta,
    void *            C,
    OPENBLAS_CONST blasint            ldc
)
{
    if(_h_openblas==NULL || _g_cblas_chemm==NULL) {
        return;
    }
    _g_cblas_chemm(
        Order,
        Side,
        Uplo,
        M,
        N,
        alpha,
        A,
        lda,
        B,
        ldb,
        beta,
        C,
        ldc    
    );
}
typedef void (CALLBACK* PFNcblas_zhemm)( /* cblas_zhemm */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_SIDE            /* Side */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST blasint            /* M */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST void *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */,
    OPENBLAS_CONST void *            /* beta */,
    void *            /* C */,
    OPENBLAS_CONST blasint            /* ldc */
);
static PFNcblas_zhemm _g_cblas_zhemm = NULL;
void cblas_zhemm(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_SIDE            Side,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST blasint            M,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST void *            B,
    OPENBLAS_CONST blasint            ldb,
    OPENBLAS_CONST void *            beta,
    void *            C,
    OPENBLAS_CONST blasint            ldc
)
{
    if(_h_openblas==NULL || _g_cblas_zhemm==NULL) {
        return;
    }
    _g_cblas_zhemm(
        Order,
        Side,
        Uplo,
        M,
        N,
        alpha,
        A,
        lda,
        B,
        ldb,
        beta,
        C,
        ldc    
    );
}
typedef void (CALLBACK* PFNcblas_cherk)( /* cblas_cherk */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* Trans */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST float            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST float            /* beta */,
    void *            /* C */,
    OPENBLAS_CONST blasint            /* ldc */
);
static PFNcblas_cherk _g_cblas_cherk = NULL;
void cblas_cherk(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            Trans,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST float            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST float            beta,
    void *            C,
    OPENBLAS_CONST blasint            ldc
)
{
    if(_h_openblas==NULL || _g_cblas_cherk==NULL) {
        return;
    }
    _g_cblas_cherk(
        Order,
        Uplo,
        Trans,
        N,
        K,
        alpha,
        A,
        lda,
        beta,
        C,
        ldc    
    );
}
typedef void (CALLBACK* PFNcblas_zherk)( /* cblas_zherk */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* Trans */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST double            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST double            /* beta */,
    void *            /* C */,
    OPENBLAS_CONST blasint            /* ldc */
);
static PFNcblas_zherk _g_cblas_zherk = NULL;
void cblas_zherk(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            Trans,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST double            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST double            beta,
    void *            C,
    OPENBLAS_CONST blasint            ldc
)
{
    if(_h_openblas==NULL || _g_cblas_zherk==NULL) {
        return;
    }
    _g_cblas_zherk(
        Order,
        Uplo,
        Trans,
        N,
        K,
        alpha,
        A,
        lda,
        beta,
        C,
        ldc    
    );
}
typedef void (CALLBACK* PFNcblas_cher2k)( /* cblas_cher2k */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* Trans */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST void *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */,
    OPENBLAS_CONST float            /* beta */,
    void *            /* C */,
    OPENBLAS_CONST blasint            /* ldc */
);
static PFNcblas_cher2k _g_cblas_cher2k = NULL;
void cblas_cher2k(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            Trans,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST void *            B,
    OPENBLAS_CONST blasint            ldb,
    OPENBLAS_CONST float            beta,
    void *            C,
    OPENBLAS_CONST blasint            ldc
)
{
    if(_h_openblas==NULL || _g_cblas_cher2k==NULL) {
        return;
    }
    _g_cblas_cher2k(
        Order,
        Uplo,
        Trans,
        N,
        K,
        alpha,
        A,
        lda,
        B,
        ldb,
        beta,
        C,
        ldc    
    );
}
typedef void (CALLBACK* PFNcblas_zher2k)( /* cblas_zher2k */
    OPENBLAS_CONST enum CBLAS_ORDER            /* Order */,
    OPENBLAS_CONST enum CBLAS_UPLO            /* Uplo */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* Trans */,
    OPENBLAS_CONST blasint            /* N */,
    OPENBLAS_CONST blasint            /* K */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* A */,
    OPENBLAS_CONST blasint            /* lda */,
    OPENBLAS_CONST void *            /* B */,
    OPENBLAS_CONST blasint            /* ldb */,
    OPENBLAS_CONST double            /* beta */,
    void *            /* C */,
    OPENBLAS_CONST blasint            /* ldc */
);
static PFNcblas_zher2k _g_cblas_zher2k = NULL;
void cblas_zher2k(
    OPENBLAS_CONST enum CBLAS_ORDER            Order,
    OPENBLAS_CONST enum CBLAS_UPLO            Uplo,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            Trans,
    OPENBLAS_CONST blasint            N,
    OPENBLAS_CONST blasint            K,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            A,
    OPENBLAS_CONST blasint            lda,
    OPENBLAS_CONST void *            B,
    OPENBLAS_CONST blasint            ldb,
    OPENBLAS_CONST double            beta,
    void *            C,
    OPENBLAS_CONST blasint            ldc
)
{
    if(_h_openblas==NULL || _g_cblas_zher2k==NULL) {
        return;
    }
    _g_cblas_zher2k(
        Order,
        Uplo,
        Trans,
        N,
        K,
        alpha,
        A,
        lda,
        B,
        ldb,
        beta,
        C,
        ldc    
    );
}
typedef void (CALLBACK* PFNcblas_saxpby)( /* cblas_saxpby */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST float            /* alpha */,
    OPENBLAS_CONST float *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    OPENBLAS_CONST float            /* beta */,
    float *            /* y */,
    OPENBLAS_CONST blasint            /* incy */
);
static PFNcblas_saxpby _g_cblas_saxpby = NULL;
void cblas_saxpby(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST float            alpha,
    OPENBLAS_CONST float *            x,
    OPENBLAS_CONST blasint            incx,
    OPENBLAS_CONST float            beta,
    float *            y,
    OPENBLAS_CONST blasint            incy
)
{
    if(_h_openblas==NULL || _g_cblas_saxpby==NULL) {
        return;
    }
    _g_cblas_saxpby(
        n,
        alpha,
        x,
        incx,
        beta,
        y,
        incy    
    );
}
typedef void (CALLBACK* PFNcblas_daxpby)( /* cblas_daxpby */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST double            /* alpha */,
    OPENBLAS_CONST double *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    OPENBLAS_CONST double            /* beta */,
    double *            /* y */,
    OPENBLAS_CONST blasint            /* incy */
);
static PFNcblas_daxpby _g_cblas_daxpby = NULL;
void cblas_daxpby(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST double            alpha,
    OPENBLAS_CONST double *            x,
    OPENBLAS_CONST blasint            incx,
    OPENBLAS_CONST double            beta,
    double *            y,
    OPENBLAS_CONST blasint            incy
)
{
    if(_h_openblas==NULL || _g_cblas_daxpby==NULL) {
        return;
    }
    _g_cblas_daxpby(
        n,
        alpha,
        x,
        incx,
        beta,
        y,
        incy    
    );
}
typedef void (CALLBACK* PFNcblas_caxpby)( /* cblas_caxpby */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    OPENBLAS_CONST void *            /* beta */,
    void *            /* y */,
    OPENBLAS_CONST blasint            /* incy */
);
static PFNcblas_caxpby _g_cblas_caxpby = NULL;
void cblas_caxpby(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            x,
    OPENBLAS_CONST blasint            incx,
    OPENBLAS_CONST void *            beta,
    void *            y,
    OPENBLAS_CONST blasint            incy
)
{
    if(_h_openblas==NULL || _g_cblas_caxpby==NULL) {
        return;
    }
    _g_cblas_caxpby(
        n,
        alpha,
        x,
        incx,
        beta,
        y,
        incy    
    );
}
typedef void (CALLBACK* PFNcblas_zaxpby)( /* cblas_zaxpby */
    OPENBLAS_CONST blasint            /* n */,
    OPENBLAS_CONST void *            /* alpha */,
    OPENBLAS_CONST void *            /* x */,
    OPENBLAS_CONST blasint            /* incx */,
    OPENBLAS_CONST void *            /* beta */,
    void *            /* y */,
    OPENBLAS_CONST blasint            /* incy */
);
static PFNcblas_zaxpby _g_cblas_zaxpby = NULL;
void cblas_zaxpby(
    OPENBLAS_CONST blasint            n,
    OPENBLAS_CONST void *            alpha,
    OPENBLAS_CONST void *            x,
    OPENBLAS_CONST blasint            incx,
    OPENBLAS_CONST void *            beta,
    void *            y,
    OPENBLAS_CONST blasint            incy
)
{
    if(_h_openblas==NULL || _g_cblas_zaxpby==NULL) {
        return;
    }
    _g_cblas_zaxpby(
        n,
        alpha,
        x,
        incx,
        beta,
        y,
        incy    
    );
}
typedef void (CALLBACK* PFNcblas_somatcopy)( /* cblas_somatcopy */
    OPENBLAS_CONST enum CBLAS_ORDER            /* CORDER */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* CTRANS */,
    OPENBLAS_CONST blasint            /* crows */,
    OPENBLAS_CONST blasint            /* ccols */,
    OPENBLAS_CONST float            /* calpha */,
    OPENBLAS_CONST float *            /* a */,
    OPENBLAS_CONST blasint            /* clda */,
    float *            /* b */,
    OPENBLAS_CONST blasint            /* cldb */
);
static PFNcblas_somatcopy _g_cblas_somatcopy = NULL;
void cblas_somatcopy(
    OPENBLAS_CONST enum CBLAS_ORDER            CORDER,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            CTRANS,
    OPENBLAS_CONST blasint            crows,
    OPENBLAS_CONST blasint            ccols,
    OPENBLAS_CONST float            calpha,
    OPENBLAS_CONST float *            a,
    OPENBLAS_CONST blasint            clda,
    float *            b,
    OPENBLAS_CONST blasint            cldb
)
{
    if(_h_openblas==NULL || _g_cblas_somatcopy==NULL) {
        return;
    }
    _g_cblas_somatcopy(
        CORDER,
        CTRANS,
        crows,
        ccols,
        calpha,
        a,
        clda,
        b,
        cldb    
    );
}
typedef void (CALLBACK* PFNcblas_domatcopy)( /* cblas_domatcopy */
    OPENBLAS_CONST enum CBLAS_ORDER            /* CORDER */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* CTRANS */,
    OPENBLAS_CONST blasint            /* crows */,
    OPENBLAS_CONST blasint            /* ccols */,
    OPENBLAS_CONST double            /* calpha */,
    OPENBLAS_CONST double *            /* a */,
    OPENBLAS_CONST blasint            /* clda */,
    double *            /* b */,
    OPENBLAS_CONST blasint            /* cldb */
);
static PFNcblas_domatcopy _g_cblas_domatcopy = NULL;
void cblas_domatcopy(
    OPENBLAS_CONST enum CBLAS_ORDER            CORDER,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            CTRANS,
    OPENBLAS_CONST blasint            crows,
    OPENBLAS_CONST blasint            ccols,
    OPENBLAS_CONST double            calpha,
    OPENBLAS_CONST double *            a,
    OPENBLAS_CONST blasint            clda,
    double *            b,
    OPENBLAS_CONST blasint            cldb
)
{
    if(_h_openblas==NULL || _g_cblas_domatcopy==NULL) {
        return;
    }
    _g_cblas_domatcopy(
        CORDER,
        CTRANS,
        crows,
        ccols,
        calpha,
        a,
        clda,
        b,
        cldb    
    );
}
typedef void (CALLBACK* PFNcblas_comatcopy)( /* cblas_comatcopy */
    OPENBLAS_CONST enum CBLAS_ORDER            /* CORDER */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* CTRANS */,
    OPENBLAS_CONST blasint            /* crows */,
    OPENBLAS_CONST blasint            /* ccols */,
    OPENBLAS_CONST float *            /* calpha */,
    OPENBLAS_CONST float *            /* a */,
    OPENBLAS_CONST blasint            /* clda */,
    float *            /* b */,
    OPENBLAS_CONST blasint            /* cldb */
);
static PFNcblas_comatcopy _g_cblas_comatcopy = NULL;
void cblas_comatcopy(
    OPENBLAS_CONST enum CBLAS_ORDER            CORDER,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            CTRANS,
    OPENBLAS_CONST blasint            crows,
    OPENBLAS_CONST blasint            ccols,
    OPENBLAS_CONST float *            calpha,
    OPENBLAS_CONST float *            a,
    OPENBLAS_CONST blasint            clda,
    float *            b,
    OPENBLAS_CONST blasint            cldb
)
{
    if(_h_openblas==NULL || _g_cblas_comatcopy==NULL) {
        return;
    }
    _g_cblas_comatcopy(
        CORDER,
        CTRANS,
        crows,
        ccols,
        calpha,
        a,
        clda,
        b,
        cldb    
    );
}
typedef void (CALLBACK* PFNcblas_zomatcopy)( /* cblas_zomatcopy */
    OPENBLAS_CONST enum CBLAS_ORDER            /* CORDER */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* CTRANS */,
    OPENBLAS_CONST blasint            /* crows */,
    OPENBLAS_CONST blasint            /* ccols */,
    OPENBLAS_CONST double *            /* calpha */,
    OPENBLAS_CONST double *            /* a */,
    OPENBLAS_CONST blasint            /* clda */,
    double *            /* b */,
    OPENBLAS_CONST blasint            /* cldb */
);
static PFNcblas_zomatcopy _g_cblas_zomatcopy = NULL;
void cblas_zomatcopy(
    OPENBLAS_CONST enum CBLAS_ORDER            CORDER,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            CTRANS,
    OPENBLAS_CONST blasint            crows,
    OPENBLAS_CONST blasint            ccols,
    OPENBLAS_CONST double *            calpha,
    OPENBLAS_CONST double *            a,
    OPENBLAS_CONST blasint            clda,
    double *            b,
    OPENBLAS_CONST blasint            cldb
)
{
    if(_h_openblas==NULL || _g_cblas_zomatcopy==NULL) {
        return;
    }
    _g_cblas_zomatcopy(
        CORDER,
        CTRANS,
        crows,
        ccols,
        calpha,
        a,
        clda,
        b,
        cldb    
    );
}
typedef void (CALLBACK* PFNcblas_simatcopy)( /* cblas_simatcopy */
    OPENBLAS_CONST enum CBLAS_ORDER            /* CORDER */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* CTRANS */,
    OPENBLAS_CONST blasint            /* crows */,
    OPENBLAS_CONST blasint            /* ccols */,
    OPENBLAS_CONST float            /* calpha */,
    float *            /* a */,
    OPENBLAS_CONST blasint            /* clda */,
    OPENBLAS_CONST blasint            /* cldb */
);
static PFNcblas_simatcopy _g_cblas_simatcopy = NULL;
void cblas_simatcopy(
    OPENBLAS_CONST enum CBLAS_ORDER            CORDER,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            CTRANS,
    OPENBLAS_CONST blasint            crows,
    OPENBLAS_CONST blasint            ccols,
    OPENBLAS_CONST float            calpha,
    float *            a,
    OPENBLAS_CONST blasint            clda,
    OPENBLAS_CONST blasint            cldb
)
{
    if(_h_openblas==NULL || _g_cblas_simatcopy==NULL) {
        return;
    }
    _g_cblas_simatcopy(
        CORDER,
        CTRANS,
        crows,
        ccols,
        calpha,
        a,
        clda,
        cldb    
    );
}
typedef void (CALLBACK* PFNcblas_dimatcopy)( /* cblas_dimatcopy */
    OPENBLAS_CONST enum CBLAS_ORDER            /* CORDER */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* CTRANS */,
    OPENBLAS_CONST blasint            /* crows */,
    OPENBLAS_CONST blasint            /* ccols */,
    OPENBLAS_CONST double            /* calpha */,
    double *            /* a */,
    OPENBLAS_CONST blasint            /* clda */,
    OPENBLAS_CONST blasint            /* cldb */
);
static PFNcblas_dimatcopy _g_cblas_dimatcopy = NULL;
void cblas_dimatcopy(
    OPENBLAS_CONST enum CBLAS_ORDER            CORDER,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            CTRANS,
    OPENBLAS_CONST blasint            crows,
    OPENBLAS_CONST blasint            ccols,
    OPENBLAS_CONST double            calpha,
    double *            a,
    OPENBLAS_CONST blasint            clda,
    OPENBLAS_CONST blasint            cldb
)
{
    if(_h_openblas==NULL || _g_cblas_dimatcopy==NULL) {
        return;
    }
    _g_cblas_dimatcopy(
        CORDER,
        CTRANS,
        crows,
        ccols,
        calpha,
        a,
        clda,
        cldb    
    );
}
typedef void (CALLBACK* PFNcblas_cimatcopy)( /* cblas_cimatcopy */
    OPENBLAS_CONST enum CBLAS_ORDER            /* CORDER */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* CTRANS */,
    OPENBLAS_CONST blasint            /* crows */,
    OPENBLAS_CONST blasint            /* ccols */,
    OPENBLAS_CONST float *            /* calpha */,
    float *            /* a */,
    OPENBLAS_CONST blasint            /* clda */,
    OPENBLAS_CONST blasint            /* cldb */
);
static PFNcblas_cimatcopy _g_cblas_cimatcopy = NULL;
void cblas_cimatcopy(
    OPENBLAS_CONST enum CBLAS_ORDER            CORDER,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            CTRANS,
    OPENBLAS_CONST blasint            crows,
    OPENBLAS_CONST blasint            ccols,
    OPENBLAS_CONST float *            calpha,
    float *            a,
    OPENBLAS_CONST blasint            clda,
    OPENBLAS_CONST blasint            cldb
)
{
    if(_h_openblas==NULL || _g_cblas_cimatcopy==NULL) {
        return;
    }
    _g_cblas_cimatcopy(
        CORDER,
        CTRANS,
        crows,
        ccols,
        calpha,
        a,
        clda,
        cldb    
    );
}
typedef void (CALLBACK* PFNcblas_zimatcopy)( /* cblas_zimatcopy */
    OPENBLAS_CONST enum CBLAS_ORDER            /* CORDER */,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            /* CTRANS */,
    OPENBLAS_CONST blasint            /* crows */,
    OPENBLAS_CONST blasint            /* ccols */,
    OPENBLAS_CONST double *            /* calpha */,
    double *            /* a */,
    OPENBLAS_CONST blasint            /* clda */,
    OPENBLAS_CONST blasint            /* cldb */
);
static PFNcblas_zimatcopy _g_cblas_zimatcopy = NULL;
void cblas_zimatcopy(
    OPENBLAS_CONST enum CBLAS_ORDER            CORDER,
    OPENBLAS_CONST enum CBLAS_TRANSPOSE            CTRANS,
    OPENBLAS_CONST blasint            crows,
    OPENBLAS_CONST blasint            ccols,
    OPENBLAS_CONST double *            calpha,
    double *            a,
    OPENBLAS_CONST blasint            clda,
    OPENBLAS_CONST blasint            cldb
)
{
    if(_h_openblas==NULL || _g_cblas_zimatcopy==NULL) {
        return;
    }
    _g_cblas_zimatcopy(
        CORDER,
        CTRANS,
        crows,
        ccols,
        calpha,
        a,
        clda,
        cldb    
    );
}
typedef void (CALLBACK* PFNcblas_sgeadd)( /* cblas_sgeadd */
    OPENBLAS_CONST enum CBLAS_ORDER            /* CORDER */,
    OPENBLAS_CONST blasint            /* crows */,
    OPENBLAS_CONST blasint            /* ccols */,
    OPENBLAS_CONST float            /* calpha */,
    float *            /* a */,
    OPENBLAS_CONST blasint            /* clda */,
    OPENBLAS_CONST float            /* cbeta */,
    float *            /* c */,
    OPENBLAS_CONST blasint            /* cldc */
);
static PFNcblas_sgeadd _g_cblas_sgeadd = NULL;
void cblas_sgeadd(
    OPENBLAS_CONST enum CBLAS_ORDER            CORDER,
    OPENBLAS_CONST blasint            crows,
    OPENBLAS_CONST blasint            ccols,
    OPENBLAS_CONST float            calpha,
    float *            a,
    OPENBLAS_CONST blasint            clda,
    OPENBLAS_CONST float            cbeta,
    float *            c,
    OPENBLAS_CONST blasint            cldc
)
{
    if(_h_openblas==NULL || _g_cblas_sgeadd==NULL) {
        return;
    }
    _g_cblas_sgeadd(
        CORDER,
        crows,
        ccols,
        calpha,
        a,
        clda,
        cbeta,
        c,
        cldc    
    );
}
typedef void (CALLBACK* PFNcblas_dgeadd)( /* cblas_dgeadd */
    OPENBLAS_CONST enum CBLAS_ORDER            /* CORDER */,
    OPENBLAS_CONST blasint            /* crows */,
    OPENBLAS_CONST blasint            /* ccols */,
    OPENBLAS_CONST double            /* calpha */,
    double *            /* a */,
    OPENBLAS_CONST blasint            /* clda */,
    OPENBLAS_CONST double            /* cbeta */,
    double *            /* c */,
    OPENBLAS_CONST blasint            /* cldc */
);
static PFNcblas_dgeadd _g_cblas_dgeadd = NULL;
void cblas_dgeadd(
    OPENBLAS_CONST enum CBLAS_ORDER            CORDER,
    OPENBLAS_CONST blasint            crows,
    OPENBLAS_CONST blasint            ccols,
    OPENBLAS_CONST double            calpha,
    double *            a,
    OPENBLAS_CONST blasint            clda,
    OPENBLAS_CONST double            cbeta,
    double *            c,
    OPENBLAS_CONST blasint            cldc
)
{
    if(_h_openblas==NULL || _g_cblas_dgeadd==NULL) {
        return;
    }
    _g_cblas_dgeadd(
        CORDER,
        crows,
        ccols,
        calpha,
        a,
        clda,
        cbeta,
        c,
        cldc    
    );
}
typedef void (CALLBACK* PFNcblas_cgeadd)( /* cblas_cgeadd */
    OPENBLAS_CONST enum CBLAS_ORDER            /* CORDER */,
    OPENBLAS_CONST blasint            /* crows */,
    OPENBLAS_CONST blasint            /* ccols */,
    OPENBLAS_CONST float *            /* calpha */,
    float *            /* a */,
    OPENBLAS_CONST blasint            /* clda */,
    OPENBLAS_CONST float *            /* cbeta */,
    float *            /* c */,
    OPENBLAS_CONST blasint            /* cldc */
);
static PFNcblas_cgeadd _g_cblas_cgeadd = NULL;
void cblas_cgeadd(
    OPENBLAS_CONST enum CBLAS_ORDER            CORDER,
    OPENBLAS_CONST blasint            crows,
    OPENBLAS_CONST blasint            ccols,
    OPENBLAS_CONST float *            calpha,
    float *            a,
    OPENBLAS_CONST blasint            clda,
    OPENBLAS_CONST float *            cbeta,
    float *            c,
    OPENBLAS_CONST blasint            cldc
)
{
    if(_h_openblas==NULL || _g_cblas_cgeadd==NULL) {
        return;
    }
    _g_cblas_cgeadd(
        CORDER,
        crows,
        ccols,
        calpha,
        a,
        clda,
        cbeta,
        c,
        cldc    
    );
}
typedef void (CALLBACK* PFNcblas_zgeadd)( /* cblas_zgeadd */
    OPENBLAS_CONST enum CBLAS_ORDER            /* CORDER */,
    OPENBLAS_CONST blasint            /* crows */,
    OPENBLAS_CONST blasint            /* ccols */,
    OPENBLAS_CONST double *            /* calpha */,
    double *            /* a */,
    OPENBLAS_CONST blasint            /* clda */,
    OPENBLAS_CONST double *            /* cbeta */,
    double *            /* c */,
    OPENBLAS_CONST blasint            /* cldc */
);
static PFNcblas_zgeadd _g_cblas_zgeadd = NULL;
void cblas_zgeadd(
    OPENBLAS_CONST enum CBLAS_ORDER            CORDER,
    OPENBLAS_CONST blasint            crows,
    OPENBLAS_CONST blasint            ccols,
    OPENBLAS_CONST double *            calpha,
    double *            a,
    OPENBLAS_CONST blasint            clda,
    OPENBLAS_CONST double *            cbeta,
    double *            c,
    OPENBLAS_CONST blasint            cldc
)
{
    if(_h_openblas==NULL || _g_cblas_zgeadd==NULL) {
        return;
    }
    _g_cblas_zgeadd(
        CORDER,
        crows,
        ccols,
        calpha,
        a,
        clda,
        cbeta,
        c,
        cldc    
    );
}

#include "lapacke_client.c"

int rindow_load_openblas_dll()
{
    if(_h_openblas!=NULL) {
        return 0;
    }
    _h_openblas = LoadLibraryA( "libopenblas.dll" );
    if(_h_openblas==NULL) {
        printf("load error: libopenblas\n");
        return -1;
    }
    LOADFUNC(openblas_set_num_threads)
    LOADFUNC(goto_set_num_threads)
    LOADFUNC(openblas_get_num_threads)
    LOADFUNC(openblas_get_num_procs)
    LOADFUNC(openblas_get_config)
    LOADFUNC(openblas_get_corename)
    LOADFUNC(openblas_get_parallel)
    LOADFUNC(cblas_sdsdot)
    LOADFUNC(cblas_dsdot)
    LOADFUNC(cblas_sdot)
    LOADFUNC(cblas_ddot)
    LOADFUNC(cblas_cdotu_sub)
    LOADFUNC(cblas_cdotc_sub)
    LOADFUNC(cblas_zdotu_sub)
    LOADFUNC(cblas_zdotc_sub)
    LOADFUNC(cblas_sasum)
    LOADFUNC(cblas_dasum)
    LOADFUNC(cblas_scasum)
    LOADFUNC(cblas_dzasum)
    LOADFUNC(cblas_ssum)
    LOADFUNC(cblas_dsum)
    LOADFUNC(cblas_scsum)
    LOADFUNC(cblas_dzsum)
    LOADFUNC(cblas_snrm2)
    LOADFUNC(cblas_dnrm2)
    LOADFUNC(cblas_scnrm2)
    LOADFUNC(cblas_dznrm2)
    LOADFUNC(cblas_isamax)
    LOADFUNC(cblas_idamax)
    LOADFUNC(cblas_icamax)
    LOADFUNC(cblas_izamax)
    LOADFUNC(cblas_isamin)
    LOADFUNC(cblas_idamin)
    LOADFUNC(cblas_icamin)
    LOADFUNC(cblas_izamin)
    LOADFUNC(cblas_ismax)
    LOADFUNC(cblas_idmax)
    LOADFUNC(cblas_icmax)
    LOADFUNC(cblas_izmax)
    LOADFUNC(cblas_ismin)
    LOADFUNC(cblas_idmin)
    LOADFUNC(cblas_icmin)
    LOADFUNC(cblas_izmin)
    LOADFUNC(cblas_saxpy)
    LOADFUNC(cblas_daxpy)
    LOADFUNC(cblas_caxpy)
    LOADFUNC(cblas_zaxpy)
    LOADFUNC(cblas_scopy)
    LOADFUNC(cblas_dcopy)
    LOADFUNC(cblas_ccopy)
    LOADFUNC(cblas_zcopy)
    LOADFUNC(cblas_sswap)
    LOADFUNC(cblas_dswap)
    LOADFUNC(cblas_cswap)
    LOADFUNC(cblas_zswap)
    LOADFUNC(cblas_srot)
    LOADFUNC(cblas_drot)
    LOADFUNC(cblas_srotg)
    LOADFUNC(cblas_drotg)
    LOADFUNC(cblas_srotm)
    LOADFUNC(cblas_drotm)
    LOADFUNC(cblas_srotmg)
    LOADFUNC(cblas_drotmg)
    LOADFUNC(cblas_sscal)
    LOADFUNC(cblas_dscal)
    LOADFUNC(cblas_cscal)
    LOADFUNC(cblas_zscal)
    LOADFUNC(cblas_csscal)
    LOADFUNC(cblas_zdscal)
    LOADFUNC(cblas_sgemv)
    LOADFUNC(cblas_dgemv)
    LOADFUNC(cblas_cgemv)
    LOADFUNC(cblas_zgemv)
    LOADFUNC(cblas_sger)
    LOADFUNC(cblas_dger)
    LOADFUNC(cblas_cgeru)
    LOADFUNC(cblas_cgerc)
    LOADFUNC(cblas_zgeru)
    LOADFUNC(cblas_zgerc)
    LOADFUNC(cblas_strsv)
    LOADFUNC(cblas_dtrsv)
    LOADFUNC(cblas_ctrsv)
    LOADFUNC(cblas_ztrsv)
    LOADFUNC(cblas_strmv)
    LOADFUNC(cblas_dtrmv)
    LOADFUNC(cblas_ctrmv)
    LOADFUNC(cblas_ztrmv)
    LOADFUNC(cblas_ssyr)
    LOADFUNC(cblas_dsyr)
    LOADFUNC(cblas_cher)
    LOADFUNC(cblas_zher)
    LOADFUNC(cblas_ssyr2)
    LOADFUNC(cblas_dsyr2)
    LOADFUNC(cblas_cher2)
    LOADFUNC(cblas_zher2)
    LOADFUNC(cblas_sgbmv)
    LOADFUNC(cblas_dgbmv)
    LOADFUNC(cblas_cgbmv)
    LOADFUNC(cblas_zgbmv)
    LOADFUNC(cblas_ssbmv)
    LOADFUNC(cblas_dsbmv)
    LOADFUNC(cblas_stbmv)
    LOADFUNC(cblas_dtbmv)
    LOADFUNC(cblas_ctbmv)
    LOADFUNC(cblas_ztbmv)
    LOADFUNC(cblas_stbsv)
    LOADFUNC(cblas_dtbsv)
    LOADFUNC(cblas_ctbsv)
    LOADFUNC(cblas_ztbsv)
    LOADFUNC(cblas_stpmv)
    LOADFUNC(cblas_dtpmv)
    LOADFUNC(cblas_ctpmv)
    LOADFUNC(cblas_ztpmv)
    LOADFUNC(cblas_stpsv)
    LOADFUNC(cblas_dtpsv)
    LOADFUNC(cblas_ctpsv)
    LOADFUNC(cblas_ztpsv)
    LOADFUNC(cblas_ssymv)
    LOADFUNC(cblas_dsymv)
    LOADFUNC(cblas_chemv)
    LOADFUNC(cblas_zhemv)
    LOADFUNC(cblas_sspmv)
    LOADFUNC(cblas_dspmv)
    LOADFUNC(cblas_sspr)
    LOADFUNC(cblas_dspr)
    LOADFUNC(cblas_chpr)
    LOADFUNC(cblas_zhpr)
    LOADFUNC(cblas_sspr2)
    LOADFUNC(cblas_dspr2)
    LOADFUNC(cblas_chpr2)
    LOADFUNC(cblas_zhpr2)
    LOADFUNC(cblas_chbmv)
    LOADFUNC(cblas_zhbmv)
    LOADFUNC(cblas_chpmv)
    LOADFUNC(cblas_zhpmv)
    LOADFUNC(cblas_sgemm)
    LOADFUNC(cblas_dgemm)
    LOADFUNC(cblas_cgemm)
    LOADFUNC(cblas_cgemm3m)
    LOADFUNC(cblas_zgemm)
    LOADFUNC(cblas_zgemm3m)
    LOADFUNC(cblas_ssymm)
    LOADFUNC(cblas_dsymm)
    LOADFUNC(cblas_csymm)
    LOADFUNC(cblas_zsymm)
    LOADFUNC(cblas_ssyrk)
    LOADFUNC(cblas_dsyrk)
    LOADFUNC(cblas_csyrk)
    LOADFUNC(cblas_zsyrk)
    LOADFUNC(cblas_ssyr2k)
    LOADFUNC(cblas_dsyr2k)
    LOADFUNC(cblas_csyr2k)
    LOADFUNC(cblas_zsyr2k)
    LOADFUNC(cblas_strmm)
    LOADFUNC(cblas_dtrmm)
    LOADFUNC(cblas_ctrmm)
    LOADFUNC(cblas_ztrmm)
    LOADFUNC(cblas_strsm)
    LOADFUNC(cblas_dtrsm)
    LOADFUNC(cblas_ctrsm)
    LOADFUNC(cblas_ztrsm)
    LOADFUNC(cblas_chemm)
    LOADFUNC(cblas_zhemm)
    LOADFUNC(cblas_cherk)
    LOADFUNC(cblas_zherk)
    LOADFUNC(cblas_cher2k)
    LOADFUNC(cblas_zher2k)
    LOADFUNC(cblas_saxpby)
    LOADFUNC(cblas_daxpby)
    LOADFUNC(cblas_caxpby)
    LOADFUNC(cblas_zaxpby)
    LOADFUNC(cblas_somatcopy)
    LOADFUNC(cblas_domatcopy)
    LOADFUNC(cblas_comatcopy)
    LOADFUNC(cblas_zomatcopy)
    LOADFUNC(cblas_simatcopy)
    LOADFUNC(cblas_dimatcopy)
    LOADFUNC(cblas_cimatcopy)
    LOADFUNC(cblas_zimatcopy)
    LOADFUNC(cblas_sgeadd)
    LOADFUNC(cblas_dgeadd)
    LOADFUNC(cblas_cgeadd)
    LOADFUNC(cblas_zgeadd)

    LOADFUNC(LAPACKE_sgesvd)
    LOADFUNC(LAPACKE_dgesvd)

    return 0;
}
void rindow_unload_openblas_dll()
{
    FreeLibrary( _h_openblas );
    _h_openblas = NULL;
}
