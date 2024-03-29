typedef lapack_int (CALLBACK* PFNLAPACKE_sgesvd)( /* LAPACKE_sgesvd */
    int            /* matrix_layout */,
    char            /* jobu */,
    char            /* jobvt */,
    lapack_int            /* m */,
    lapack_int            /* n */,
    float *            /* a */,
    lapack_int            /* lda */,
    float *            /* s */,
    float *            /* u */,
    lapack_int            /* ldu */,
    float *            /* vt */,
    lapack_int            /* ldvt */,
    float *            /* superb */
);
static PFNLAPACKE_sgesvd _g_LAPACKE_sgesvd = NULL;
lapack_int LAPACKE_sgesvd(
    int            matrix_layout,
    char            jobu,
    char            jobvt,
    lapack_int            m,
    lapack_int            n,
    float *            a,
    lapack_int            lda,
    float *            s,
    float *            u,
    lapack_int            ldu,
    float *            vt,
    lapack_int            ldvt,
    float *            superb
)
{
    if(_h_openblas==NULL || _g_LAPACKE_sgesvd==NULL) {
        return 0;
    }
    return _g_LAPACKE_sgesvd(
        matrix_layout,
        jobu,
        jobvt,
        m,
        n,
        a,
        lda,
        s,
        u,
        ldu,
        vt,
        ldvt,
        superb    
    );
}
typedef lapack_int (CALLBACK* PFNLAPACKE_dgesvd)( /* LAPACKE_dgesvd */
    int            /* matrix_layout */,
    char            /* jobu */,
    char            /* jobvt */,
    lapack_int            /* m */,
    lapack_int            /* n */,
    double *            /* a */,
    lapack_int            /* lda */,
    double *            /* s */,
    double *            /* u */,
    lapack_int            /* ldu */,
    double *            /* vt */,
    lapack_int            /* ldvt */,
    double *            /* superb */
);
static PFNLAPACKE_dgesvd _g_LAPACKE_dgesvd = NULL;
lapack_int LAPACKE_dgesvd(
    int            matrix_layout,
    char            jobu,
    char            jobvt,
    lapack_int            m,
    lapack_int            n,
    double *            a,
    lapack_int            lda,
    double *            s,
    double *            u,
    lapack_int            ldu,
    double *            vt,
    lapack_int            ldvt,
    double *            superb
)
{
    if(_h_openblas==NULL || _g_LAPACKE_dgesvd==NULL) {
        return 0;
    }
    return _g_LAPACKE_dgesvd(
        matrix_layout,
        jobu,
        jobvt,
        m,
        n,
        a,
        lda,
        s,
        u,
        ldu,
        vt,
        ldvt,
        superb    
    );
}
