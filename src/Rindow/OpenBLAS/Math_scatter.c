/*
   A(X[i]) := Y(i)

   Method Rindow\OpenBLAS\Math::
    public function scatterAxis0(
        int $m,
        int $n,
        int $k,
        Buffer $A, int $offsetA, int $ldA,
        Buffer $X, int $offsetX, int $incX,
        Buffer $Y, int $offsetY, int $ldY,
        int $addMode
        ) : void
 {{{ */
static PHP_METHOD(Math, scatterAxis0)
{
    php_rindow_openblas_buffer_t* bufferA;
    php_rindow_openblas_buffer_t* bufferX;
    php_rindow_openblas_buffer_t* bufferY;
    zend_bool addMode;
    zend_long m;
    zend_long n;
    zend_long k;
    zval* a=NULL;
    zend_long offsetA;
    zend_long ldA;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zval* y=NULL;
    zend_long offsetY;
    zend_long ldY;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 13, 13)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_LONG(k)
        Z_PARAM_OBJECT_OF_CLASS(a,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_OBJECT_OF_CLASS(y,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetY)
        Z_PARAM_LONG(ldY)
        Z_PARAM_BOOL(addMode)
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
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_A, bufferA,m,n,offsetA,ldA)) {
        return;
    }

    // Check Buffer X
    bufferX = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(x);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, bufferX,k,offsetX,incX)) {
        return;
    }

    // Check Buffer Y
    bufferY = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(y);
    if(php_rindow_openblas_assert_matrix_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_Y, bufferY,k,n,offsetY,ldY)) {
        return;
    }

    // Check Buffer A and Y
    if(bufferA->dtype!=bufferY->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for A and Y", 0);
        return;
    }
    if(bufferX->dtype==php_rindow_openblas_dtype_bool) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Data type of BufferX must not be bool", 0);
        return;
    }

    switch (bufferA->dtype) {
        case php_rindow_openblas_dtype_float32:
            {
                float *a = &(((float *)bufferA->data)[offsetA]);
                float *y = &(((float *)bufferY->data)[offsetY]);
                zend_long i,selector;
                for(i=0; i<k; i++,y+=ldY) {
                    if(rindow_openblas_math_get_integer(
                                bufferX->dtype, bufferX->data, offsetX,incX,
                                i, &selector)) {
                        zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type of label number.", 0);
                        return;
                    }
                    if(selector<0||selector>=m) {
                        zend_throw_exception(spl_ce_RuntimeException, "Label number is out of bounds.", 0);
                        return;
                    }
                    if(n==1) {
                        if(addMode){
                            a[selector*ldA] += *y;
                        }else{
                            a[selector*ldA] = *y;
                        }
                    } else {
                        if(addMode){
                            cblas_saxpy((blasint)n,
                            1.0,
                            y,(blasint)1,
                            &(a[selector*ldA]), (blasint)1);
                        } else {
                            cblas_scopy((blasint)n,
                            y, (blasint)1,
                            &(a[selector*ldA]), (blasint)1);
                        }
                    }
                }
            }
            break;
        case php_rindow_openblas_dtype_float64:
            {
                double *a = &(((double *)bufferA->data)[offsetA]);
                double *y = &(((double *)bufferY->data)[offsetY]);
                zend_long i,selector;
                for(i=0; i<k; i++,y+=ldY) {
                    if(rindow_openblas_math_get_integer(
                                bufferX->dtype, bufferX->data, offsetX,incX,
                                i, &selector)) {
                        zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type of label number.", 0);
                        return;
                    }
                    if(selector<0||selector>=m) {
                        zend_throw_exception(spl_ce_RuntimeException, "Label number is out of bounds.", 0);
                        return;
                    }
                    if(n==1) {
                        if(addMode){
                            a[selector*ldA] += *y;
                        } else {
                            a[selector*ldA] = *y;
                        }
                    } else {
                        if(addMode){
                            cblas_daxpy((blasint)n,
                            1.0,
                            y,(blasint)1,
                            &(a[selector*ldA]), (blasint)1);
                        } else {
                            cblas_dcopy((blasint)n,
                            y, (blasint)1,
                            &(a[selector*ldA]), (blasint)1);
                        }
                    }
                }
            }
            break;
        default:
            if(!php_rindow_openblas_dtype_is_int(bufferA->dtype)) {
                if(!php_rindow_openblas_dtype_is_bool(bufferA->dtype)||addMode){
                    zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type.", 0);
                    return;
                }
            }
            {
                int valueSize;
                uint8_t *a, *y;
                valueSize = php_rindow_openblas_dtype_to_valuesize(bufferA->dtype);

                a = php_rindow_openblas_get_address(bufferA,offsetA,valueSize);
                y = php_rindow_openblas_get_address(bufferY,offsetY,valueSize);
                zend_long i,selector,incrementY,selectSizeA,copySize;
                incrementY = ldY*valueSize;
                selectSizeA = ldA*valueSize;
                copySize = n*valueSize;
                for(i=0; i<k; i++,y+=incrementY) {
                    if(rindow_openblas_math_get_integer(
                                bufferX->dtype, bufferX->data, offsetX,incX,
                                i, &selector)) {
                        zend_throw_exception(spl_ce_RuntimeException, "Unsupported data type of label number.", 0);
                        return;
                    }
                    if(selector<0||selector>=m) {
                        zend_throw_exception(spl_ce_RuntimeException, "Label number is out of bounds.", 0);
                        return;
                    }
                    if(addMode){
                        php_rindow_openblas_math_add(n,bufferA->dtype,y,1,&(a[selector*selectSizeA]),1);
                    } else {
                        memcpy(&(a[selector*selectSizeA]),y,  copySize);
                    }
                }
            }
            break;
    }
}

/*
   A(X[i]) := Y(i)

   Method Rindow\OpenBLAS\Math::
    public function scatterArxis1(
        int $m,
        int $n,
        Buffer $A, int $offsetA, int $ldA,
        Buffer $X, int $offsetX, int $incX,
        Buffer $Y, int $offsetY, int $incY,
        bool $addMode
        ) : void
 {{{ */
static PHP_METHOD(Math, scatterAxis1)
{
    php_rindow_openblas_buffer_t* bufferA;
    php_rindow_openblas_buffer_t* bufferX;
    php_rindow_openblas_buffer_t* bufferY;
    zend_bool addMode;
    zend_long m;
    zend_long n;
    zval* a=NULL;
    zend_long offsetA;
    zend_long ldA;
    zval* x=NULL;
    zend_long offsetX;
    zend_long incX;
    zval* y=NULL;
    zend_long offsetY;
    zend_long incY;

    //if (zend_parse_parameters(ZEND_NUM_ARGS(), "llOllOllOll",
    //        &m,
    //        &n,
    //        &a,php_rindow_openblas_buffer_ce,&offsetA,&ldA,
    //        &x,php_rindow_openblas_buffer_ce,&offsetX,&incX,
    //        &y,php_rindow_openblas_buffer_ce,&offsetY,&incY) == FAILURE) {
    //    zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid Arguments", 0);
    //    return;
    //}
    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 12, 12)
        Z_PARAM_LONG(m)
        Z_PARAM_LONG(n)
        Z_PARAM_OBJECT_OF_CLASS(a,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetA)
        Z_PARAM_LONG(ldA)
        Z_PARAM_OBJECT_OF_CLASS(x,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetX)
        Z_PARAM_LONG(incX)
        Z_PARAM_OBJECT_OF_CLASS(y,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(offsetY)
        Z_PARAM_LONG(incY)
        Z_PARAM_BOOL(addMode)
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
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_X, bufferX,m,offsetX,incX)) {
        return;
    }

    // Check Buffer Y
    bufferY = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(y);
    if(php_rindow_openblas_assert_vector_buffer_spec(
        PHP_RINDOW_OPENBLAS_ASSERT_Y, bufferY,m,offsetY,incY)) {
        return;
    }

    // Check Buffer A and Y
    if(bufferA->dtype!=bufferY->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch data type for A and Y", 0);
        return;
    }
    if(bufferX->dtype==php_rindow_openblas_dtype_bool) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Data type of BufferX must not be bool", 0);
        return;
    }

    switch (bufferA->dtype) {
        case php_rindow_openblas_dtype_float32:
            {
                float *a = &(((float *)bufferA->data)[offsetA]);
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
                    if(addMode){
                        a[i*ldA+selector] += y[i*incY];
                    } else {
                        a[i*ldA+selector] = y[i*incY];
                    }
                }
            }
            break;
        case php_rindow_openblas_dtype_float64:
            {
                double *a = &(((double *)bufferA->data)[offsetA]);
                double *y = &(((double *)bufferY->data)[offsetX]);
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
                    if(addMode){
                        a[i*ldA+selector] += y[i*incY];
                    } else {
                        a[i*ldA+selector] = y[i*incY];
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
