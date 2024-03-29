
/*
static inline int im2col1d_copyCell(
    zend_bool reverse,
    php_interop_polite_math_matrix_linear_buffer_t *images,
    zend_long images_pos,
    zend_long im_w,
    zend_long channels,
    zend_long channel_step,
    zend_long filter_w_step,
    zend_long vim_x,
    zend_long vfilter_w,
    zend_long dilation_w,
    php_interop_polite_math_matrix_linear_buffer_t *out,
    zend_long out_pos,
    zend_long out_filter_step,
    zend_long out_channel_step
    )
{
    zend_long vfilter_x;
    zend_long filter_w_pos;
    zend_long out_filter_pos;
    zend_long input_x;
    zend_long channel_pos;
    zend_long out_channel_pos;
    zend_long c;

    filter_w_pos = images_pos;
    out_filter_pos = out_pos;
    for(vfilter_x=0; vfilter_x<vfilter_w; vfilter_x+=dilation_w) {
        channel_pos = filter_w_pos;
        out_channel_pos = out_filter_pos;
        input_x = vim_x+vfilter_x;
        for(c=0; c<channels; c++) {
            if(input_x<0 || input_x>=im_w) {
                if(out_channel_pos<0 ||out_channel_pos>=out->size) {
                   zend_throw_exception(spl_ce_RuntimeException, "cols data out of range", 0);
                    return -1;
                }
                if(!reverse) {
                    if(images->dtype== php_interop_polite_math_matrix_dtype_float32) {
                        ((float*)(out->data))[out_channel_pos]
                            = 0;
                    } else {
                        ((double*)(out->data))[out_channel_pos]
                            = 0;
                    }
                }
            } else {
                if(channel_pos<0 ||channel_pos>=images->size) {
                   zend_throw_exception(spl_ce_RuntimeException, "images data out of range", 0);
                    return -1;
                }
                if(out_channel_pos<0 ||out_channel_pos>=out->size) {
                   zend_throw_exception(spl_ce_RuntimeException, "cols data out of range", 0);
                    return -1;
                }
                if(!reverse) {
                    if(images->dtype== php_interop_polite_math_matrix_dtype_float32) {
                        ((float*)(out->data))[out_channel_pos]
                            = ((float*)(images->data))[channel_pos];
                    } else {
                        ((double*)(out->data))[out_channel_pos]
                            = ((double*)(images->data))[channel_pos];
                    }
                } else {
                    // Sum for Back propagation
                    if(images->dtype== php_interop_polite_math_matrix_dtype_float32) {
                        ((float*)(images->data))[channel_pos]
                            += ((float*)(out->data))[out_channel_pos];
                    } else {
                        ((double*)(images->data))[channel_pos]
                            += ((double*)(out->data))[out_channel_pos];
                    }
                }
            }
            out_channel_pos += out_channel_step;
            channel_pos += channel_step;
        }
        out_filter_pos += out_filter_step;
            filter_w_pos += filter_w_step;
    }
    return 0;
}


static inline int im2col1d_execute(
    zend_bool reverse,
    php_interop_polite_math_matrix_linear_buffer_t* images,
    zend_long images_offset,
    zend_long images_size,
    zend_long batches,

    zend_long im_w,
    zend_long channels,
    zend_long filter_w,
    zend_long stride_w,
    zend_bool padding,

    zend_bool channels_first,
    zend_long dilation_w,
    zend_bool cols_channels_first,
    php_interop_polite_math_matrix_linear_buffer_t* cols,
    zend_long cols_offset,

    zend_long cols_size
    )
{
    zend_long out_w;
    zend_long stride_w_step;
    zend_long batch_step;
    zend_long channel_step;
    zend_long filter_w_step;
    zend_long out_filter_step;
    zend_long out_channel_step;
    zend_long out_cell_step;
    zend_long out_pos;
    zend_long batch_pos;
    zend_long padding_w;
    zend_long im_w_step;

    zend_long batch;
    zend_long stride_w_pos;
    zend_long vim_x;
    zend_long vim_w;
    zend_long vfilter_w;

    if((batches*im_w*channels)
        !=images_size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch images buffer size and images shape", 0);
        return -1;
    }
    out_w = ((im_w-(filter_w-1)*dilation_w-1)/stride_w)+1;
    if(out_w<=0) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid shape or parameters.", 0);
        return -1;
    }
    if(padding) {
        if((batches*
            im_w*filter_w*
            channels)!=cols_size) {
            zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch cols buffer size and images shape", 0);
            return -1;
        }
        padding_w = ((im_w-1)*stride_w-im_w+(filter_w-1)*dilation_w+1)/2;
        out_w = im_w;
    } else {
        if((batches*
            out_w*filter_w*
            channels)!=cols_size) {
            zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch cols buffer size and images shape", 0);
            return -1;
        }
        padding_w = 0;
    }
    if(channels_first) {
        im_w_step =    1;
        channel_step = im_w;
        batch_step =   im_w*channels;
    } else {
        channel_step = 1;
        im_w_step =  channels;
        batch_step = channels*im_w;
    }
    stride_w_step = im_w_step*stride_w;
    filter_w_step = im_w_step*dilation_w;

    if(cols_channels_first) {
        out_filter_step = 1;
        out_channel_step = filter_w;
    } else {
        out_filter_step = channels;
        out_channel_step = 1;
    }
    out_cell_step = filter_w*channels;

    batch_pos = images_offset-im_w_step*padding_w;
    out_pos = cols_offset;

    vim_w = out_w*stride_w;
    vfilter_w = filter_w*dilation_w;

    for(batch=0; batch<batches;batch++) {
        stride_w_pos = batch_pos;
        for(vim_x=0;vim_x<vim_w;vim_x+=stride_w) {
            int rc;
            rc = im2col1d_copyCell(
                reverse,
                images,
                stride_w_pos,
                im_w,
                channels,
                channel_step,
                filter_w_step,
                vim_x-padding_w,
                vfilter_w,
                dilation_w,
                cols,
                out_pos,
                out_filter_step,
                out_channel_step
            );
            if(rc) {
                return rc;
            }
            stride_w_pos += stride_w_step;
            out_pos += out_cell_step;
        }
        batch_pos += batch_step;
    }
    return 0;
}
*/

/*
    cols := im2col1d(images)

    Method Rindow\OpenBLAS\Math::
    public function im2col1d(
        bool $reverse,
        Buffer $images,
        int $images_offset,
        int $images_size,
        int $batches,

        int $im_w,
        int $channels,
        int $filter_w,
        int $stride_w,
        bool $padding,

        bool $channels_first,
        int $dilation_w,
        bool $cols_channels_first,
        Buffer $cols,
        int $cols_offset,

        int $cols_size
    ) : void
 {{{ */
static PHP_METHOD(Math, im2col1d)
{
    zend_bool reverse;
    php_interop_polite_math_matrix_linear_buffer_t* images;
    zend_long images_offset;
    zend_long images_size;
    zend_long batches;
    zend_long im_w;
    zend_long channels;
    zend_long filter_w;
    zend_long stride_w;
    zend_bool padding;
    zend_bool channels_first;
    zend_long dilation_w;
    zend_bool cols_channels_first;
    php_interop_polite_math_matrix_linear_buffer_t* cols;
    zend_long cols_offset;
    zend_long cols_size;

    zval* images_obj=NULL;
    zval* cols_obj=NULL;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 16, 16)
        Z_PARAM_BOOL(reverse)
        Z_PARAM_OBJECT(images_obj) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(images_offset)
        Z_PARAM_LONG(images_size)
        Z_PARAM_LONG(batches)

        Z_PARAM_LONG(im_w)
        Z_PARAM_LONG(channels)
        Z_PARAM_LONG(filter_w)
        Z_PARAM_LONG(stride_w)
        Z_PARAM_BOOL(padding)

        Z_PARAM_BOOL(channels_first)
        Z_PARAM_LONG(dilation_w)
        Z_PARAM_BOOL(cols_channels_first)
        Z_PARAM_OBJECT(cols_obj) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(cols_offset)

        Z_PARAM_LONG(cols_size)
    ZEND_PARSE_PARAMETERS_END();

    images = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(images_obj);
    if(php_rindow_openblas_assert_buffer_type(images,"images")) {
        return;
    }
    if(php_rindow_openblas_assert_buffer_size(
        images, images_offset, images_size,
        "Invalid images buffer offset or size")) {
        return;
    }
    cols = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(cols_obj);
    if(php_rindow_openblas_assert_buffer_type(cols,"cols")) {
        return;
    }
    if(php_rindow_openblas_assert_buffer_size(
        cols, cols_offset, cols_size,
        "Invalid cols buffer offset or size")) {
        return;
    }
    if(images->dtype!=php_interop_polite_math_matrix_dtype_float32 &&
        images->dtype!=php_interop_polite_math_matrix_dtype_float64) {
        zend_throw_exception(spl_ce_InvalidArgumentException,
            "Unsupported data type", 0);
        return;
    }
    // Check dtype and Buffer Y
    if(images->dtype!=cols->dtype) {
        zend_throw_exception(spl_ce_InvalidArgumentException,
            "Unmatch data type of images and cols", 0);
        return;
    }
    if(images->size<images_offset+images_size) {
        zend_throw_exception(spl_ce_InvalidArgumentException,
            "Images size is out of range", 0);
        return;
    }

    void *pDataImages = rindow_matlib_common_get_address(images->dtype, images->data,images_offset);
    void *pDataCols = rindow_matlib_common_get_address(cols->dtype, cols->data,cols_offset);

    int32_t rc = rindow_matlib_im2col1d(
        images->dtype,reverse,
        pDataImages,
        images_size,
        batches,
        im_w,
        channels,
        filter_w,
        stride_w,
        padding,
        channels_first,
        dilation_w,
        cols_channels_first,
        pDataCols,
        cols_size
    );
    switch(rc) {
        case 0: {
            break;
        }
        case RINDOW_MATLIB_E_UNSUPPORTED_DATA_TYPE: {
            zend_throw_exception(spl_ce_InvalidArgumentException, "Unsupported data type.", 0);
            return;
        }
        case RINDOW_MATLIB_E_INVALID_SHAPE_OR_PARAM: {
            zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid shape or parameters.", 0);
            return;
        }
        case RINDOW_MATLIB_E_UNMATCH_COLS_BUFFER_SIZE: {
            zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch cols buffer size and images shape", 0);
            return;
        }
        case RINDOW_MATLIB_E_UNMATCH_IMAGE_BUFFER_SIZE: {
            zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch images buffer size and images shape", 0);
            return;
        }
        case RINDOW_MATLIB_E_IMAGES_OUT_OF_RANGE: {
            zend_throw_exception(spl_ce_RuntimeException, "Images data out of range", 0);
            return;
        }
        case RINDOW_MATLIB_E_COLS_OUT_OF_RANGE: {
            zend_throw_exception(spl_ce_RuntimeException, "Cols data out of range", 0);
            return;
        }
        default: {
            zend_throw_exception_ex(spl_ce_RuntimeException, 0, "Unkown Error (%d)", rc);
            return;
        }
    }
    /*
    im2col1d_execute(
        reverse,
        images,
        images_offset,
        images_size,
        batches,

        im_w,
        channels,
        filter_w,
        stride_w,
        padding,

        channels_first,
        dilation_w,
        cols_channels_first,
        cols,
        cols_offset,

        cols_size
    );
    */
    return;
}
/* }}} */
