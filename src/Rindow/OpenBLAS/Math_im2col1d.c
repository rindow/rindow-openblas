static inline int im1d_copyCell(
    zend_bool reverse,
    php_rindow_openblas_buffer_t *images,
    zend_long images_pos,
    zend_long filter_w,
    zend_long channels,
    zend_long channel_step,
    zend_long filter_w_step,
    zend_long vim_x,
    zend_long vim_w,
    php_rindow_openblas_buffer_t *out,
    zend_long out_pos,
    zend_long out_filter_step,
    zend_long out_channel_step
    )
{
    zend_long x;
    zend_long filter_w_pos;
    zend_long out_filter_pos;
    zend_long xx;
    zend_long channel_pos;
    zend_long out_channel_pos;
    zend_long c;

    filter_w_pos = images_pos;
    out_filter_pos = out_pos;
    for(x=0; x<filter_w; x++) {
        channel_pos = filter_w_pos;
        out_channel_pos = out_filter_pos;
        xx = x+vim_x;
        for(c=0; c<channels; c++) {
            if(xx<0 || xx>=vim_w) {
                if(out_channel_pos<0 ||out_channel_pos>=out->size) {
                   zend_throw_exception(spl_ce_RuntimeException, "cols data out of range", 0);
                    return -1;
                }
                if(!reverse) {
                    if(images->dtype== php_rindow_openblas_dtype_float32) {
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
                    if(images->dtype== php_rindow_openblas_dtype_float32) {
                        ((float*)(out->data))[out_channel_pos]
                            = ((float*)(images->data))[channel_pos];
                    } else {
                        ((double*)(out->data))[out_channel_pos]
                            = ((double*)(images->data))[channel_pos];
                    }
                } else {
                    // Sum for Back propagation
                    if(images->dtype== php_rindow_openblas_dtype_float32) {
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


static inline int im1d_stride(
    zend_long batches,
    zend_long batch_pos,
    zend_long batch_step,
    zend_long start_w,
    zend_long end_w,
    zend_long stride_w_step,
    zend_long start_vim_x,
    zend_long stride_w,
    zend_bool reverse,
    php_rindow_openblas_buffer_t* images,
    zend_long filter_w,
    zend_long channels,
    zend_long channel_step,
    zend_long filter_w_step,
    zend_long vim_w,
    php_rindow_openblas_buffer_t* cols,
    zend_long out_pos,
    zend_long out_cell_step,
    zend_long out_filter_step,
    zend_long out_channel_step
    )
{
    zend_long batch;
    zend_long x;
    zend_long stride_w_pos;
    zend_long vim_x;

    for(batch=0; batch<batches;batch++) {
        stride_w_pos = batch_pos+(start_w*stride_w_step);
        vim_x = start_vim_x;
        for(x=start_w;x<end_w;x++) {
            int rc;
            rc = im1d_copyCell(
                reverse,
                images,
                stride_w_pos,
                filter_w,
                channels,
                channel_step,
                filter_w_step,
                vim_x,
                vim_w,
                cols,
                out_pos,
                out_filter_step,
                out_channel_step
            );
            if(rc) {
                return rc;
            }
            stride_w_pos += stride_w_step;
            vim_x += stride_w;
            out_pos += out_cell_step;
        }
        batch_pos += batch_step;
    }
    return 0;
}

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
        bool $cols_channels_first,

        Buffer $cols,
        int $cols_offset,
        int $cols_size
    ) : void
 {{{ */
static PHP_METHOD(Math, im2col1d)
{
    zend_bool reverse;
    php_rindow_openblas_buffer_t* images;
    zend_long images_offset;
    zend_long images_size;
    zend_long batches;
    zend_long im_w;
    zend_long channels;
    zend_long filter_w;
    zend_long stride_w;
    zend_bool padding;
    zend_bool channels_first;
    zend_bool cols_channels_first;
    php_rindow_openblas_buffer_t* cols;
    zend_long cols_offset;
    zend_long cols_size;
    zval* images_obj=NULL;
    zval* cols_obj=NULL;
    zend_long out_w;
    zend_long start_w;
    zend_long end_w;
    zend_long stride_w_step;
    zend_long batch_step;
    zend_long channel_step;
    zend_long filter_w_step;
    zend_long out_filter_step;
    zend_long out_channel_step;
    zend_long out_cell_step;
    zend_long out_pos;
    zend_long batch_pos;
    zend_long start_vim_x;
    zend_long vim_w;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 15, 15)
        Z_PARAM_BOOL(reverse)
        Z_PARAM_OBJECT_OF_CLASS(images_obj,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(images_offset)
        Z_PARAM_LONG(images_size)
        Z_PARAM_LONG(batches)

        Z_PARAM_LONG(im_w)
        Z_PARAM_LONG(channels)
        Z_PARAM_LONG(filter_w)
        Z_PARAM_LONG(stride_w)
        Z_PARAM_BOOL(padding)

        Z_PARAM_BOOL(channels_first)
        Z_PARAM_BOOL(cols_channels_first)
        Z_PARAM_OBJECT_OF_CLASS(cols_obj,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(cols_offset)
        Z_PARAM_LONG(cols_size)
    ZEND_PARSE_PARAMETERS_END();

    images = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(images_obj);
    if(php_rindow_openblas_assert_buffer_size(
        images, images_offset, images_size,
        "Invalid images buffer offset or size")) {
        return;
    }
    cols = Z_RINDOW_OPENBLAS_BUFFER_OBJ_P(cols_obj);
    if(php_rindow_openblas_assert_buffer_size(
        cols, cols_offset, cols_size,
        "Invalid cols buffer offset or size")) {
        return;
    }
    if(images->dtype!=php_rindow_openblas_dtype_float32 &&
        images->dtype!=php_rindow_openblas_dtype_float64) {
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

    if((batches*im_w*channels)
        !=images_size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch images buffer size and images shape", 0);
        return;
    }
    out_w = ((im_w-filter_w)/stride_w)+1;
    if(padding) {
        if((batches*
            im_w*filter_w*
            channels)!=cols_size) {
            zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch cols buffer size and images shape", 0);
            return;
        }
        start_w = -((im_w-out_w)/2);
        end_w = start_w+im_w;
    } else {
        start_w = 0;
        end_w = out_w;
        if((batches*
            out_w*filter_w*
            channels)!=cols_size) {
            zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch cols buffer size and images shape", 0);
            return;
        }
    }
    if(channels_first) {
        // stride parameters
        stride_w_step = stride_w;
        batch_step = channels*im_w;
        // copy parameters
        channel_step = im_w;
        filter_w_step = 1;
    } else {
        // stride parameters
        stride_w_step = channels*stride_w;
        batch_step = channels*im_w;
        // copy parameters
        channel_step = 1;
        filter_w_step = channels;
    }
    if(cols_channels_first) {
        out_filter_step = 1;
        out_channel_step = filter_w;
    } else {
        out_filter_step = channels;
        out_channel_step = 1;
    }
    out_cell_step = filter_w*channels;

    out_pos = cols_offset;
    batch_pos = images_offset;

    start_vim_x = start_w*stride_w;
    vim_w = (out_w-1)*stride_w+filter_w;

    im1d_stride(
        batches,
        batch_pos,
        batch_step,
        start_w,
        end_w,
        stride_w_step,
        start_vim_x,
        stride_w,

        reverse,
        images,
        filter_w,
        channels,
        channel_step,
        filter_w_step,
        vim_w,
        cols,
        out_pos,
        out_cell_step,
        out_filter_step,
        out_channel_step
    );
    return;
}
/* }}} */
