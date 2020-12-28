static inline int im2col2d_copyCell(
    zend_bool reverse,
    php_interop_polite_math_matrix_linear_buffer_t *images,
    zend_long images_pos,
    zend_long im_h,
    zend_long im_w,
    zend_long channels,
    zend_long channel_step,
    zend_long filter_h_step,
    zend_long filter_w_step,
    zend_long vim_y,
    zend_long vim_x,
    zend_long vfilter_h,
    zend_long vfilter_w,
    zend_long dilation_h,
    zend_long dilation_w,
    php_interop_polite_math_matrix_linear_buffer_t *out,
    zend_long out_pos,
    zend_long out_filter_step,
    zend_long out_channel_step
    )
{
    zend_long vfilter_y;
    zend_long vfilter_x;
    zend_long filter_h_pos;
    zend_long filter_w_pos;
    zend_long out_filter_pos;
    zend_long input_y;
    zend_long input_x;
    zend_long channel_pos;
    zend_long out_channel_pos;
    zend_long c;

    filter_h_pos = images_pos;
    out_filter_pos = out_pos;
    for(vfilter_y=0; vfilter_y<vfilter_h; vfilter_y+=dilation_h) {
        filter_w_pos = filter_h_pos;
        for(vfilter_x=0; vfilter_x<vfilter_w; vfilter_x+=dilation_w) {
            channel_pos = filter_w_pos;
            out_channel_pos = out_filter_pos;
            input_y = vim_y+vfilter_y;
            input_x = vim_x+vfilter_x;
            for(c=0; c<channels; c++) {
                if(input_y<0 || input_y>=im_h ||
                   input_x<0 || input_x>=im_w) {
                    if(out_channel_pos<0 ||out_channel_pos>=out->size) {
                       zend_throw_exception(spl_ce_RuntimeException, "cols data out of range", 0);
                        return -1;
                    }
                    if(!reverse) {
                        if(images->dtype == php_interop_polite_math_matrix_dtype_float32) {
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
        filter_h_pos += filter_h_step;
    }
    return 0;
}


static inline int im2col2d_execute(
    zend_bool reverse,
    php_interop_polite_math_matrix_linear_buffer_t* images,
    zend_long images_offset,
    zend_long images_size,
    zend_long batches,

    zend_long im_h,
    zend_long im_w,
    zend_long channels,
    zend_long filter_h,
    zend_long filter_w,

    zend_long stride_h,
    zend_long stride_w,
    zend_bool padding,
    zend_bool channels_first,
    zend_long dilation_h,

    zend_long dilation_w,
    zend_bool cols_channels_first,
    php_interop_polite_math_matrix_linear_buffer_t* cols,
    zend_long cols_offset,
    zend_long cols_size
    )
{
    zend_long out_h;
    zend_long out_w;
    zend_long stride_h_step;
    zend_long stride_w_step;
    zend_long batch_step;
    zend_long channel_step;
    zend_long filter_h_step;
    zend_long filter_w_step;
    zend_long out_filter_step;
    zend_long out_channel_step;
    zend_long out_cell_step;
    zend_long out_pos;
    zend_long batch_pos;
    zend_long padding_h;
    zend_long padding_w;
    zend_long im_w_step;
    zend_long im_h_step;

    zend_long batch;
    zend_long stride_h_pos;
    zend_long stride_w_pos;
    zend_long vim_y;
    zend_long vim_x;
    zend_long vim_h;
    zend_long vim_w;
    zend_long vfilter_h;
    zend_long vfilter_w;

    if((batches*im_h*im_w*channels)
        !=images_size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch images buffer size and images shape", 0);
        return -1;
    }
    out_h = ((im_h-(filter_h-1)*dilation_h-1)/stride_h)+1;
    out_w = ((im_w-(filter_w-1)*dilation_w-1)/stride_w)+1;
    if(out_h<=0 || out_w<=0) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Invalid shape or parameters.", 0);
        return -1;
    }
    if(padding) {
        if((batches*
            im_h*filter_h*
            im_w*filter_w*
            channels)!=cols_size) {
            zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch cols buffer size and images shape", 0);
            return -1;
        }
        padding_h = ((im_h-1)*stride_h-im_h+(filter_h-1)*dilation_h+1)/2;
        padding_w = ((im_w-1)*stride_w-im_w+(filter_w-1)*dilation_w+1)/2;
        out_h = im_h;
        out_w = im_w;
    } else {
        if((batches*
            out_h*filter_h*
            out_w*filter_w*
            channels)!=cols_size) {
            zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch cols buffer size and images shape", 0);
            return -1;
        }
        padding_h = 0;
        padding_w = 0;
    }
    if(channels_first) {
        im_w_step =    1;
        im_h_step =    im_w;
        channel_step = im_w*im_h;
        batch_step =   im_w*im_h*channels;
    } else {
        channel_step = 1;
        im_w_step =  channels;
        im_h_step =  channels*im_w;
        batch_step = channels*im_w*im_h;
    }
    stride_w_step = im_w_step*stride_w;
    stride_h_step = im_h_step*stride_h;
    filter_w_step = im_w_step*dilation_w;
    filter_h_step = im_h_step*dilation_h;

    if(cols_channels_first) {
        out_filter_step = 1;
        out_channel_step = filter_h*filter_w;
    } else {
        out_filter_step = channels;
        out_channel_step = 1;
    }
    out_cell_step = filter_h*filter_w*channels;

    batch_pos = images_offset-im_h_step*padding_h-im_w_step*padding_w;
    out_pos = cols_offset;

    vim_h = out_h*stride_h;
    vim_w = out_w*stride_w;
    vfilter_h = filter_h*dilation_h;
    vfilter_w = filter_w*dilation_w;

    for(batch=0; batch<batches;batch++) {
        stride_h_pos = batch_pos;
        for(vim_y=0;vim_y<vim_h;vim_y+=stride_h){
            stride_w_pos = stride_h_pos;
            for(vim_x=0;vim_x<vim_w;vim_x+=stride_w) {
                int rc;
                rc = im2col2d_copyCell(
                    reverse,
                    images,
                    stride_w_pos,
                    im_h,
                    im_w,
                    channels,
                    channel_step,
                    filter_h_step,
                    filter_w_step,
                    vim_y-padding_h,
                    vim_x-padding_w,
                    vfilter_h,
                    vfilter_w,
                    dilation_h,
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
            stride_h_pos += stride_h_step;
        }
        batch_pos += batch_step;
    }

    return 0;
}

/*
    cols := im2col2d(images)

    Method Rindow\OpenBLAS\Math::
    public function im2col2d(
        bool $reverse,
        Buffer $images,
        int $images_offset,
        int $images_size,
        int $batches,

        int $im_h,
        int $im_w,
        int $channels,
        int $filter_h,
        int $filter_w,

        int $stride_h,
        int $stride_w,
        bool $padding,
        bool $channels_first,
        int $dilation_h,

        int $dilation_w,
        bool $cols_channels_first,
        Buffer $cols,
        int $cols_offset,
        int $cols_size
    ) : void
 {{{ */
static PHP_METHOD(Math, im2col2d)
{
    zend_bool reverse;
    php_interop_polite_math_matrix_linear_buffer_t* images;
    zend_long images_offset;
    zend_long images_size;
    zend_long batches;
    zend_long im_h;
    zend_long im_w;
    zend_long channels;
    zend_long filter_h;
    zend_long filter_w;
    zend_long stride_h;
    zend_long stride_w;
    zend_bool padding;
    zend_bool channels_first;
    zend_long dilation_h;
    zend_long dilation_w;
    zend_bool cols_channels_first;
    php_interop_polite_math_matrix_linear_buffer_t* cols;
    zend_long cols_offset;
    zend_long cols_size;

    zval* images_obj=NULL;
    zval* cols_obj=NULL;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 20, 20)
        Z_PARAM_BOOL(reverse)
        Z_PARAM_ZVAL(images_obj) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(images_offset)
        Z_PARAM_LONG(images_size)
        Z_PARAM_LONG(batches)

        Z_PARAM_LONG(im_h)
        Z_PARAM_LONG(im_w)
        Z_PARAM_LONG(channels)
        Z_PARAM_LONG(filter_h)
        Z_PARAM_LONG(filter_w)

        Z_PARAM_LONG(stride_h)
        Z_PARAM_LONG(stride_w)
        Z_PARAM_BOOL(padding)
        Z_PARAM_BOOL(channels_first)
        Z_PARAM_LONG(dilation_h)

        Z_PARAM_LONG(dilation_w)
        Z_PARAM_BOOL(cols_channels_first)
        Z_PARAM_ZVAL(cols_obj) // Interop\Polite\Math\Matrix\LinearBuffer
        Z_PARAM_LONG(cols_offset)
        Z_PARAM_LONG(cols_size)
    ZEND_PARSE_PARAMETERS_END();

    images = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(images_obj);
    if(php_rindow_openblas_assert_buffer_size(
        images, images_offset, images_size,
        "Invalid images buffer offset or size")) {
        return;
    }
    cols = Z_INTEROP_POLITE_MATH_MATRIX_LINEAR_BUFFER_OBJ_P(cols_obj);
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

    im2col2d_execute(
        reverse,
        images,
        images_offset,
        images_size,
        batches,

        im_h,
        im_w,
        channels,
        filter_h,
        filter_w,

        stride_h,
        stride_w,
        padding,
        channels_first,
        dilation_h,

        dilation_w,
        cols_channels_first,
        cols,
        cols_offset,
        cols_size
    );
    return;
}
/* }}} */
