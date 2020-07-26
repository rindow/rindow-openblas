static inline int im3d_copyCell(
    zend_bool reverse,
    php_rindow_openblas_buffer_t *images,
    zend_long images_pos,
    zend_long filter_d,
    zend_long filter_h,
    zend_long filter_w,
    zend_long channels,
    zend_long channel_step,
    zend_long filter_d_step,
    zend_long filter_h_step,
    zend_long filter_w_step,
    zend_long vim_z,
    zend_long vim_y,
    zend_long vim_x,
    zend_long vim_d,
    zend_long vim_h,
    zend_long vim_w,
    php_rindow_openblas_buffer_t *out,
    zend_long out_pos,
    zend_long out_filter_step,
    zend_long out_channel_step
    )
{
    zend_long z;
    zend_long y;
    zend_long x;
    zend_long filter_d_pos;
    zend_long filter_h_pos;
    zend_long filter_w_pos;
    zend_long out_filter_pos;
    zend_long zz;
    zend_long yy;
    zend_long xx;
    zend_long channel_pos;
    zend_long out_channel_pos;
    zend_long c;

    filter_d_pos = images_pos;
    out_filter_pos = out_pos;
    for(z=0; z<filter_d; z++) {
        zz = z+vim_z;
        filter_h_pos = filter_d_pos;
        for(y=0; y<filter_h; y++) {
            yy = y+vim_y;
            filter_w_pos = filter_h_pos;
            for(x=0; x<filter_w; x++) {
                channel_pos = filter_w_pos;
                out_channel_pos = out_filter_pos;
                xx = x+vim_x;
                for(c=0; c<channels; c++) {
                    if(zz<0 || zz>=vim_d ||
                        yy<0 || yy>=vim_h ||
                        xx<0 || xx>=vim_w) {
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
            filter_h_pos += filter_h_step;
        }
        filter_d_pos += filter_d_step;
    }
    return 0;
}


static inline int im3d_stride(
    zend_long batches,
    zend_long batch_pos,
    zend_long batch_step,
    zend_long start_d,
    zend_long start_h,
    zend_long start_w,
    zend_long end_d,
    zend_long end_h,
    zend_long end_w,
    zend_long stride_d_step,
    zend_long stride_h_step,
    zend_long stride_w_step,
    zend_long start_vim_z,
    zend_long start_vim_y,
    zend_long start_vim_x,
    zend_long stride_d,
    zend_long stride_h,
    zend_long stride_w,
    zend_long reverse,
    php_rindow_openblas_buffer_t* images,
    zend_long filter_d,
    zend_long filter_h,
    zend_long filter_w,
    zend_long channels,
    zend_long channel_step,
    zend_long filter_d_step,
    zend_long filter_h_step,
    zend_long filter_w_step,
    zend_long vim_d,
    zend_long vim_h,
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
    zend_long y;
    zend_long z;
    zend_long stride_d_pos;
    zend_long stride_h_pos;
    zend_long stride_w_pos;
    zend_long vim_z;
    zend_long vim_y;
    zend_long vim_x;

    for(batch=0; batch<batches;batch++) {
        stride_d_pos = batch_pos+(start_d*stride_d_step);
        vim_z = start_vim_z;
        for(z=start_d;z<end_d;z++){
            stride_h_pos = stride_d_pos+(start_h*stride_h_step);
            vim_y = start_vim_y;
            for(y=start_h;y<end_h;y++){
                stride_w_pos = stride_h_pos+(start_w*stride_w_step);
                vim_x = start_vim_x;
                for(x=start_w;x<end_w;x++) {
                    int rc;
                    rc = im3d_copyCell(
                        reverse,
                        images,
                        stride_w_pos,
                        filter_d,
                        filter_h,
                        filter_w,
                        channels,
                        channel_step,
                        filter_d_step,
                        filter_h_step,
                        filter_w_step,
                        vim_z,
                        vim_y,
                        vim_x,
                        vim_d,
                        vim_h,
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
                stride_h_pos += stride_h_step;
                vim_y += stride_h;
            }
            stride_d_pos += stride_d_step;
                vim_z += stride_d;
        }
        batch_pos += batch_step;
    }
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
        bool $cols_channels_first,

        Buffer $cols,
        int $cols_offset,
        int $cols_size
    ) : void
 {{{ */
static PHP_METHOD(Math, im2col3d)
{
    zend_bool reverse;
    php_rindow_openblas_buffer_t* images;
    zend_long images_offset;
    zend_long images_size;
    zend_long batches;
    zend_long im_d;
    zend_long im_h;
    zend_long im_w;
    zend_long channels;
    zend_long filter_d;
    zend_long filter_h;
    zend_long filter_w;
    zend_long stride_d;
    zend_long stride_h;
    zend_long stride_w;
    zend_bool padding;
    zend_bool channels_first;
    zend_bool cols_channels_first;
    php_rindow_openblas_buffer_t* cols;
    zend_long cols_offset;
    zend_long cols_size;
    zval* images_obj=NULL;
    zval* cols_obj=NULL;
    zend_long out_d;
    zend_long out_h;
    zend_long out_w;
    zend_long start_d;
    zend_long start_h;
    zend_long start_w;
    zend_long end_d;
    zend_long end_h;
    zend_long end_w;
    zend_long stride_d_step;
    zend_long stride_h_step;
    zend_long stride_w_step;
    zend_long batch_step;
    zend_long channel_step;
    zend_long filter_d_step;
    zend_long filter_h_step;
    zend_long filter_w_step;
    zend_long out_filter_step;
    zend_long out_channel_step;
    zend_long out_cell_step;
    zend_long out_pos;
    zend_long batch_pos;
    zend_long start_vim_z;
    zend_long start_vim_y;
    zend_long start_vim_x;
    zend_long vim_d;
    zend_long vim_h;
    zend_long vim_w;

    ZEND_PARSE_PARAMETERS_START_EX(ZEND_PARSE_PARAMS_THROW, 21, 21)
        Z_PARAM_BOOL(reverse)
        Z_PARAM_OBJECT_OF_CLASS(images_obj,php_rindow_openblas_buffer_ce)
        Z_PARAM_LONG(images_offset)
        Z_PARAM_LONG(images_size)
        Z_PARAM_LONG(batches)

        Z_PARAM_LONG(im_d)
        Z_PARAM_LONG(im_h)
        Z_PARAM_LONG(im_w)
        Z_PARAM_LONG(channels)
        Z_PARAM_LONG(filter_d)
        
        Z_PARAM_LONG(filter_h)
        Z_PARAM_LONG(filter_w)
        Z_PARAM_LONG(stride_d)
        Z_PARAM_LONG(stride_h)
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

    if((batches*im_d*im_h*im_w*channels)
        !=images_size) {
        zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch images buffer size and images shape", 0);
        return;
    }
    out_d = ((im_d-filter_d)/stride_d)+1;
    out_h = ((im_h-filter_h)/stride_h)+1;
    out_w = ((im_w-filter_w)/stride_w)+1;
    if(padding) {
        if((batches*
            im_d*filter_d*
            im_h*filter_h*
            im_w*filter_w*
            channels)!=cols_size) {
            zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch cols buffer size and images shape", 0);
            return;
        }
        start_d = -((im_d-out_d)/2);
        start_h = -((im_h-out_h)/2);
        start_w = -((im_w-out_w)/2);
        end_d = start_d+im_d;
        end_h = start_h+im_h;
        end_w = start_w+im_w;
    } else {
        start_d = start_h = start_w = 0;
        end_d = out_d;
        end_h = out_h;
        end_w = out_w;
        if((batches*
            out_d*filter_d*
            out_h*filter_h*
            out_w*filter_w*
            channels)!=cols_size) {
            zend_throw_exception(spl_ce_InvalidArgumentException, "Unmatch cols buffer size and images shape", 0);
            return;
        }
    }
    if(channels_first) {
        // stride parameters
        stride_w_step = stride_w;
        stride_h_step = im_w*stride_h;
        stride_d_step = im_w*im_h*stride_d;
        batch_step = channels*im_w*im_h*im_d;
        // copy parameters
        channel_step = im_d*im_h*im_w;
        filter_w_step = 1;
        filter_h_step = im_w;
        filter_d_step = filter_h_step*im_h;
    } else {
        // stride parameters
        stride_w_step = channels*stride_w;
        stride_h_step = channels*im_w*stride_h;
        stride_d_step = channels*im_w*im_h*stride_d;
        batch_step = channels*im_w*im_h*im_d;
        // copy parameters
        channel_step = 1;
        filter_w_step = channels;
        filter_h_step = filter_w_step*im_w;
        filter_d_step = filter_h_step*im_h;
    }
    if(cols_channels_first) {
        out_filter_step = 1;
        out_channel_step = filter_d*filter_h*filter_w;
    } else {
        out_filter_step = channels;
        out_channel_step = 1;
    }
    out_cell_step = filter_d*filter_h*filter_w*channels;
    
    out_pos = cols_offset;
    batch_pos = images_offset;
    
    start_vim_z = start_d*stride_d;
    start_vim_y = start_h*stride_h;
    start_vim_x = start_w*stride_w;
    vim_d = (out_d-1)*stride_d+filter_d;
    vim_h = (out_h-1)*stride_h+filter_h;
    vim_w = (out_w-1)*stride_w+filter_w;
    
    im3d_stride(
        batches,
        batch_pos,
        batch_step,
        start_d,
        start_h,
        start_w,
        end_d,
        end_h,
        end_w,
        stride_d_step,
        stride_h_step,
        stride_w_step,
        start_vim_z,
        start_vim_y,
        start_vim_x,
        stride_d,
        stride_h,
        stride_w,
        
        reverse,
        images,
        filter_d,
        filter_h,
        filter_w,
        channels,
        channel_step,
        filter_d_step,
        filter_h_step,
        filter_w_step,
        vim_d,
        vim_h,
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