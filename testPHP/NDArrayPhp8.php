<?php
namespace RindowTest\OpenBLAS;

use RindowTest\OpenBLAS\AbstractNDArrayPhp;

class NDArrayPhp extends AbstractNDArrayPhp
{
    public function offsetGet( $offset ) : mixed
    {
        return $this->doOffsetGet( $offset );
    }
}