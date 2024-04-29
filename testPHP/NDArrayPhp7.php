<?php
namespace RindowTest\OpenBLAS;

class NDArrayPhp extends AbstractNDArrayPhp
{
    public function offsetGet( $offset )
    {
        return $this->doOffsetGet( $offset );
    }
}