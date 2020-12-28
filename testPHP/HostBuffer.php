<?php
namespace RindowTest\OpenBLAS;

use Rindow\OpenBLAS\Buffer as BufferImplement;
use Interop\Polite\Math\Matrix\LinearBuffer;

class HostBuffer extends BufferImplement implements LinearBuffer
{
}
