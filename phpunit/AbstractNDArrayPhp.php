<?php
namespace RindowTest\OpenBLAS;

use Interop\Polite\Math\Matrix\NDArray;
use InvalidArgumentException;
use ArrayObject;
use ArrayAccess;
use Traversable;
use SplFixedArray;
use Rindow\OpenBLAS\Buffer;
use RindowTest\OpenBLAS\HostBuffer;

abstract class AbstractNDArrayPhp implements NDArray
{
    protected $buffer;
    protected $size;
    protected $dtype;
    protected $offset;
    protected $shape;

    public function __construct(
        $array=null, int $dtype=null, array $shape=null, int $offset=null
    )
    {
        $dtype = $dtype ?? NDArray::float32;
        $offset = $offset ?? 0;
        if(is_array($array)||$array instanceof ArrayObject) {
            $dummyBuffer = new ArrayObject();
            $idx = 0;
            $this->array2Flat($array,$dummyBuffer,$idx,$prepare=true);
            $buffer = $this->newBuffer($idx,$dtype);
            $idx = 0;
            $this->array2Flat($array,$buffer,$idx,$prepare=false);
            if($shape===null) {
                $shape = $this->genShape($array);
            }
        } elseif(is_numeric($array)||is_bool($array)) {
            if(is_bool($array)&&$dtype!=NDArray::bool) {
                throw new InvalidArgumentException("unmatch dtype with bool value");
            }
            $buffer = $this->newBuffer(1,$dtype);
            $buffer[0] = $array;
            if($shape===null) {
                $shape = [];
            }
            $this->checkShape($shape);
            if(array_product($shape)!=1)
                throw new InvalidArgumentException("Invalid dimension size");
        } elseif($array===null && $shape!==null) {
            $this->checkShape($shape);
            $size = (int)array_product($shape);
            $buffer = $this->newBuffer($size,$dtype);
        } elseif($array===null && $shape!==null) {
            $this->checkShape($shape);
            $size = (int)array_product($shape);
            $buffer = $this->newBuffer($size,$dtype);
        } elseif($this->isBuffer($array)) {
            if($offset===null||!is_int($offset))
                throw new InvalidArgumentException("Must specify offset with the buffer");
            if($shape===null)
                throw new InvalidArgumentException("Invalid dimension size");
            $buffer = $array;
        } else {
            var_dump($array);var_dump($shape);
            throw new \Exception("Illegal array type");
        }
        $this->buffer = $buffer;
        $this->size = $buffer->count();
        $this->dtype = $buffer->dtype();
        $this->shape = $shape;
        $this->offset = $offset;
    }

    function newBuffer($size,$dtype)
    {
        return new HostBuffer($size,$dtype);
    }
    
    protected function isBuffer($buffer)
    {
        if($buffer instanceof SplFixedArray || $buffer instanceof Buffer) {
            return true;
        } else {
            return false;
        }
    }

    protected function array2Flat($A, $F, &$idx, $prepare)
    {
        if(is_array($A)) {
            ksort($A);
        } elseif($A instanceof ArrayObject) {
            $A->ksort();
        }

        $num = null;
        foreach ($A as $key => $value) {
            if(!is_int($key))
                throw new InvalidArgumentException("Dimension must be integer");
            if(is_array($value)||$value instanceof ArrayObject) {
                $num2 = $this->array2Flat($value, $F, $idx, $prepare);
                if($num===null) {
                    $num = $num2;
                } else {
                    if($num!=$num2)
                        throw new InvalidArgumentException("The shape of the dimension is broken");
                }
            } else {
                if($num!==null)
                    throw new InvalidArgumentException("The shape of the dimension is broken");
                if(!$prepare)
                    $F[$idx] = $value;
                $idx++;
            }
        }
        return count($A);
    }

    protected function flat2Array($F, &$idx, array $shape)
    {
        $size = array_shift($shape);
        if(count($shape)) {
            $A = [];
            for($i=0; $i<$size; $i++) {
                $A[$i] = $this->flat2Array($F,$idx,$shape);
            }
        }  else {
            $A = [];
            for($i=0; $i<$size; $i++) {
                $A[$i] = $F[$idx];
                $idx++;
            }
        }
        return $A;
    }
        
    protected function genShape($A)
    {
        $shape = [];
        while(is_array($A) || $A instanceof ArrayObject) {
            $shape[] = count($A);
            $A = $A[0];
        }
        return $shape;
    }

    protected function checkShape(array $shape)
    {
        foreach($shape as $num) {
            if(!is_int($num)) {
                throw new InvalidArgumentException(
                    "Invalid shape numbers. It gives ".gettype($num));
            }
            if($num<=0) {
                throw new InvalidArgumentException(
                    "Invalid shape numbers. It gives ".$num);
            }
        }
    }

    public function toArray()
    {
        if(count($this->shape)==0) {
            return $this->buffer[$this->offset];
        }
        $idx = $this->offset;
        return $this->flat2Array($this->buffer, $idx, $this->shape);
    }

    public function shape() : array { return $this->shape; }

    public function ndim() : int { return count($this->shape); }

    public function dtype() { return $this->dtype; }

    public function buffer() : ArrayAccess { return $this->buffer; }

    public function offset() : int { return $this->offset; }

    public function size() : int { return $this->buffer->count(); }

    public function reshape(array $shape) : NDArray
    {
        if(array_product($shape)==array_product($this->shape)) {
            $this->shape = $shape;
        } else {
            throw new \Exception("unmatch shape");
        }
        return $this;
    }
    public function offsetExists( $offset ) : bool { throw new \Excpetion('not implement'); }
    protected function doOffsetGet( $offset )
    {
        if(is_array($offset)) {
            throw new InvalidArgumentException("offset style is old renge style.");
        }
        // for single index specification
        if(is_numeric($offset)) {
            $shape = $this->shape;
            $max = array_shift($shape);
            if(count($shape)==0) {
                return $this->buffer[$this->offset+$offset];
            }
            $size = (int)array_product($shape);
            $new = new static($this->buffer,$this->dtype,$shape,$this->offset+$offset*$size);
            return $new;
        }

        // for range spesification
        $shape = $this->shape;
        array_shift($shape);
        if(is_array($offset)) {
            $start = $offset[0];
            $limit = $offset[1]+1;
        } else {
            $start = $offset->start();
            $limit = $offset->limit();
        }
        $rowsCount = $limit-$start;
        if(count($shape)>0) {
            $itemSize = (int)array_product($shape);
        } else {
            $itemSize = 1;
        }
        if($rowsCount<0) {
            throw new OutOfRangeException('Invalid range');
        }
        array_unshift($shape,$rowsCount);
        $size = (int)array_product($shape);
        $new = new static($this->buffer,$this->dtype,$shape,$this->offset+$start*$itemSize);
        return $new;
    }
    public function offsetSet( $offset , $value ) : void { throw new \Exception('not implement'); }
    public function offsetUnset( $offset ) : void { throw new \Exception('not implement'); }
    public function count() : int
    {
        return $this->buffer->count();
    }
    public function  getIterator() : Traversable  { throw new \Exception('not implement'); }
}
