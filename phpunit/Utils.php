<?php
namespace RindowTest\OpenBLAS;

use Interop\Polite\Math\Matrix\NDArray;
use ArrayAccess;
use ArrayObject;
use Rindow\OpenBLAS\Buffer;

function R(
    int $start,
    int $limit,
) : Range
{
    if(func_num_args()!=2) {
        throw new InvalidArgumentException('R must have only two arguments: "start" and "limit".');
    }
    return new Range(start:$start,limit:$limit);
}

class Range
{
    protected mixed $start;
    protected mixed $limit;
    protected mixed $delta;

    public function __construct(
        int|float $limit,
        int|float $start=null,
        int|float $delta=null)
    {
        $this->limit = $limit;
        $this->start = $start ?? 0;
        $this->delta = $delta ?? (($limit>=$start)? 1 : -1);
    }

    public function start() : mixed
    {
        return $this->start;
    }

    public function limit() : mixed
    {
        return $this->limit;
    }

    public function delta() : mixed
    {
        return $this->delta;
    }
}


trait Utils
{
    public function zeros(array $shape,int $dtype=null)
    {
        $ndarray = $this->array(null,$dtype,$shape);
        $buffer = $ndarray->buffer();
        $size = $buffer->count();
        for($i=0;$i<$size;$i++) {
            $buffer[$i] = 0;
        }
        return $ndarray;
    }

    public function ones(array $shape,int $dtype=null)
    {
        $ndarray = $this->array(null,$dtype,$shape);
        $buffer = $ndarray->buffer();
        $size = $buffer->count();
        for($i=0;$i<$size;$i++) {
            $buffer[$i] = 1;
        }
        return $ndarray;
    }

    public function zerosLike(object $ndarray)
    {
        return $this->zeros($ndarray->shape(),$ndarray->dtype());
    }

    public function arange(int $count ,$start=null, $step=null, $dtype=null)
    {
        if($start===null)
            $start = 0;
        if($step===null)
            $step = 1;
        if($dtype===null) {
            if(is_int($start))
                $dtype = NDArray::int32;
            else
                $dtype = NDArray::float32;
        }
        $array = $this->zeros([$count], $dtype);
        $buffer = $array->buffer();
        $n = $start;
        for($i=0; $i<$count; $i++) {
            $buffer[$i] = $n;
            $n += $step;
        }
        return $array;
    }

    public function array(int|float|array $array=null, int $dtype=null, array $shape=null) : object
    {
        $ndarray = new class ($array, $dtype, $shape) implements NDArray {
            protected object $buffer;
            protected int $size;
            protected int $dtype;
            protected int $offset;
            protected array $shape;
            public function __construct(int|float|array|Buffer $array=null, int $dtype=null, array $shape=null, int $offset=null) {
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

            function newBuffer($size,$dtype) : object
            {
                return new Buffer($size,$dtype);
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
            public function offsetGet( $offset ) : mixed
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
                    $new = new self($this->buffer,$this->dtype,$shape,$this->offset+$offset*$size);
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
                $new = new self($this->buffer,$this->dtype,$shape,$this->offset+$start*$itemSize);
                return $new;
            }
            public function offsetSet( $offset , $value ) : void { throw new \Exception('not implement'); }
            public function offsetUnset( $offset ) : void { throw new \Exception('not implement'); }
            public function count() : int
            {
                return $this->buffer->count();
            }
            public function  getIterator() : Traversable  { throw new \Exception('not implement'); }
        };
        return $ndarray;
    }


    public function getMatlibVersion($matlib)
    {
        $config = $matlib->getConfig();
        if(strpos($config,'Matlib')===0) {
            $config = explode(' ',$config);
            return $config[1];
        } else {
            return '0.0.0';
        }
    }

    public function checkSkip($mark)
    {
        if(!in_array($mark,[
            //'multiply',
            //'duplicate'
            ])) {
            return false;
        }
        $this->markTestSkipped($mark);
        return true;
    }


    public function sum($n,$X,$offsetX,$incX)
    {
        $a = 0;
        for($i=0;$i<$n;$i++) {
            $a += $X[$offsetX+$i*$incX];
        }
        return $a;
    }

    protected function printableShapes($values)
    {
        if(!is_array($values)) {
            if($values instanceof NDArray)
                return '('.implode(',',$values->shape()).')';
            if(is_object($values))
                return '"'.get_class($values).'"';
            if(is_numeric($values) || is_string($values))
                return strval($values);
            return gettype($values);
        }
        $string = '[';
        foreach($values as $value) {
            if($string!='[') {
                $string .= ',';
            }
            $string .= $this->printableShapes($value);
        }
        $string .= ']';
        return $string;
    }

}