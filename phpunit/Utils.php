<?php
namespace RindowTest\OpenBLAS;
require_once __DIR__.'/AbstractNDArrayPhp.php';
require_once __DIR__.'/../testPHP/HostBuffer.php';
if(version_compare(phpversion(),'8.0.0','<')) {
    require_once __DIR__.'/../testPHP/NDArrayPhp7.php';
} else {
    require_once __DIR__.'/../testPHP/NDArrayPhp8.php';
}

use Interop\Polite\Math\Matrix\NDArray;
use Interop\Polite\Math\Matrix\BLAS;
use TypeError;
use InvalidArgumentException;
use RindowTest\OpenBLAS\NDArrayPhp;

function R(
    int $start,
    int $limit
) : Range
{
    if(func_num_args()!=2) {
        throw new InvalidArgumentException('R must have only two arguments: "start" and "limit".');
    }
    return new Range($limit,$start);
}

class Range
{
    protected $start;
    protected $limit;
    protected $delta;

    public function __construct(
        $limit,
        $start=null,
        $delta=null)
    {
        $this->limit = $limit;
        $this->start = $start ?? 0;
        $this->delta = $delta ?? (($limit>=$start)? 1 : -1);
    }

    public function start()
    {
        return $this->start;
    }

    public function limit()
    {
        return $this->limit;
    }

    public function delta()
    {
        return $this->delta;
    }
}

trait Utils
{
    public function alloc(array $shape,int $dtype=null)
    {
        $ndarray = $this->array(null,$dtype,$shape);
        return $ndarray;
    }

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

    public function array($array=null, int $dtype=null, array $shape=null) : object
    {
        $ndarray = new NDArrayPhp($array, $dtype, $shape);
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

    protected function copy(NDArray $x,NDArray $y=null) : NDArray
    {
        $blas = $this->getBlas();

        if($y==null) {
            $y = $this->zeros($x->shape(),$x->dtype());
        }
        $N = $x->size();
        $XX = $x->buffer();
        $offX = $x->offset();
        $YY = $y->buffer();
        $offY = $y->offset();
        $blas->copy($N,$XX,$offX,1,$YY,$offY,1);
        return $y;
    }

    protected function axpy(NDArray $x,NDArray $y=null,$alpha=null) : NDArray
    {
        $blas = $this->getBlas();

        if($y==null) {
            $y = $this->zeros($x->shape(),$x->dtype());
        }
        if($alpha===null) {
            $alpha = 1.0;
        }
        $N = $x->size();
        $XX = $x->buffer();
        $offX = $x->offset();
        $YY = $y->buffer();
        $offY = $y->offset();

        $blas->axpy($N,$alpha,$XX,$offX,1,$YY,$offY,1);
        return $y;
    }

    protected function iamax(NDArray $x) : int
    {
        $blas = $this->getBlas();

        $N = $x->size();
        $XX = $x->buffer();
        $offX = $x->offset();

        $y = $blas->iamax($N,$XX,$offX,1);
        return $y;
    }

    protected function scal(float $a,NDArray $x) : NDArray
    {
        $blas = $this->getBlas();

        $N = $x->size();
        $XX = $x->buffer();
        $offX = $x->offset();

        $blas->scal($N,$a,$XX,$offX,1);
        return $x;
    }

    protected function isComplex($dtype) : bool
    {
        return $dtype==NDArray::complex64||$dtype==NDArray::complex128;
    }

    protected function buildValByType($value, int $dtype)
    {
        if($this->isComplex($dtype)) {
            throw new InvalidArgumentException('complex value is not supported.');
        }
        return $value;
    }

    protected function transToCode(bool $trans,bool $conj) : int
    {
        if($trans) {
            return $conj ? BLAS::ConjTrans : BLAS::Trans;
        } else {
            return $conj ? BLAS::ConjNoTrans : BLAS::NoTrans;
        }
    }

    protected function complementTrans(?bool $trans,?bool $conj,int $dtype) : array
    {
        $trans = $trans ?? false;
        if($this->isComplex($dtype)) {
            $conj = $conj ?? $trans;
        } else {
            $conj = $conj ?? false;
        }
        return [$trans,$conj];
    }

    protected function abs($value) : float
    {
        if(is_numeric($value)) {
            return abs($value);
        }
        throw new InvalidArgumentException('complex value is not supported.');
        //$abs = sqrt(($value->real)**2+($value->imag)**2);
        //return $abs;
    }

    protected function isclose(NDArray $a, NDArray $b, $rtol=null, $atol=null) : bool
    {
        $blas = $this->getBlas();

        $isCpx = $this->isComplex($a->dtype());
        if($rtol===null) {
            //$rtol = $isCpx?C(1e-04):1e-04;
            $rtol = 1e-04;
        }
        if($atol===null) {
            $atol = 1e-07;
        }
        if($a->shape()!=$b->shape()) {
            return false;
        }
        // diff = b - a
        //$alpha =  $isCpx?C(-1):-1;
        $alpha = -1;
        $diffs = $this->copy($b);
        $this->axpy($a,$diffs,$alpha);
        $iDiffMax = $this->iamax($diffs);
        $diff = $this->abs($diffs->buffer()[$iDiffMax]);

        // close = atol + rtol * b
        $scalB = $this->copy($b);
        $this->scal($rtol,$scalB);
        $iCloseMax = $this->iamax($scalB);
        $close = $atol+$this->abs($scalB->buffer()[$iCloseMax]);

        return $diff < $close;
    }

    protected function argExpectException($class)
    {
        if(version_compare(phpversion(),'8.0.0','<')) {
            $this->expectException(TypeError::class);
        } else {
            $this->expectException($class);
        }
    }

    protected function argExpectExceptionMessage($message)
    {
        if(version_compare(phpversion(),'8.0.0','<')) {
            $this->expectExceptionMessage('Argument ');
        } else {
            $this->expectExceptionMessage($message);
        }
    }
}