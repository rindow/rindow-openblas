<?php
namespace RindowTest\OpenBlas\MathTest;

use PHPUnit\Framework\TestCase;
use Interop\Polite\Math\Matrix\NDArray;
use Interop\Polite\Math\Matrix\BLAS;
use Rindow\Math\Matrix\MatrixOperator;
use Rindow\OpenBLAS\Math as Math;
use InvalidArgumentException;
use RuntimeException;
use TypeError;

/**
 * @requires extension rindow_openblas
 */
class Test extends TestCase
{
    public function getMath($mo)
    {
        $math = new Math();
        return $math;
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

    public function translate_amin(
        NDArray $X) : array
    {
        $N = $X->size();
        $XX = $X->buffer();
        $offX = $X->offset();
        return [$N,$XX,$offX,1];
    }

    public function translate_increment(
        NDArray $X,
        float $beta=null,
        float $alpha=null) : array
    {
        $N = $X->size();
        $XX = $X->buffer();
        $offX = $X->offset();
        if($alpha===null) {
            $alpha = 1.0;
        }
        if($beta===null) {
            $beta = 0.0;
        }
        return [$N,$alpha,$XX,$offX,1,$beta];
    }

    public function translate_maximum(
        NDArray $A,
        NDArray $X,
        ) : array
    {
        [$m,$n] = $A->shape();
        $AA = $A->buffer();
        $offA = $A->offset();
        $XX = $X->buffer();
        $offX = $X->offset();

        return [$m,$n,$AA,$offA,$n,$XX,$offX,1];
    }

    public function translate_multiply(
       NDArray $X,
       NDArray $A,
       bool $trans=null
       ) : array
    {
        if($trans===null)
            $trans = false;
        $shapeX = $X->shape();
        $shapeA = $A->shape();
        if($trans)
            $shapeA = array_reverse($shapeA);
        while(true) {
            $xd = array_pop($shapeX);
            if($xd===null)
                break;
            $ad = array_pop($shapeA);
            if($xd!==$ad)
                throw new InvalidArgumentException('Unmatch dimension size for broadcast.: '.
                    '['.implode(',',$X->shape()).'] ['.implode(',',$A->shape()).']');
        }
        $n = $X->size();
        $XX = $X->buffer();
        $offX = $X->offset();
        $m = $A->size()/$n;
        $AA = $A->buffer();
        $offA = $A->offset();
        if($trans) {
            [$m,$n] = [$n,$m];
        }

        return [
            $trans,
            $m,
            $n,
            $XX,$offX,1,
            $AA,$offA,$n,
        ];
    }

    public function translate_add(
       NDArray $X,
       NDArray $A,
       float $alpha=null,
       bool $trans=null
       ) : array
    {
        if($trans===null)
            $trans = false;
        if($alpha===null)
            $alpha = 1.0;
        $shapeX = $X->shape();
        $shapeA = $A->shape();
        if($trans)
            $shapeA = array_reverse($shapeA);
        while(true) {
            $xd = array_pop($shapeX);
            if($xd===null)
                break;
            $ad = array_pop($shapeA);
            if($xd!==$ad)
                throw new InvalidArgumentException('Unmatch dimension size for broadcast.: '.
                    '['.implode(',',$X->shape()).'] ['.implode(',',$A->shape()).']');
        }
        $n = $X->size();
        $XX = $X->buffer();
        $offX = $X->offset();
        $m = $A->size()/$n;
        $AA = $A->buffer();
        $offA = $A->offset();
        if($trans) {
            [$m,$n] = [$n,$m];
        }

        return [
            $trans,
            $m,
            $n,
            $alpha,
            $XX,$offX,1,
            $AA,$offA,$n,
        ];
    }

    public function translate_duplicate(NDArray $X, int $n=null, bool $trans=null,NDArray $A=null) : array
    {
        if($trans===null)
            $trans = false;
        if($A===null) {
            if(!$trans) {
                $A = $this->alloc(array_merge([$n],$X->shape()));
            } else {
                $A = $this->alloc(array_merge($X->shape(),[$n]));
            }
        } else {
            $shapeX = $X->shape();
            $shapeA = $A->shape();
            if($trans)
                $shapeA = array_reverse($shapeA);
            while(true) {
                $xd = array_pop($shapeX);
                if($xd===null)
                    break;
                $ad = array_pop($shapeA);
                if($xd!==$ad)
                    throw new InvalidArgumentException('Unmatch dimension size for broadcast.: '.
                        '['.implode(',',$X->shape()).'] => ['.implode(',',$A->shape()).']');
            }
        }

        $n = $X->size();
        $XX = $X->buffer();
        $offX = $X->offset();
        $m = $A->size()/$n;
        $AA = $A->buffer();
        $offA = $A->offset();
        if($trans) {
            [$m,$n] = [$n,$m];
        }

        return [
            $trans,
            $m,
            $n,
            $XX,$offX,1,
            $AA,$offA,$n
        ];
    }

    public function translate_square(
        NDArray $X
        ) : array
    {
        $n = $X->size();
        $XX = $X->buffer();
        $offX = $X->offset();

        return [
            $n,
            $XX,$offX,1
        ];
    }
/*
    public function translate_selectAxis0(
        NDArray $A,
        NDArray $X,
        NDArray $Y=null) : array
    {
        if($X->ndim()!=1) {
            throw new InvalidArgumentException('"X" must be 1D-NDArray.');
        }
        $countX = $X->shape()[0];
        if($A->ndim()==1) {
            $shape = $X->shape();
            $m = $A->shape()[0];
            $n = 1;
        } else {
            $shape = $A->shape();
            $m = $shape[0];
            $n = (int)($A->size()/$m);
            array_shift($shape);
            array_unshift($shape,$countX);
        }
        if($Y===null) {
            $Y = $this->alloc($shape,$A->dtype());
        } else {
            if($Y->shape()!=$shape) {
                throw new InvalidArgumentException('Unmatch size "Y" with "X" and "A" .');
            }
        }

        //if($A->ndim()==1) {
        //    $A = $A->reshape([$n,1]);
        //}
        //if($Y->ndim()==1) {
        //    $newY = $Y->reshape([$n,1]);
        //} else {
        //    $newY = $Y;
        //}
        //for($i=0;$i<$n;$i++) {
        //    $this->copy($A[$X[$i]],$newY[$i]);
        //}
        //return $Y;

        $AA = $A->buffer();
        $offA = $A->offset();
        $ldA = $n;
        $XX = $X->buffer();
        $offX = $X->offset();
        $YY = $Y->buffer();
        $offY = $Y->offset();
        $ldY = $n;

        return [
            $m,
            $n,
            $countX,
            $AA,$offA,$ldA,
            $XX,$offX,1,
            $YY,$offY,$ldY];
    }

    public function translate_selectAxis1(
        NDArray $A,
        NDArray $X,
        NDArray $Y=null) : array
    {
        if($A->ndim()!=2) {
            throw new InvalidArgumentException('"A" must be 2D-NDArray.');
        }
        if($X->ndim()!=1) {
            throw new InvalidArgumentException('"X" must be 1D-NDArray.');
        }
        [$m,$n] = $A->shape();
        if($X->size()!=$m) {
            throw new InvalidArgumentException('Unmatch size "X" with rows of "A".');
        }
        if($Y==null) {
            $Y = $this->alloc([$m],$A->dtype());
        }
        if($Y->ndim()!=1) {
            throw new InvalidArgumentException('"Y" must be 1D-NDArray.');
        }
        if($Y->size()!=$m) {
            throw new InvalidArgumentException('Unmatch size "Y" with rows of "A".');
        }

        $AA = $A->buffer();
        $offA = $A->offset();
        $ldA = $n;
        $XX = $X->buffer();
        $offX = $X->offset();
        $YY = $Y->buffer();
        $offY = $Y->offset();

        return [
            $m,
            $n,
            $AA,$offA,$ldA,
            $XX,$offX,1,
            $YY,$offY,1,
        ];
    }
*/
    public function translate_gather(
        bool $scatterAdd,
        NDArray $A,
        NDArray $X,
        int $axis=null,
        NDArray $B=null,
        $dtype=null) : array
    {
//echo "shapeX=[".implode(',',$X->shape())."],shapeA=[".implode(',',$A->shape())."]\n";
        if($axis===null) {
            $postfixShape = $A->shape();
            $prefixShape = $X->shape();
            $numClass = array_shift($postfixShape);
            $m = 1;
            $n = array_product($prefixShape);
            $k = array_product($postfixShape);
            $reductionDims = false;
            $outputShape = array_merge($prefixShape,$postfixShape);
        } else {
            $ndim = $A->ndim();
            $orgAxis = $axis;
            if($axis<0) {
                $axis = $ndim+$axis;
            }
            $postfixShape = $A->shape();
            $prefixShape = [];
            for($i=0;$i<$axis;$i++) {
                $prefixShape[] = array_shift($postfixShape);
            }
            $numClass = array_shift($postfixShape);
            $m = array_product($prefixShape);
            $n = array_product($postfixShape);
            $k = 1;
            $reductionDims = true;
            $outputShape = array_merge($prefixShape,$postfixShape);
            if($X->shape()!=$outputShape) {
                throw new InvalidArgumentException('Unmatch Shape:'.
                                        $this->printableShapes([$A,$X]));
            }
        }
//echo "outputShape=[".implode(',',$outputShape)."]\n";
        if($dtype===null) {
            $dtype = $A->dtype();
        }
        if($B==null) {
            $B = $this->alloc($outputShape,$dtype);
            $this->zeros($B);
        } else {
            if($B->shape()!=$outputShape) {
                throw new InvalidArgumentException("Unmatch output shape of dimension: ".
                                            $this->printableShapes([$outputShape,$B]));
            }
        }

        $AA = $A->buffer();
        $offA = $A->offset();
        $XX = $X->buffer();
        $offX = $X->offset();
        $BB = $B->buffer();
        $offB = $B->offset();

        if($scatterAdd) {
            $reverse=true;
            $addMode=true;
        } else {
            $reverse=false;
            $addMode=false;
        }
        if($reductionDims) {
            return [ $reduce=true,
                $reverse,
                $addMode,
                $m,
                $n,
                $numClass,
                $XX,$offX,
                $AA,$offA,
                $BB,$offB];
        } else {
            return [ $reduce=false,
                $reverse,
                $addMode,
                $n,
                $k,
                $numClass,
                $XX,$offX,
                $AA,$offA,
                $BB,$offB];
        }
    }

    public function translate_scatter(
        NDArray $X,
        NDArray $A,
        int $numClass,
        int $axis=null,
        NDArray $B=null,
        $dtype=null) : array
    {
//echo "shapeX=[".implode(',',$X->shape())."],shapeA=[".implode(',',$A->shape())."]\n";
//echo "axis=$axis,numClass=$numClass\n";
        if($axis===null) {
            $postfixShape = $A->shape();
            $prefixShape = $X->shape();
            //$numClass
            $ndimX = $X->ndim();
            $tmpShape = [];
            for($i=0;$i<$ndimX;$i++) {
                $tmpShape[] = array_shift($postfixShape);
            }
            if($tmpShape!=$prefixShape) {
                throw new InvalidArgumentException('Unmatch Shape:'.
                                        $this->printableShapes([$X,$A]));
            }
            $n = array_product($prefixShape);
            $k = array_product($postfixShape);
            $m = 1;
            $expandDims = false;
            $outputShape = array_merge([$numClass],$postfixShape);
        } else {
            $ndim = $A->ndim();
            $orgAxis = $axis;
            if($axis<0) {
                $axis = $ndim+$axis;
            }
            //if($axis<0 || $axis>$ndim-1) {
            //    throw new InvalidArgumentException("Invalid axis: ".$orgAxis);
            //}
            $postfixShape = $A->shape();
            $postfixX = $X->shape();
            if($postfixShape!=$postfixX) {
                throw new InvalidArgumentException('Unmatch Shape:'.
                                        $this->printableShapes([$X,$A]));
            }
            $prefixShape = [];
            for($i=0;$i<$axis;$i++) {
                $prefixShape[] = array_shift($postfixShape);
                array_shift($postfixX);
            }
            $m = array_product($prefixShape);
            $n = array_product($postfixShape);
            $k = 1;
            $expandDims = true;
            $outputShape = array_merge($prefixShape,[$numClass],$postfixShape);
        }
//echo "outputShape=[".implode(',',$outputShape)."]\n";
        if($dtype===null) {
            $dtype = $A->dtype();
        }
        if($B==null) {
            $B = $this->alloc($outputShape,$dtype);
            $this->zeros($B);
        } else {
            if($B->shape()!=$outputShape) {
                $shapeError = '('.implode(',',$A->shape()).'),('.implode(',',$B->shape()).')';
                throw new InvalidArgumentException("Unmatch shape of dimension: ".$shapeError);
            }
        }

        $AA = $A->buffer();
        $offA = $A->offset();
        $XX = $X->buffer();
        $offX = $X->offset();
        $BB = $B->buffer();
        $offB = $B->offset();

        if($expandDims) {
            return [ $reduce=true,
                $reverse=true,
                $addMode=false,
                $m,
                $n,
                $numClass,
                $XX,$offX,
                $BB,$offB,
                $AA,$offA];

        } else {
            return [ $reduce=false,
                $reverse=true,
                $addMode=false,
                $n,
                $k,
                $numClass,
                $XX,$offX,
                $BB,$offB,
                $AA,$offA];
        }
    }

    public function translate_onehot(
        NDArray $X,
        int $numClass,
        float $a=null,
        NDArray $Y=null) : array
    {
        if($X->ndim()!=1) {
            throw new InvalidArgumentException('"X" must be 1D-NDArray.');
        }
        $sizeX = $X->size();
        if($Y===null) {
            $Y = $this->zeros($this->alloc([$sizeX,$numClass]));
        }
        if($Y->ndim()!=2) {
            throw new InvalidArgumentException('"Y" must be 2D-NDArray.');
        }
        [$m,$n] = $Y->shape();
        if($m!=$sizeX || $n!=$numClass) {
            $shapeError = '('.implode(',',$X->shape()).'),('.implode(',',$Y->shape()).')';
            throw new InvalidArgumentException('Unmatch shape of dimension "X" and "Y" and "numClass": '.$shapeError);
        }
        if($a===null) {
            $a = 1.0;
        }
        $XX = $X->buffer();
        $offX = $X->offset();
        $YY = $Y->buffer();
        $offY = $Y->offset();
        $ldY = $n;

        return [
            $m,
            $n,
            $a,
            $XX,$offX,1,
            $YY,$offY,$ldY
        ];
    }

    public function translate_reduceSum(
        NDArray $A,
        int $axis=null,
        NDArray $B=null,
        $dtype=null) : array
    {
        $ndim = $A->ndim();
        if($axis<0) {
            $axis = $ndim+$axis;
        }
        if($axis<0 || $axis>$ndim-1) {
            throw new InvalidArgumentException("Invalid axis");
        }
        $postfixShape = $A->shape();
        $prefixShape = [];
        for($i=0;$i<$axis;$i++) {
            $prefixShape[] = array_shift($postfixShape);
        }
        $n = array_shift($postfixShape);
        $m = array_product($prefixShape);
        $k = array_product($postfixShape);
        $outputShape = array_merge($prefixShape,$postfixShape);
        if($dtype===null) {
            $dtype = $A->dtype();
        }
        if($B==null) {
            $B = $this->alloc($outputShape,$dtype);
        } else {
            if($B->shape()!=$outputShape) {
                $shapeError = '('.implode(',',$A->shape()).'),('.implode(',',$B->shape()).')';
                throw new InvalidArgumentException("Unmatch shape of dimension: ".$shapeError);
            }
        }
        $AA = $A->buffer();
        $offA = $A->offset();
        $BB = $B->buffer();
        $offB = $B->offset();
        return [
            $m,
            $n,
            $k,
            $AA,$offA,
            $BB,$offB
        ];
    }

   public function translate_astype(NDArray $X, $dtype, NDArray $Y) : array
   {
       $n = $X->size();
       $XX = $X->buffer();
       $offX = $X->offset();
       $YY = $Y->buffer();
       $offY = $Y->offset();

       return [
           $n,
           $dtype,
           $XX,$offX,1,
           $YY,$offY,1
       ];
   }

   public function translate_searchsorted(
       NDArray $A,
       NDArray $X,
       bool $right=null,
       $dtype=null,
       NDArray $Y=null
       ) : array
   {
       if($A->ndim()!=1) {
           throw new InvalidArgumentException('A must be 1D NDArray.');
       }
       if($right===null) {
           $right = false;
       }
       if($dtype===null) {
           $dtype = NDArray::uint32;
       }
       if($Y===null) {
           $Y = $this->alloc($X->shape(),$dtype);
       }
       $dtype = $Y->dtype();
       if($dtype!=NDArray::uint32&&$dtype!=NDArray::int32&&
           $dtype!=NDArray::uint64&&$dtype!=NDArray::int64) {
           throw new InvalidArgumentException('dtype of Y must be int32 or int64');
       }
       if($X->shape()!=$Y->shape()) {
           $shapeError = '('.implode(',',$X->shape()).'),('.implode(',',$Y->shape()).')';
           throw new InvalidArgumentException("Unmatch shape of dimension: ".$shapeError);
       }
       $m = $A->size();
       $AA = $A->buffer();
       $offA = $A->offset();
       $n = $X->size();
       $XX = $X->buffer();
       $offX = $X->offset();
       $YY = $Y->buffer();
       $offY = $Y->offset();

       return [
           $m,
           $AA,$offA,1,
           $n,
           $XX,$offX,1,
           $right,
           $YY,$offY,1
       ];
   }

   public function translate_cumsum(
       NDArray $X,
       bool $exclusive=null,
       bool $reverse=null,
       NDArray $Y=null
       ) : array
   {
       if($exclusive===null) {
           $exclusive = false;
       }
       if($reverse===null) {
           $reverse = false;
       }
       if($Y===null) {
           $Y = $this->alloc($X->shape(),$X->dtype());
       }
       if($X->shape()!=$Y->shape()) {
           $shapeError = '('.implode(',',$X->shape()).'),('.implode(',',$Y->shape()).')';
           throw new InvalidArgumentException("Unmatch shape of dimension: ".$shapeError);
       }
       $n = $X->size();
       $XX = $X->buffer();
       $offX = $X->offset();
       $YY = $Y->buffer();
       $offY = $Y->offset();

       return [
           $n,
           $XX,$offX,1,
           $exclusive,
           $reverse,
           $YY,$offY,1
       ];
   }

   public function testSumNormal()
   {
       $mo = new MatrixOperator();
       $math = $this->getMath($mo);

       $X = $mo->array([100,-10,-1000],NDArray::float32);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);
       $min = $math->sum($N,$XX,$offX,$incX);
       $this->assertEquals(-910,$min);

       $X = $mo->array([100,-10,-1000],NDArray::float64);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);
       $min = $math->sum($N,$XX,$offX,$incX);
       $this->assertEquals(-910,$min);

       $X = $mo->array([-100,-100,-120],NDArray::int8);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);
       $min = $math->sum($N,$XX,$offX,$incX);
       $this->assertEquals(-320,$min);

       $X = $mo->array([-1,-2,-3],NDArray::uint8);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);
       $min = $math->sum($N,$XX,$offX,$incX);
       $this->assertEquals(256*3-1-2-3,$min);

       $X = $mo->array([-100,-100,-120],NDArray::int16);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);
       $min = $math->sum($N,$XX,$offX,$incX);
       $this->assertEquals(-320,$min);

       $X = $mo->array([-1,-2,-3],NDArray::uint16);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);
       $min = $math->sum($N,$XX,$offX,$incX);
       $this->assertEquals(65536*3-1-2-3,$min);

       $X = $mo->array([-100,-100,-120],NDArray::int32);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);
       $min = $math->sum($N,$XX,$offX,$incX);
       $this->assertEquals(-320,$min);

       $X = $mo->array([-1,-2,-3],NDArray::uint32);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);
       $min = $math->sum($N,$XX,$offX,$incX);
       $this->assertEquals((2**32)*3-1-2-3,$min);

       $X = $mo->array([-100,-100,-120],NDArray::int64);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);
       $min = $math->sum($N,$XX,$offX,$incX);
       $this->assertEquals(-320,$min);

       $X = $mo->array([-1,-2,-3],NDArray::uint64);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);
       $min = $math->sum($N,$XX,$offX,$incX);
       $this->assertEquals((2**64)*3-1-2-3,$min);

       $X = $mo->array([true,false,true],NDArray::bool);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);
       $min = $math->sum($N,$XX,$offX,$incX);
       $this->assertEquals(2,$min);
   }

   public function testSumMinusN()
   {
       $mo = new MatrixOperator();
       $math = $this->getMath($mo);

       $X = $mo->array([100,-10,1]);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);

       $N = 0;
       $this->expectException(RuntimeException::class);
       $this->expectExceptionMessage('Argument n must be greater than 0.');
       $min = $math->sum($N,$XX,$offX,$incX);
   }

   public function testSumMinusOffsetX()
   {
       $mo = new MatrixOperator();
       $math = $this->getMath($mo);

       $X = $mo->array([100,-10,1]);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);

       $offX = -1;
       $this->expectException(RuntimeException::class);
       $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
       $min = $math->sum($N,$XX,$offX,$incX);
   }

   public function testSumMinusIncX()
   {
       $mo = new MatrixOperator();
       $math = $this->getMath($mo);

       $X = $mo->array([100,-10,1]);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);

       $incX = 0;
       $this->expectException(RuntimeException::class);
       $this->expectExceptionMessage('Argument incX must be greater than 0.');
       $min = $math->sum($N,$XX,$offX,$incX);
   }

   public function testSumIllegalBufferX()
   {
       $mo = new MatrixOperator();
       $math = $this->getMath($mo);

       $X = $mo->array([100,-10,1]);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);

       $XX = new \stdClass();
       $this->expectException(TypeError::class);
       $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
       $min = $math->sum($N,$XX,$offX,$incX);
   }

   public function testSumOverflowBufferXwithSize()
   {
       $mo = new MatrixOperator();
       $math = $this->getMath($mo);

       $X = $mo->array([100,-10,1]);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);

       $XX = $mo->array([100,-10])->buffer();
       $this->expectException(RuntimeException::class);
       $this->expectExceptionMessage('Vector specification too large for bufferX');
       $min = $math->sum($N,$XX,$offX,$incX);
   }

   public function testSumOverflowBufferXwithOffsetX()
   {
       $mo = new MatrixOperator();
       $math = $this->getMath($mo);

       $X = $mo->array([100,-10,1]);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);

       $offX = 1;
       $this->expectException(RuntimeException::class);
       $this->expectExceptionMessage('Vector specification too large for bufferX');
       $min = $math->sum($N,$XX,$offX,$incX);
   }

   public function testSumOverflowBufferXwithIncX()
   {
       $mo = new MatrixOperator();
       $math = $this->getMath($mo);

       $X = $mo->array([100,-10,1]);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);

       $incX = 2;
       $this->expectException(RuntimeException::class);
       $this->expectExceptionMessage('Vector specification too large for bufferX');
       $min = $math->sum($N,$XX,$offX,$incX);
   }

   public function testMinNormal()
   {
       $mo = new MatrixOperator();
       $math = $this->getMath($mo);

       $X = $mo->array([100,-10,1]);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);

       $min = $math->imin($N,$XX,$offX,$incX);
       $this->assertEquals(1,$min);
   }

   public function testMinMinusN()
   {
       $mo = new MatrixOperator();
       $math = $this->getMath($mo);

       $X = $mo->array([100,-10,1]);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);

       $N = 0;
       $this->expectException(RuntimeException::class);
       $this->expectExceptionMessage('Argument n must be greater than 0.');
       $min = $math->imin($N,$XX,$offX,$incX);
   }

   public function testMinMinusOffsetX()
   {
       $mo = new MatrixOperator();
       $math = $this->getMath($mo);

       $X = $mo->array([100,-10,1]);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);

       $offX = -1;
       $this->expectException(RuntimeException::class);
       $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
       $min = $math->imin($N,$XX,$offX,$incX);
   }

   public function testMinMinusIncX()
   {
       $mo = new MatrixOperator();
       $math = $this->getMath($mo);

       $X = $mo->array([100,-10,1]);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);

       $incX = 0;
       $this->expectException(RuntimeException::class);
       $this->expectExceptionMessage('Argument incX must be greater than 0.');
       $min = $math->imin($N,$XX,$offX,$incX);
   }

   public function testMinIllegalBufferX()
   {
       $mo = new MatrixOperator();
       $math = $this->getMath($mo);

       $X = $mo->array([100,-10,1]);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);

       $XX = new \stdClass();
       $this->expectException(TypeError::class);
       $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
       $min = $math->imin($N,$XX,$offX,$incX);
   }

   public function testMinOverflowBufferXwithSize()
   {
       $mo = new MatrixOperator();
       $math = $this->getMath($mo);

       $X = $mo->array([100,-10,1]);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);

       $XX = $mo->array([100,-10])->buffer();
       $this->expectException(RuntimeException::class);
       $this->expectExceptionMessage('Vector specification too large for bufferX');
       $min = $math->imin($N,$XX,$offX,$incX);
   }

   public function testMinOverflowBufferXwithOffsetX()
   {
       $mo = new MatrixOperator();
       $math = $this->getMath($mo);

       $X = $mo->array([100,-10,1]);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);

       $offX = 1;
       $this->expectException(RuntimeException::class);
       $this->expectExceptionMessage('Vector specification too large for bufferX');
       $min = $math->imin($N,$XX,$offX,$incX);
   }

   public function testMinOverflowBufferXwithIncX()
   {
       $mo = new MatrixOperator();
       $math = $this->getMath($mo);

       $X = $mo->array([100,-10,1]);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);

       $incX = 2;
       $this->expectException(RuntimeException::class);
       $this->expectExceptionMessage('Vector specification too large for bufferX');
       $min = $math->imin($N,$XX,$offX,$incX);
   }

   public function testMaxNormal()
   {
       $mo = new MatrixOperator();
       $math = $this->getMath($mo);

       $X = $mo->array([100,-10,-1000]);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);

       $min = $math->imax($N,$XX,$offX,$incX);
       $this->assertEquals(0,$min);
   }

   public function testMaxMinusN()
   {
       $mo = new MatrixOperator();
       $math = $this->getMath($mo);

       $X = $mo->array([100,-10,1]);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);

       $N = 0;
       $this->expectException(RuntimeException::class);
       $this->expectExceptionMessage('Argument n must be greater than 0.');
       $min = $math->imax($N,$XX,$offX,$incX);
   }

   public function testMaxMinusOffsetX()
   {
       $mo = new MatrixOperator();
       $math = $this->getMath($mo);

       $X = $mo->array([100,-10,1]);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);

       $offX = -1;
       $this->expectException(RuntimeException::class);
       $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
       $min = $math->imax($N,$XX,$offX,$incX);
   }

   public function testMaxMinusIncX()
   {
       $mo = new MatrixOperator();
       $math = $this->getMath($mo);

       $X = $mo->array([100,-10,1]);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);

       $incX = 0;
       $this->expectException(RuntimeException::class);
       $this->expectExceptionMessage('Argument incX must be greater than 0.');
       $min = $math->imax($N,$XX,$offX,$incX);
   }

   public function testMaxIllegalBufferX()
   {
       $mo = new MatrixOperator();
       $math = $this->getMath($mo);

       $X = $mo->array([100,-10,1]);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);

       $XX = new \stdClass();
       $this->expectException(TypeError::class);
       $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
       $min = $math->imax($N,$XX,$offX,$incX);
   }

   public function testMaxOverflowBufferXwithSize()
   {
       $mo = new MatrixOperator();
       $math = $this->getMath($mo);

       $X = $mo->array([100,-10,1]);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);

       $XX = $mo->array([100,-10])->buffer();
       $this->expectException(RuntimeException::class);
       $this->expectExceptionMessage('Vector specification too large for bufferX');
       $min = $math->imax($N,$XX,$offX,$incX);
   }

   public function testMaxOverflowBufferXwithOffsetX()
   {
       $mo = new MatrixOperator();
       $math = $this->getMath($mo);

       $X = $mo->array([100,-10,1]);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);

       $offX = 1;
       $this->expectException(RuntimeException::class);
       $this->expectExceptionMessage('Vector specification too large for bufferX');
       $min = $math->imax($N,$XX,$offX,$incX);
   }

   public function testMaxOverflowBufferXwithIncX()
   {
       $mo = new MatrixOperator();
       $math = $this->getMath($mo);

       $X = $mo->array([100,-10,1]);
       [$N,$XX,$offX,$incX] =
           $this->translate_amin($X);

       $incX = 2;
       $this->expectException(RuntimeException::class);
       $this->expectExceptionMessage('Vector specification too large for bufferX');
       $min = $math->imax($N,$XX,$offX,$incX);
   }


    public function testIncrementNormal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,10,2);

        $math->increment($N,$alpha,$XX,$offX,$incX,$beta);
        $this->assertEquals([12,14,16],$X->toArray());
    }

    public function testIncrementInvalidArgments()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,10,2);

        $this->expectException(TypeError::class);
        if(version_compare(PHP_VERSION, '8.0.0')<0) {
            $this->expectExceptionMessage('parameter 2 to be float');
        } else {
            $this->expectExceptionMessage('Argument #2 ($alpha) must be of type float');
        }
        $math->increment($N,new \stdClass(),$XX,$offX,$incX,$beta);
    }

    public function testIncrementMinusN()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,10,2);

        $N = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->increment($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testIncrementMinusOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,10,2);

        $offX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->increment($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testIncrementMinusIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,10,2);

        $incX = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->increment($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testIncrementIllegalBufferX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,10,2);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->increment($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testIncrementOverflowBufferXwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,10,2);

        $XX = $mo->array([1,2])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->increment($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testIncrementOverflowBufferXwithOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,10,2);

        $offX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->increment($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testIncrementOverflowBufferXwithIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,10,2);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->increment($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testReciprocalNormal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $math->reciprocal($N,$alpha,$XX,$offX,$incX,$beta);
        // X := 1 / ( alpha * X + beta )
        //    = [1/1, 1/2, 1/4]
        $this->assertEquals([1,0.5,0.25],$X->toArray());
    }

    public function testReciprocalZeroDivide()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([4,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,$beta=0,$alpha=1);

        // X := 1 / ( alpha * X + beta )

        // *** CAUTION ***
        // disable checking for INFINITY values
        //$this->expectException(RuntimeException::class);
        //$this->expectExceptionMessage('Zero divide.');

        $math->reciprocal($N,$alpha,$XX,$offX,$incX,$beta);
        $this->assertEquals(
            [0.25,0.5,INF],
            $X->toArray());
    }

    public function testReciprocalMinusN()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $N = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->reciprocal($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testReciprocalMinusOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $offX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->reciprocal($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testReciprocalMinusIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $incX = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->reciprocal($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testReciprocalIllegalBufferX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->reciprocal($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testReciprocalOverflowBufferXwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $XX = $mo->array([3,2])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->reciprocal($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testReciprocalOverflowBufferXwithOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $offX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->reciprocal($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testReciprocalOverflowBufferXwithIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->reciprocal($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testMaximumNormal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
        $this->assertEquals([[2,3],[2,3],[3,4]],$A->toArray());
    }

    public function testMaximumMinusM()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $M = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumMinusN()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $N = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumMinusOffsetA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offA = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetA must be greater than equals 0.');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumMinusLdA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $ldA = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument ldA must be greater than 0.');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumIllegalBufferA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $AA = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumOverflowBufferAwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $AA = $mo->array([1,2,2,3,3])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumOverflowBufferXwithOffsetA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offA = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumOverflowBufferXwithLdA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $ldA = 3;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumMinusOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumMinusIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $incX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumIllegalBufferX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumOverflowBufferXwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $XX = $mo->array([1])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumOverflowBufferXwithOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumOverflowBufferXwithIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumNormal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
        $this->assertEquals([[1,2],[2,3],[2,3]],$A->toArray());
    }

    public function testMinimumMinusM()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $M = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumMinusN()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $N = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumMinusOffsetA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offA = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetA must be greater than equals 0.');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumMinusLdA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $ldA = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument ldA must be greater than 0.');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumIllegalBufferA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $AA = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumOverflowBufferAwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $AA = $mo->array([1,2,2,3,3])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumOverflowBufferXwithOffsetA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offA = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumOverflowBufferXwithLdA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $ldA = 3;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumMinusOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumMinusIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $incX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumIllegalBufferX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumOverflowBufferXwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $XX = $mo->array([1])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumOverflowBufferXwithOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumOverflowBufferXwithIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }
    
    public function testGreaterNormal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
        $this->assertEquals([[0,0],[0,0],[1,1]],$A->toArray());
    }

    public function testGreaterMinusM()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $M = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterMinusN()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $N = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterMinusOffsetA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offA = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetA must be greater than equals 0.');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterMinusLdA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $ldA = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument ldA must be greater than 0.');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterIllegalBufferA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $AA = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterOverflowBufferAwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $AA = $mo->array([1,2,2,3,3])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterOverflowBufferXwithOffsetA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offA = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterOverflowBufferXwithLdA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $ldA = 3;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterMinusOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterMinusIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $incX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterIllegalBufferX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterOverflowBufferXwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $XX = $mo->array([1])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterOverflowBufferXwithOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterOverflowBufferXwithIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterEqualNormal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $math->greaterEqual($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
        $this->assertEquals([[0,0],[1,1],[1,1]],$A->toArray());
    }

    public function testLessNormal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
        $this->assertEquals([[1,1],[0,0],[0,0]],$A->toArray());
    }

    public function testLessMinusM()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $M = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessMinusN()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $N = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessMinusOffsetA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offA = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetA must be greater than equals 0.');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessMinusLdA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $ldA = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument ldA must be greater than 0.');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessIllegalBufferA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $AA = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessOverflowBufferAwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $AA = $mo->array([1,2,2,3,3])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessOverflowBufferXwithOffsetA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offA = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessOverflowBufferXwithLdA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $ldA = 3;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessMinusOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessMinusIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $incX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessIllegalBufferX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessOverflowBufferXwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $XX = $mo->array([1])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessOverflowBufferXwithOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessOverflowBufferXwithIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessEqualNormal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2],[2,3],[3,4]]);
        $X = $mo->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $math->lessEqual($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
        $this->assertEquals([[1,1],[1,1],[0,0]],$A->toArray());
    }

    public function testMultiplySameSizeNormal()
    {
        if($this->checkSkip('multiply')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
        $this->assertEquals([10,200,3000],$A->toArray());
    }

    public function testMultiplyBroadcastNormal()
    {
        if($this->checkSkip('multiply')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([[10,100,1000],[-1,-1,-1]]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
        $this->assertEquals([[10,200,3000],[-1,-2,-3]],$A->toArray());
    }

    public function testMultiplyBroadcastTranspose()
    {
        if($this->checkSkip('multiply')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([[10,100],[1000,10000],[-1,-1]]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A,true);

        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
        $this->assertEquals([[10,100],[2000,20000],[-3,-3]],$A->toArray());
    }

    public function testMultiplyMinusM()
    {
        if($this->checkSkip('multiply')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $M = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testMultiplyMinusN()
    {
        if($this->checkSkip('multiply')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $N = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testMultiplyMinusOffsetX()
    {
        if($this->checkSkip('multiply')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $offX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testMultiplyMinusIncX()
    {
        if($this->checkSkip('multiply')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $incX = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testMultiplyIllegalBufferX()
    {
        if($this->checkSkip('multiply')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testMultiplyOverflowBufferXwithSize()
    {
        if($this->checkSkip('multiply')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $XX = $mo->array([1,2])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testMultiplyOverflowBufferXwithOffsetX()
    {
        if($this->checkSkip('multiply')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $offX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testMultiplyOverflowBufferXwithIncX()
    {
        if($this->checkSkip('multiply')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testMultiplyMinusOffsetA()
    {
        if($this->checkSkip('multiply')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $offA = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetA must be greater than equals 0.');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testMultiplyMinusIncA()
    {
        if($this->checkSkip('multiply')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $ldA = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument ldA must be greater than 0.');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testMultiplyIllegalBufferA()
    {
        if($this->checkSkip('multiply')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $AA = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testMultiplyOverflowBufferAwithSize()
    {
        if($this->checkSkip('multiply')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $AA = $mo->array([1,2])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testMultiplyOverflowBufferXwithOffsetA()
    {
        if($this->checkSkip('multiply')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $offA = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testMultiplyOverflowBufferXwithLdA()
    {
        if($this->checkSkip('multiply')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([[10,100,1000],[10,100,1000]]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $ldA = 4;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testaddSameSizeNormal()
    {
        if($this->checkSkip('add')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A,-1);

        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
        $this->assertEquals([9,98,997],$A->toArray());
    }

    public function testaddBroadcastNormal()
    {
        if($this->checkSkip('add')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([[10,100,1000],[-1,-1,-1]]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
        $this->assertEquals([[11,102,1003],[0,1,2]],$A->toArray());
    }

    public function testaddBroadcastTranspose()
    {
        if($this->checkSkip('add')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([[10,100],[1000,10000],[-1,-1]]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A,null,true);

        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
        $this->assertEquals([[11,101],[1002,10002],[2,2]],$A->toArray());
    }

    public function testaddMinusM()
    {
        if($this->checkSkip('add')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $M = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testaddMinusN()
    {
        if($this->checkSkip('add')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $N = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testaddMinusOffsetX()
    {
        if($this->checkSkip('add')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $offX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testaddMinusIncX()
    {
        if($this->checkSkip('add')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $incX = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testaddIllegalBufferX()
    {
        if($this->checkSkip('add')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testaddOverflowBufferXwithSize()
    {
        if($this->checkSkip('add')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $XX = $mo->array([1,2])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testaddOverflowBufferXwithOffsetX()
    {
        if($this->checkSkip('add')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $offX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testaddOverflowBufferXwithIncX()
    {
        if($this->checkSkip('add')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testaddMinusOffsetA()
    {
        if($this->checkSkip('add')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $offA = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetA must be greater than equals 0.');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testaddMinusIncA()
    {
        if($this->checkSkip('add')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $ldA = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument ldA must be greater than 0.');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testaddIllegalBufferA()
    {
        if($this->checkSkip('add')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $AA = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testaddOverflowBufferAwithSize()
    {
        if($this->checkSkip('add')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $AA = $mo->array([1,2])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testaddOverflowBufferXwithOffsetA()
    {
        if($this->checkSkip('add')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $offA = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testaddOverflowBufferXwithLdA()
    {
        if($this->checkSkip('add')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([[10,100,1000],[10,100,1000]]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $ldA = 4;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testDuplicateSameSizeNormal()
    {
        if($this->checkSkip('duplicate')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
        $this->assertEquals([1,2,3],$A->toArray());
    }

    public function testDuplicateBroadcastNormal()
    {
        if($this->checkSkip('duplicate')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([[10,100,1000],[-1,-1,-1]]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
        $this->assertEquals([[1,2,3],[1,2,3]],$A->toArray());
    }

    public function testDuplicateBroadcastTranspose()
    {
        if($this->checkSkip('duplicate')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2]);
        $A = $mo->array([[10,100,1000],[-1,-1,-1]]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,true,$A);

        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
        $this->assertEquals([[1,1,1],[2,2,2]],$A->toArray());
    }

    public function testDuplicateMinusM()
    {
        if($this->checkSkip('duplicate')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $M = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testDuplicateMinusN()
    {
        if($this->checkSkip('duplicate')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $N = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testDuplicateMinusOffsetX()
    {
        if($this->checkSkip('duplicate')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $offX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testDuplicateMinusIncX()
    {
        if($this->checkSkip('duplicate')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $incX = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testDuplicateIllegalBufferX()
    {
        if($this->checkSkip('duplicate')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testDuplicateOverflowBufferXwithSize()
    {
        if($this->checkSkip('duplicate')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $XX = $mo->array([2])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testDuplicateOverflowBufferXwithOffsetX()
    {
        if($this->checkSkip('duplicate')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $offX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testDuplicateOverflowBufferXwithIncX()
    {
        if($this->checkSkip('duplicate')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testDuplicateMinusOffsetA()
    {
        if($this->checkSkip('duplicate')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $offA = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetA must be greater than equals 0.');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testDuplicateMinusLdA()
    {
        if($this->checkSkip('duplicate')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $ldA = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument ldA must be greater than 0.');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testDuplicateIllegalBufferA()
    {
        if($this->checkSkip('duplicate')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $AA = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testDuplicateOverflowBufferAwithSize()
    {
        if($this->checkSkip('duplicate')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $AA = $mo->array([1,2])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testDuplicateOverflowBufferAwithOffsetA()
    {
        if($this->checkSkip('duplicate')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $offA = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testDuplicateOverflowBufferAwithLdA()
    {
        if($this->checkSkip('duplicate')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $A = $mo->array([[10,100,1000],[10,100,1000]]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $ldA = 4;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testSquareNormal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $math->square($N,$XX,$offX,$incX);
        $this->assertEquals([1,4,9],$X->toArray());
    }

    public function testsquareMinusN()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $N = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->square($N,$XX,$offX,$incX);
    }

    public function testsquareMinusOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $offX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->square($N,$XX,$offX,$incX);
    }

    public function testsquareMinusIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $incX = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->square($N,$XX,$offX,$incX);
    }

    public function testsquareIllegalBufferX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->square($N,$XX,$offX,$incX);
    }

    public function testsquareOverflowBufferXwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $XX = $mo->array([1,2])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->square($N,$XX,$offX,$incX);
    }

    public function testsquareOverflowBufferXwithOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $offX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->square($N,$XX,$offX,$incX);
    }

    public function testsquareOverflowBufferXwithIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->square($N,$XX,$offX,$incX);
    }

    public function testsqrtNormal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,1,4,9]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $math->sqrt($N,$XX,$offX,$incX);
        $this->assertEquals([0,1,2,3],$X->toArray());
    }

    public function testsqrtIllegalValue()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,4,-1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        // *** CAUTION ***
        // disable checking for INFINITY values
        //$this->expectException(RuntimeException::class);
        //$this->expectExceptionMessage('Invalid value in sqrt.');

        $math->sqrt($N,$XX,$offX,$incX);
        $this->assertEquals(1.0, $X[0]);
        $this->assertEquals(2.0, $X[1]);
        $this->assertTrue(is_nan($X[2])); // -NAN
    }

    public function testsqrtMinusN()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $N = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->sqrt($N,$XX,$offX,$incX);
    }

    public function testsqrtMinusOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $offX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->sqrt($N,$XX,$offX,$incX);
    }

    public function testsqrtMinusIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $incX = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->sqrt($N,$XX,$offX,$incX);
    }

    public function testssqrtIllegalBufferX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->sqrt($N,$XX,$offX,$incX);
    }

    public function testsqrtOverflowBufferXwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $XX = $mo->array([1,2])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->sqrt($N,$XX,$offX,$incX);
    }

    public function testsqrtOverflowBufferXwithOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $offX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->sqrt($N,$XX,$offX,$incX);
    }

    public function testsqrtOverflowBufferXwithIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->sqrt($N,$XX,$offX,$incX);
    }

    public function testrsqrtNormal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,4,16]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,1,2);

        $math->rsqrt($N,$alpha,$XX,$offX,$incX,$beta);
        // X := 1 / (a * sqrt(X) + b)

        $math->reciprocal($N,1.0,$XX,$offX,$incX,0.0);
        $math->increment($N,1.0,$XX,$offX,$incX,-1.0);
        $math->increment($N,0.5,$XX,$offX,$incX,0);
        $math->square($N,$XX,$offX,$incX);

        $this->assertEquals([1,4,16],$X->toArray());
    }

    public function testrsqrtZeroDivide()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([4,1,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,$beta=0,$alpha=1);

        // X := 1 / (a * sqrt(X) + b)

        // *** CAUTION ***
        // disable checking for INFINITY values
        //$this->expectException(RuntimeException::class);
        //$this->expectExceptionMessage('Zero divide.');

        $math->rsqrt($N,$alpha,$XX,$offX,$incX,$beta);
        $this->assertEquals(
            [0.5, 1.0, INF],
            $X->toArray());
    }

    public function testrsqrtInvalidSqrt()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([4,1,-1]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,0,1);

        // X := 1 / (a * sqrt(X) + b)
        // *** CAUTION ***
        // disable checking for INFINITY values
        //$this->expectException(RuntimeException::class);
        //$this->expectExceptionMessage('Invalid value in sqrt.');

        $math->rsqrt($N,$alpha,$XX,$offX,$incX,$beta);
        $this->assertEquals(0.5, $X[0]);
        $this->assertEquals(1.0, $X[1]);
        $this->assertTrue(is_nan($X[2]));
    }

    public function testrsqrtMinusN()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $N = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->rsqrt($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testrsqrtMinusOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $offX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->rsqrt($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testrsqrtMinusIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $incX = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->rsqrt($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testrsqrtIllegalBufferX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->rsqrt($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testrsqrtOverflowBufferXwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $XX = $mo->array([3,2])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->rsqrt($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testrsqrtOverflowBufferXwithOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $offX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->rsqrt($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testrsqrtOverflowBufferXwithIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->rsqrt($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testpowNormal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX,$alpha] =
            $this->translate_maximum($X,2);

        $math->pow($N,$XX,$offX,$incX,$alpha);
        $this->assertEquals([1,4,9],$X->toArray());
    }

    public function testpowMinusN()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX,$alpha] =
            $this->translate_maximum($X,2);

        $N = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->pow($N,$XX,$offX,$incX,$alpha);
    }

    public function testpowMinusOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX,$alpha] =
            $this->translate_maximum($X,2);

        $offX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->pow($N,$XX,$offX,$incX,$alpha);
    }

    public function testpowMinusIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX,$alpha] =
            $this->translate_maximum($X,2);

        $incX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->pow($N,$XX,$offX,$incX,$alpha);
    }

    public function testpowIllegalBufferX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX,$alpha] =
            $this->translate_maximum($X,2);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->pow($N,$XX,$offX,$incX,$alpha);
    }

    public function testpowOverflowBufferXwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX,$alpha] =
            $this->translate_maximum($X,2);

        $XX = $mo->array([1,2])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->pow($N,$XX,$offX,$incX,$alpha);
    }

    public function testpowOverflowBufferXwithOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX,$alpha] =
            $this->translate_maximum($X,2);

        $offX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->pow($N,$XX,$offX,$incX,$alpha);
    }

    public function testpowOverflowBufferXwithIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX,$alpha] =
            $this->translate_maximum($X,2);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->pow($N,$XX,$offX,$incX,$alpha);
    }

    public function testexpNormal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,2,4,9]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $math->exp($N,$XX,$offX,$incX);
        $math->log($N,$XX,$offX,$incX);

        $this->assertEquals([0,2,4,9],$X->toArray());
    }

    public function testexpMinusN()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $N = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->exp($N,$XX,$offX,$incX);
    }

    public function testexpMinusOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $offX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->exp($N,$XX,$offX,$incX);
    }

    public function testexpMinusIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $incX = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->exp($N,$XX,$offX,$incX);
    }

    public function testexpIllegalBufferX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->exp($N,$XX,$offX,$incX);
    }

    public function testexpOverflowBufferXwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $XX = $mo->array([1,2])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->exp($N,$XX,$offX,$incX);
    }

    public function testexpOverflowBufferXwithOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $offX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->exp($N,$XX,$offX,$incX);
    }

    public function testexpOverflowBufferXwithIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->exp($N,$XX,$offX,$incX);
    }

    public function testlogNormal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,4,9]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $math->log($N,$XX,$offX,$incX);
        $math->exp($N,$XX,$offX,$incX);

        $this->assertEquals([1,2,4,9],$X->toArray());
    }

    public function testlogInvalidValue()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,0,-1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        // *** CAUTION ***
        // disable checking for INFINITY values
        //$this->expectException(RuntimeException::class);
        //$this->expectExceptionMessage('Invalid value in log.');

        $math->log($N,$XX,$offX,$incX);
        $this->assertEquals(0.0,  $X[0]);
        $this->assertEquals(-INF, $X[1]);
        $this->assertTrue(is_nan($X[2]));
    }

    public function testlogMinusN()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $N = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->log($N,$XX,$offX,$incX);
    }

    public function testlogMinusOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $offX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->log($N,$XX,$offX,$incX);
    }

    public function testlogMinusIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $incX = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->log($N,$XX,$offX,$incX);
    }

    public function testlogIllegalBufferX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->log($N,$XX,$offX,$incX);
    }

    public function testlogOverflowBufferXwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $XX = $mo->array([1,2])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->log($N,$XX,$offX,$incX);
    }

    public function testlogOverflowBufferXwithOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $offX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->log($N,$XX,$offX,$incX);
    }

    public function testlogOverflowBufferXwithIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->log($N,$XX,$offX,$incX);
    }

    public function testnan2numNormal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([NAN,2,4,NAN]);
        [$N,$XX,$offX,$incX,$alpha] =
            $this->translate_maximum($X,0.0);
        $math->nan2num($N,$XX,$offX,$incX,$alpha);
        $this->assertEquals([0,2,4,0],$X->toArray());

        $X = $mo->array([NAN,2,4,NAN]);
        [$N,$XX,$offX,$incX,$alpha] =
            $this->translate_maximum($X,1.0);
        $math->nan2num($N,$XX,$offX,$incX,$alpha);
        $this->assertEquals([1,2,4,1],$X->toArray());
    }

    public function testisnanNormal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([NAN,2,4,NAN]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $math->isnan($N,$XX,$offX,$incX);

        $this->assertEquals([1,0,0,1],$X->toArray());
    }

    public function testsearchsortedNormal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([0,2,4]);
        $X = $mo->array([-1,1,2,5]);
        $Y = $mo->zeros([4],NDArray::int32);
        [$m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $math->searchsorted($m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
        $this->assertEquals([0,1,1,3],$Y->toArray());

        $math->searchsorted($m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,true,$YY,$offsetY,$incY);
        $this->assertEquals([0,1,2,3],$Y->toArray());
    }

    public function testsearchsortedMinusM()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([0,2,4]);
        $X = $mo->array([-1,1,2,5]);
        $Y = $mo->zeros([4],NDArray::int32);
        [$m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $m = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->searchsorted($m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedMinusOffsetA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([0,2,4]);
        $X = $mo->array([-1,1,2,5]);
        $Y = $mo->zeros([4],NDArray::int32);
        [$m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $offsetA = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetA must be greater than equals 0.');
        $math->searchsorted($m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedMinusIncA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([0,2,4]);
        $X = $mo->array([-1,1,2,5]);
        $Y = $mo->zeros([4],NDArray::int32);
        [$m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $incA = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incA must be greater than 0.');
        $math->searchsorted($m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedIllegalBufferA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([0,2,4]);
        $X = $mo->array([-1,1,2,5]);
        $Y = $mo->zeros([4],NDArray::int32);
        [$m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $AA = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('A must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->searchsorted($m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedOverflowBufferAwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([0,2,4]);
        $X = $mo->array([-1,1,2,5]);
        $Y = $mo->zeros([4],NDArray::int32);
        [$m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $AA = $mo->array([1,2])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferA');
        $math->searchsorted($m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedOverflowBufferAwithOffsetA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([0,2,4]);
        $X = $mo->array([-1,1,2,5]);
        $Y = $mo->zeros([4],NDArray::int32);
        [$m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $offsetA = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferA');
        $math->searchsorted($m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedOverflowBufferAwithIncA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([0,2,4]);
        $X = $mo->array([-1,1,2,5]);
        $Y = $mo->zeros([4],NDArray::int32);
        [$m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $incA = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferA');
        $math->searchsorted($m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedMinusN()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([0,2,4]);
        $X = $mo->array([-1,1,2,5]);
        $Y = $mo->zeros([4],NDArray::int32);
        [$m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $n = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->searchsorted($m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedMinusOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([0,2,4]);
        $X = $mo->array([-1,1,2,5]);
        $Y = $mo->zeros([4],NDArray::int32);
        [$m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $offsetX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->searchsorted($m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedMinusIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([0,2,4]);
        $X = $mo->array([-1,1,2,5]);
        $Y = $mo->zeros([4],NDArray::int32);
        [$m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $incX = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->searchsorted($m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedIllegalBufferX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([0,2,4]);
        $X = $mo->array([-1,1,2,5]);
        $Y = $mo->zeros([4],NDArray::int32);
        [$m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('X must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->searchsorted($m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedOverflowBufferXwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([0,2,4]);
        $X = $mo->array([-1,1,2,5]);
        $Y = $mo->zeros([4],NDArray::int32);
        [$m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $XX = $mo->array([1,2])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->searchsorted($m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedOverflowBufferXwithOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([0,2,4]);
        $X = $mo->array([-1,1,2,5]);
        $Y = $mo->zeros([4],NDArray::int32);
        [$m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $offsetX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->searchsorted($m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedOverflowBufferXwithIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([0,2,4]);
        $X = $mo->array([-1,1,2,5]);
        $Y = $mo->zeros([4],NDArray::int32);
        [$m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->searchsorted($m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedMinusOffsetY()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([0,2,4]);
        $X = $mo->array([-1,1,2,5]);
        $Y = $mo->zeros([4],NDArray::int32);
        [$m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $offsetY = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetY must be greater than equals 0.');
        $math->searchsorted($m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedMinusIncY()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([0,2,4]);
        $X = $mo->array([-1,1,2,5]);
        $Y = $mo->zeros([4],NDArray::int32);
        [$m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $incY = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incY must be greater than 0.');
        $math->searchsorted($m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedIllegalBufferY()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([0,2,4]);
        $X = $mo->array([-1,1,2,5]);
        $Y = $mo->zeros([4],NDArray::int32);
        [$m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $YY = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('Y must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->searchsorted($m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedOverflowBufferYwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([0,2,4]);
        $X = $mo->array([-1,1,2,5]);
        $Y = $mo->zeros([4],NDArray::int32);
        [$m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $YY = $mo->array([1,2])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferY');
        $math->searchsorted($m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedOverflowBufferYwithOffsetY()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([0,2,4]);
        $X = $mo->array([-1,1,2,5]);
        $Y = $mo->zeros([4],NDArray::int32);
        [$m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $offsetX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->searchsorted($m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedOverflowBufferYwithIncY()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([0,2,4]);
        $X = $mo->array([-1,1,2,5]);
        $Y = $mo->zeros([4],NDArray::int32);
        [$m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $incY = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferY');
        $math->searchsorted($m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedUnmatchDataType()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([0,2,4],NDArray::float32);
        $X = $mo->array([-1,1,2,5],NDArray::float64);
        $Y = $mo->zeros([4],NDArray::int32);
        [$m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Unmatch data type for A and X');
        $math->searchsorted($m,$AA,$offsetA,$incA,$n,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

//=========================================================================

    public function testcumsumNormal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $Y = $mo->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
        $this->assertEquals([1,3,6],$Y->toArray());

        $math->cumsum($n,$XX,$offsetX,$incX,true,$reverse,$YY,$offsetY,$incY);
        $this->assertEquals([0,1,3],$Y->toArray());

        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,true,$YY,$offsetY,$incY);
        $this->assertEquals([6,3,1],$Y->toArray());

        $math->cumsum($n,$XX,$offsetX,$incX,true,true,$YY,$offsetY,$incY);
        $this->assertEquals([3,1,0],$Y->toArray());

        $X = $mo->array([1,NAN,3]);
        $Y = $mo->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
        $this->assertEquals(1,$Y[0]);
        $this->assertTrue(is_nan($Y[1]));
        $this->assertTrue(is_nan($Y[2]));

        $math->cumsum($n,$XX,$offsetX,$incX,true,$reverse,$YY,$offsetY,$incY);
        $this->assertEquals(0,$Y[0]);
        $this->assertEquals(1,$Y[1]);
        $this->assertTrue(is_nan($Y[2]));

        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,true,$YY,$offsetY,$incY);
        $this->assertTrue(is_nan($Y[0]));
        $this->assertTrue(is_nan($Y[1]));
        $this->assertEquals(1,$Y[2]);
    }

    public function testcumsumMinusN()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $Y = $mo->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $n = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumMinusOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $Y = $mo->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $offsetX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumMinusIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $Y = $mo->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $incX = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumIllegalBufferX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $Y = $mo->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('X must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumOverflowBufferXwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $Y = $mo->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $XX = $mo->array([1,2])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumOverflowBufferXwithOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $Y = $mo->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $offsetX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumOverflowBufferXwithIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $Y = $mo->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumMinusOffsetY()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $Y = $mo->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $offsetY = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetY must be greater than equals 0.');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumMinusIncY()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $Y = $mo->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $incY = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incY must be greater than 0.');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumIllegalBufferY()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $Y = $mo->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $YY = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('Y must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumOverflowBufferYwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $Y = $mo->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $YY = $mo->array([1,2])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferY');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumOverflowBufferYwithOffsetY()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $Y = $mo->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $offsetX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumOverflowBufferYwithIncY()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $Y = $mo->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $incY = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferY');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumUnmatchDataType()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        $Y = $mo->zeros([3],NDArray::float64);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Unmatch data type for X and Y');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

//=========================================================================

    public function testzerosNormal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,4,9]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $math->zeros($N,$XX,$offX,$incX);

        $this->assertEquals([0,0,0,0],$X->toArray());
    }

    public function testzerosMinusN()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $N = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->zeros($N,$XX,$offX,$incX);
    }

    public function testzerosMinusOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $offX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->zeros($N,$XX,$offX,$incX);
    }

    public function testzerosMinusIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $incX = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->zeros($N,$XX,$offX,$incX);
    }

    public function testzerosIllegalBufferX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->zeros($N,$XX,$offX,$incX);
    }

    public function testzerosOverflowBufferXwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $XX = $mo->array([1,2])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->zeros($N,$XX,$offX,$incX);
    }

    public function testzerosOverflowBufferXwithOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $offX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->zeros($N,$XX,$offX,$incX);
    }

    public function testzerosOverflowBufferXwithIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->zeros($N,$XX,$offX,$incX);
    }

######################################################################

    public function testGatherAxisNullNormal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,2],NDArray::int32);
        $B = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([[1,2,3],[7,8,9]],$B->toArray());

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float64);
        $X = $mo->array([0,2],NDArray::int64);
        $B = $mo->array([[0,0,0],[0,0,0]],NDArray::float64);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([[1,2,3],[7,8,9]],$B->toArray());

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::int64);
        $X = $mo->array([0,2],NDArray::int64);
        $B = $mo->array([[0,0,0],[0,0,0]],NDArray::int64);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([[1,2,3],[7,8,9]],$B->toArray());

        $A = $mo->array([1,2,3,4],NDArray::float32);
        $X = $mo->array([0,2],NDArray::int32);
        $B = $mo->array([0,0],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([1,3],$B->toArray());

        $A = $mo->array([1,2,3,4],NDArray::float64);
        $X = $mo->array([0,2],NDArray::int64);
        $B = $mo->array([0,0],NDArray::float64);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([1,3],$B->toArray());

        $A = $mo->array([1,2,3,4],NDArray::int64);
        $X = $mo->array([0,2],NDArray::int64);
        $B = $mo->array([0,0],NDArray::int64);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([1,3],$B->toArray());
    }


    public function testGatherAxisNullLabelNumberOutOfBounds1()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,4],NDArray::int32);
        $B = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Label number is out of bounds.');
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testGatherAxisNullLabelNumberOutOfBounds2()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,-1],NDArray::int32);
        $B = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Label number is out of bounds.');
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testGatherAxisNullMinusN()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $B = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $n = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testGatherAxisNullMinusK()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $B = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $k = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument k must be greater than 0.');
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testGatherAxisNullMinusNumClass()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $B = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $numClass = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument numClass must be greater than or equal 0.');
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testGatherAxisNullMinusOffsetA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $B = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $offA = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetA must be greater than or equal 0.');
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testGatherAxisNullIllegalBufferA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $B = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $AA = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testGatherAxisNullOverflowBufferAwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $B = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $AA = $mo->array([1,2,3,4,5,6,7,8,9,10,11])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix A specification too large for buffer');
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testGatherAxisNullOverflowBufferAwithOffset()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $B = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $offA = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix A specification too large for buffer');
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testGatherAxisNullMinusOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $B = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equal 0.');
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testGatherAxisNullIllegalBufferX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $B = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testGatherAxisNullOverflowBufferXwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $B = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $XX = $mo->array([1])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix X specification too large for buffer.');
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testGatherAxisNullOverflowBufferXwithOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $B = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix X specification too large for buffer.');
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testGatherAxisNullMinusOffsetB()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $B = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $offB = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetB must be greater than or equal 0.');
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testGatherAxisNullIllegalBufferB()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $B = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $BB = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testGatherAxisNullOverflowBufferBwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $B = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $BB = $mo->array([0])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix B specification too large for buffer.');
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testGatherAxisNullOverflowBufferBwithOffsetB()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $B = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $offB = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix B specification too large for buffer.');
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

######################################################################
    public function testReduceGatherAxis1Normal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,2]);
        $B = $mo->array([0,0]);
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=1,$B,$A->dtype());
        $this->assertTrue($reduce);

        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([2,6],$B->toArray());
    }

    public function testReduceGatherAxis1LabelNumberOutOfBounds1()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,3]);
        $B = $mo->array([0,0]);
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=1,$B,$A->dtype());
        $this->assertTrue($reduce);

        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Label number is out of bounds.');
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testReduceGatherAxis1LabelNumberOutOfBounds2()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $B = $mo->array([0,0]);
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=1,$B,$A->dtype());
        $this->assertTrue($reduce);

        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Label number is out of bounds.');
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testReduceGatherAxis1MinusM()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $B = $mo->array([0,0]);
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=1,$B,$A->dtype());
        $this->assertTrue($reduce);

        $m = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testReduceGatherAxis1MinusN()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $B = $mo->array([0,0]);
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=1,$B,$A->dtype());
        $this->assertTrue($reduce);

        $n = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testReduceGatherAxis1MinusOffsetA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $B = $mo->array([0,0]);
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=1,$B,$A->dtype());
        $this->assertTrue($reduce);

        $offA = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetA must be greater than or equal 0.');
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testReduceGatherAxis1IllegalBufferA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $B = $mo->array([0,0]);
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=1,$B,$A->dtype());
        $this->assertTrue($reduce);

        $AA = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testReduceGatherAxis1OverflowBufferAwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $B = $mo->array([0,0]);
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=1,$B,$A->dtype());
        $this->assertTrue($reduce);

        $AA = $mo->array([1,2,3,4,5])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix A specification too large for buffer.');
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testReduceGatherAxis1OverflowBufferAwithOffset()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $B = $mo->array([0,0]);
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=1,$B,$A->dtype());
        $this->assertTrue($reduce);

        $offA = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix A specification too large for buffer.');
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testReduceGatherAxis1MinusOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $B = $mo->array([0,0]);
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=1,$B,$A->dtype());
        $this->assertTrue($reduce);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equal 0.');
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testReduceGatherAxis1IllegalBufferX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $B = $mo->array([0,0]);
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=1,$B,$A->dtype());
        $this->assertTrue($reduce);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testReduceGatherAxis1OverflowBufferXwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $B = $mo->array([0,0]);
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=1,$B,$A->dtype());
        $this->assertTrue($reduce);
        $XX = $mo->array([1])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix X specification too large for buffer.');
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testReduceGatherAxis1OverflowBufferXwithOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $B = $mo->array([0,0]);
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=1,$B,$A->dtype());
        $this->assertTrue($reduce);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix X specification too large for buffer.');
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testReduceGatherAxis1MinusOffsetB()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $B = $mo->array([0,0]);
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=1,$B,$A->dtype());
        $this->assertTrue($reduce);

        $offB = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetB must be greater than or equal 0.');
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testReduceGatherAxis1IllegalBufferB()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $B = $mo->array([0,0]);
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=1,$B,$A->dtype());
        $this->assertTrue($reduce);

        $BB = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testReduceGatherAxis1OverflowBufferBwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $B = $mo->array([0,0]);
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=1,$B,$A->dtype());
        $this->assertTrue($reduce);

        $BB = $mo->array([0])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix B specification too large for buffer.');
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testReduceGatherAxis1OverflowBufferXwithOffsetB()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $B = $mo->array([0,0]);
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=1,$B,$A->dtype());
        $this->assertTrue($reduce);

        $offB = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix B specification too large for buffer.');
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testScatterAxisNull()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);
        // float32
        $numClass = 4;
        $X = $mo->array([0,2],NDArray::int64);
        $A = $mo->array([[1,2,3],[7,8,9]],NDArray::float32);
        $B = $mo->zeros([4,3],$A->dtype());
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA]
            = $this->translate_scatter($X,$A,$numClass,$axis=null,$B);
        $this->assertFalse($reduce);
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA);

        $this->assertEquals(
           [[1,2,3],
            [0,0,0],
            [7,8,9],
            [0,0,0]],
            $B->toArray()
        );

        // float64
        $X = $mo->array([0,2],NDArray::int64);
        $A = $mo->array([[1,2,3],[7,8,9]],NDArray::float64);
        $B = $mo->zeros([4,3],$A->dtype());
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA]
            = $this->translate_scatter($X,$A,$numClass,$axis=null,$B);
        $this->assertFalse($reduce);
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA);

        $this->assertEquals(
           [[1,2,3],
            [0,0,0],
            [7,8,9],
            [0,0,0]],
            $B->toArray()
        );
        // int64
        $X = $mo->array([0,2],NDArray::int64);
        $A = $mo->array([[1,2,3],[7,8,9]],NDArray::int64);
        $B = $mo->zeros([4,3],$A->dtype());
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA]
            = $this->translate_scatter($X,$A,$numClass,$axis=null,$B);
        $this->assertFalse($reduce);
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA);
        $this->assertEquals(
           [[1,2,3],
            [0,0,0],
            [7,8,9],
            [0,0,0]],
            $B->toArray()
        );
        // uint8
        $X = $mo->array([0,2],NDArray::int64);
        $A = $mo->array([[1,2,3],[7,8,9]],NDArray::uint8);
        $B = $mo->zeros([4,3],$A->dtype());
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA]
            = $this->translate_scatter($X,$A,$numClass,$axis=null,$B);
        $this->assertFalse($reduce);
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA);
        $this->assertEquals(
           [[1,2,3],
            [0,0,0],
            [7,8,9],
            [0,0,0]],
            $B->toArray()
        );
        // float32
        $X = $mo->array([0,2],NDArray::int64);
        $A = $mo->array([1,3],NDArray::float32);
        $B = $mo->zeros([4],$A->dtype());
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA]
            = $this->translate_scatter($X,$A,$numClass,$axis=null,$B);
        $this->assertFalse($reduce);
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA);
        $this->assertEquals(
           [1,0,3,0],
            $B->toArray()
        );
        // int32
        $X = $mo->array([0,2],NDArray::int64);
        $A = $mo->array([1,3],NDArray::int32);
        $B = $mo->zeros([4],$A->dtype());
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA]
            = $this->translate_scatter($X,$A,$numClass,$axis=null,$B);
        $this->assertFalse($reduce);
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA);
        $this->assertEquals(
           [1,0,3,0],
            $B->toArray()
        );
        // float64
        $X = $mo->array([0,2],NDArray::int64);
        $A = $mo->array([1,3],NDArray::float64);
        $B = $mo->zeros([4],$A->dtype());
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA]
            = $this->translate_scatter($X,$A,$numClass,$axis=null,$B);
        $this->assertFalse($reduce);
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA);
        $this->assertEquals(
           [1,0,3,0],
            $B->toArray()
        );
        // int64
        $X = $mo->array([0,2],NDArray::int64);
        $A = $mo->array([1,3],NDArray::int64);
        $B = $mo->zeros([4],$A->dtype());
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA]
            = $this->translate_scatter($X,$A,$numClass,$axis=null,$B);
        $this->assertFalse($reduce);
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA);
        $this->assertEquals(
           [1,0,3,0],
            $B->toArray()
        );
        // uint8
        $X = $mo->array([0,2],NDArray::int64);
        $A = $mo->array([252,254],NDArray::uint8);
        $B = $mo->zeros([4],$A->dtype());
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA]
            = $this->translate_scatter($X,$A,$numClass,$axis=null,$B);
        $this->assertFalse($reduce);
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA);
        $this->assertEquals(
           [252,0,254,0],
            $B->toArray()
        );
        // x=uint8
        $X = $mo->array([0,255],NDArray::uint8);
        $A = $mo->array([252,254],NDArray::uint8);
        $B = $mo->zeros([256],$A->dtype());
        $numClass = 256;
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA]
            = $this->translate_scatter($X,$A,$numClass,$axis=null,$B);
        $this->assertFalse($reduce);
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA);
        $this->assertEquals(252,$B[0]);
        $this->assertEquals(254,$B[255]);
    }

    public function testScatterAxis1()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);
        $numClass = 3;
        $X = $mo->array([0,1,2,0],NDArray::int32);
        $A = $mo->array([1,5,9,10],NDArray::float32);
        $B = $mo->zeros([4,3],$A->dtype());
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$BB,$offB,$AA,$offA]
            = $this->translate_scatter($X,$A,$numClass,$axis=1,$B);
        $this->assertTrue($reduce);
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$BB,$offB,$AA,$offA);

        $this->assertEquals(
           [[1,0,0],
            [0,5,0],
            [0,0,9],
            [10,0,0]],
            $B->toArray());

        $X = $mo->array([0,1,2,0],NDArray::int64);
        $B = $mo->zeros([4,3],$A->dtype());
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$BB,$offB,$AA,$offA]
            = $this->translate_scatter($X,$A,$numClass,$axis=1,$B);
        $this->assertTrue($reduce);
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$BB,$offB,$AA,$offA);

        $this->assertEquals(
           [[1,0,0],
            [0,5,0],
            [0,0,9],
            [10,0,0]],
            $B->toArray());

        $X = $mo->array([0,1,2,0],NDArray::float32);
        $B = $mo->zeros([4,3],$A->dtype());
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$BB,$offB,$AA,$offA]
            = $this->translate_scatter($X,$A,$numClass,$axis=1,$B);
        $this->assertTrue($reduce);
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$BB,$offB,$AA,$offA);
        $this->assertEquals(
           [[1,0,0],
            [0,5,0],
            [0,0,9],
            [10,0,0]],
            $B->toArray());

        $X = $mo->array([0,1,2,0],NDArray::float64);
        $B = $mo->zeros([4,3],$A->dtype());
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$BB,$offB,$AA,$offA]
            = $this->translate_scatter($X,$A,$numClass,$axis=1,$B);
        $this->assertTrue($reduce);
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$BB,$offB,$AA,$offA);
        $this->assertEquals(
           [[1,0,0],
            [0,5,0],
            [0,0,9],
            [10,0,0]],
            $B->toArray());
    }

    public function testupdateAddOnehotNormal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);
        $X = $mo->array([1, 2]);
        $Y = $mo->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
        $this->assertEquals([[10,9,10],[10,10,9]],$Y->toArray());
    }

    public function testupdateAddOnehotOutOfboundsLabelNumber1()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);
        $X = $mo->array([1, 3]);
        $Y = $mo->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Label number is out of bounds.');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotOutOfboundsLabelNumber2()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);
        $X = $mo->array([1, -1]);
        $Y = $mo->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Label number is out of bounds.');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotMinusM()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);
        $X = $mo->array([1, -1]);
        $Y = $mo->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $m = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

     public function testupdateAddOnehotMinusN()
     {
         $mo = new MatrixOperator();
         $math = $this->getMath($mo);
         $X = $mo->array([1, -1]);
         $Y = $mo->array([[10,10,10],[10,10,10]]);
         $numClass = 3;
         [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
             $X,$numClass,-1,$Y);

         $n = 0;
         $this->expectException(RuntimeException::class);
         $this->expectExceptionMessage('Argument n must be greater than 0.');
         $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
     }

    public function testupdateAddOnehotMinusOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);
        $X = $mo->array([1, -1]);
        $Y = $mo->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $offX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotMinusIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);
        $X = $mo->array([1, -1]);
        $Y = $mo->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $incX = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotIllegalBufferX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);
        $X = $mo->array([1, -1]);
        $Y = $mo->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotOverflowBufferXwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);
        $X = $mo->array([1, -1]);
        $Y = $mo->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $XX = $mo->array([1])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotOverflowBufferXwithOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);
        $X = $mo->array([1, -1]);
        $Y = $mo->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $offX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotOverflowBufferXwithIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);
        $X = $mo->array([1, -1]);
        $Y = $mo->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotMinusOffsetY()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);
        $X = $mo->array([1, -1]);
        $Y = $mo->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $offY = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetY must be greater than equals 0.');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotMinusLdY()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);
        $X = $mo->array([1, -1]);
        $Y = $mo->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $ldY = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument ldY must be greater than 0.');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotIllegalBufferY()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);
        $X = $mo->array([1, -1]);
        $Y = $mo->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $YY = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotOverflowBufferYwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);
        $X = $mo->array([1, -1]);
        $Y = $mo->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $YY = $mo->array([1])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferY');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotOverflowBufferYwithOffsetY()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);
        $X = $mo->array([1, -1]);
        $Y = $mo->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $offY = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferY');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotOverflowBufferYwithIncY()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);
        $X = $mo->array([1, -1]);
        $Y = $mo->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $ldY = 4;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferY');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testreduceSumSameSizeNormal()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([0,0]);
        [$m,$n,$k,$AA,$offA,$BB,$offB] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $math->reduceSum($m,$n,$k,$AA,$offA,$BB,$offB);
        $this->assertEquals([6,15],$X->toArray());
    }

    public function testreduceSumBroadcastTranspose()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([0,0,0]);
        [$m,$n,$k,$AA,$offA,$BB,$offB] =
            $this->translate_reduceSum($A,$axis=0,$X);

        $math->reduceSum($m,$n,$k,$AA,$offA,$BB,$offB);
        $this->assertEquals([5,7,9],$X->toArray());
    }

    public function testreduceSumMinusM()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,0]);
        $A = $mo->array([[1,2,3],[4,5,6]]);
        [$m,$n,$k,$AA,$offA,$BB,$offB] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $m = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->reduceSum($m,$n,$k,$AA,$offA,$BB,$offB);
    }

    public function testreduceSumMinusK()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,0]);
        $A = $mo->array([[1,2,3],[4,5,6]]);
        [$m,$n,$k,$AA,$offA,$BB,$offB] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $k = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument k must be greater than 0.');
        $math->reduceSum($m,$n,$k,$AA,$offA,$BB,$offB);
    }

    public function testreduceSumMinusN()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,0]);
        $A = $mo->array([[1,2,3],[4,5,6]]);
        [$m,$n,$k,$AA,$offA,$BB,$offB] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $n = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->reduceSum($m,$n,$k,$AA,$offA,$BB,$offB);
    }

    public function testreduceSumMinusOffsetB()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,0]);
        $A = $mo->array([[1,2,3],[4,5,6]]);
        [$m,$n,$k,$AA,$offA,$BB,$offB] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $offB = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetB must be greater than or equals 0.');
        $math->reduceSum($m,$n,$k,$AA,$offA,$BB,$offB);
    }

    public function testreduceSumIllegalBufferB()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,0]);
        $A = $mo->array([[1,2,3],[4,5,6]]);
        [$m,$n,$k,$AA,$offA,$BB,$offB] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $BB = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->reduceSum($m,$n,$k,$AA,$offA,$BB,$offB);
    }

    public function testreduceSumOverflowBufferBwithSize()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,0]);
        $A = $mo->array([[1,2,3],[4,5,6]]);
        [$m,$n,$k,$AA,$offA,$BB,$offB] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $BB = $mo->array([1])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferB');
        $math->reduceSum($m,$n,$k,$AA,$offA,$BB,$offB);
    }

    public function testreduceSumOverflowBufferBwithOffsetB()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,0]);
        $A = $mo->array([[1,2,3],[4,5,6]]);
        [$m,$n,$k,$AA,$offA,$BB,$offB] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $offB = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferB');
        $math->reduceSum($m,$n,$k,$AA,$offA,$BB,$offB);
    }

    public function testreduceSumMinusOffsetA()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,0]);
        $A = $mo->array([[1,2,3],[4,5,6]]);
        [$m,$n,$k,$AA,$offA,$BB,$offB] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $offA = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetA must be greater than or equals 0.');
        $math->reduceSum($m,$n,$k,$AA,$offA,$BB,$offB);
    }

    public function testreduceSumIllegalBufferA()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,0]);
        $A = $mo->array([[1,2,3],[4,5,6]]);
        [$m,$n,$k,$AA,$offA,$BB,$offB] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $AA = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->reduceSum($m,$n,$k,$AA,$offA,$BB,$offB);
    }

    public function testreduceSumOverflowBufferAwithSize()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,0]);
        $A = $mo->array([[1,2,3],[4,5,6]]);
        [$m,$n,$k,$AA,$offA,$BB,$offB] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $AA = $mo->array([1,2,3,4,5])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->reduceSum($m,$n,$k,$AA,$offA,$BB,$offB);
    }

    public function testreduceSumOverflowBufferXwithOffsetA()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,0]);
        $A = $mo->array([[1,2,3],[4,5,6]]);
        [$m,$n,$k,$AA,$offA,$BB,$offB] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $offA = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->reduceSum($m,$n,$k,$AA,$offA,$BB,$offB);
    }

    public function testsoftmax()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);
        $X = $mo->array([-1.0,-0.5,0.0,0.5,1.0]);
        $m = 1;
        $n = 5;
        $XX = $X->buffer();
        $offX = 0;
        $ldX = 5;
        $math->softmax($m,$n,$XX,$offX,$ldX);

        $this->assertTrue($X[0]>0.0);
        $this->assertTrue($X[0]<$X[1]);
        $this->assertTrue($X[1]<$X[2]);
        $this->assertTrue($X[2]<$X[3]);
        $this->assertTrue($X[3]<$X[4]);
        $this->assertTrue($X[4]<1.0);
        $this->assertTrue(1.0e-5>abs(1.0-$this->sum($n,$X,$offX,1)));
        $single = $X->toArray();

        // batch mode
        $y = $mo->array([
            [-1.0,-0.5,0.0,0.5,1.0],
            [-1.0,-0.5,0.0,0.5,1.0],
            [-1.0,-0.5,0.0,0.5,1.0],
            [-1.0,-0.5,0.0,0.5,1.0],
            [-1.0,-0.5,0.0,0.5,1.0],
        ]);

        $m = 5;
        $n = 5;
        $YY = $y->buffer();
        $offY = 0;
        $ldY = 5;
        $math->softmax($m,$n,$YY,$offY,$ldY);

        $this->assertEquals($single,$y[0]->toArray());
        $this->assertEquals($single,$y[1]->toArray());
        $this->assertEquals($single,$y[2]->toArray());
        $this->assertEquals($single,$y[3]->toArray());
        $this->assertEquals($single,$y[4]->toArray());
    }

    public function testequal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);
        $n = 5;
        $offX = 0;
        $incX = 1;
        $offY = 0;
        $incY = 1;

        $X = $mo->array([-1.0,-0.5,0.0,0.5,-1.0],NDArray::float32);
        $Y = $mo->array([1.0,-0.5,0.0,0.5,1.0],NDArray::float32);
        $XX = $X->buffer();
        $YY = $Y->buffer();
        $math->equal($n,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([0,1,1,1,0],$Y->toArray());

        $X = $mo->array([-1.0,-0.5,0.0,0.5,-1.0],NDArray::float64);
        $Y = $mo->array([1.0,-0.5,0.0,0.5,1.0],NDArray::float64);
        $XX = $X->buffer();
        $YY = $Y->buffer();
        $math->equal($n,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([0,1,1,1,0],$Y->toArray());

        $X = $mo->array([-1.0,-0.5,0.0,0.5,-1.0],NDArray::int8);
        $Y = $mo->array([1.0,-0.5,0.0,0.5,1.0],NDArray::int8);
        $XX = $X->buffer();
        $YY = $Y->buffer();
        $math->equal($n,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([0,1,1,1,0],$Y->toArray());

        $X = $mo->array([-1.0,-0.5,0.0,0.5,-1.0],NDArray::int16);
        $Y = $mo->array([1.0,-0.5,0.0,0.5,1.0],NDArray::int16);
        $XX = $X->buffer();
        $YY = $Y->buffer();
        $math->equal($n,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([0,1,1,1,0],$Y->toArray());

        $X = $mo->array([-1.0,-0.5,0.0,0.5,-1.0],NDArray::int32);
        $Y = $mo->array([1.0,-0.5,0.0,0.5,1.0],NDArray::int32);
        $XX = $X->buffer();
        $YY = $Y->buffer();
        $math->equal($n,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([0,1,1,1,0],$Y->toArray());

        $X = $mo->array([-1.0,-0.5,0.0,0.5,-1.0],NDArray::int64);
        $Y = $mo->array([1.0,-0.5,0.0,0.5,1.0],NDArray::int64);
        $XX = $X->buffer();
        $YY = $Y->buffer();
        $math->equal($n,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([0,1,1,1,0],$Y->toArray());

        $X = $mo->array([false,false,true ,true,true ],NDArray::bool);
        $Y = $mo->array([true ,false,true ,true,false],NDArray::bool);
        $XX = $X->buffer();
        $YY = $Y->buffer();
        $math->equal($n,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([false,true,true,true,false],$Y->toArray());
    }


    public function testastype()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        #### int to any
        $X = $mo->array([-1,0,1,2,3],NDArray::int32);
        $dtype = NDArray::float32;
        $Y = $mo->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());

        $dtype = NDArray::float64;
        $Y = $mo->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());

        $dtype = NDArray::int8;
        $Y = $mo->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());

        $dtype = NDArray::int16;
        $Y = $mo->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());

        $dtype = NDArray::int32;
        $Y = $mo->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());

        $dtype = NDArray::int64;
        $Y = $mo->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());

        $dtype = NDArray::bool;
        $Y = $mo->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([true,false,true,true,true],$Y->toArray());

        #### float to any ######
        $X = $mo->array([-1,0,1,2,3],NDArray::float32);
        $dtype = NDArray::float32;
        $Y = $mo->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());

        $dtype = NDArray::float64;
        $Y = $mo->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());

        $dtype = NDArray::int8;
        $Y = $mo->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());

        $dtype = NDArray::int16;
        $Y = $mo->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());

        $dtype = NDArray::int32;
        $Y = $mo->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());

        $dtype = NDArray::int64;
        $Y = $mo->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());

        $dtype = NDArray::bool;
        $Y = $mo->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([true,false,true,true,true],$Y->toArray());

        #### bool to any ######
        $X = $mo->array([true,false,true,true,true],NDArray::bool);
        $dtype = NDArray::float32;
        $Y = $mo->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([1,0,1,1,1],$Y->toArray());

        $dtype = NDArray::float64;
        $Y = $mo->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([1,0,1,1,1],$Y->toArray());

        $dtype = NDArray::int8;
        $Y = $mo->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([1,0,1,1,1],$Y->toArray());

        $dtype = NDArray::int16;
        $Y = $mo->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([1,0,1,1,1],$Y->toArray());

        $dtype = NDArray::int32;
        $Y = $mo->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([1,0,1,1,1],$Y->toArray());

        $dtype = NDArray::int64;
        $Y = $mo->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([1,0,1,1,1],$Y->toArray());

        $dtype = NDArray::bool;
        $Y = $mo->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([true,false,true,true,true],$Y->toArray());

        #### float to unsigned ######
        $X = $mo->array([-1,0,1,2,3],NDArray::float32);
        $dtype = NDArray::uint8;
        $Y = $mo->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([255,0,1,2,3],$Y->toArray());

        $dtype = NDArray::uint16;
        $Y = $mo->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([65535,0,1,2,3],$Y->toArray());

        $dtype = NDArray::uint32;
        $Y = $mo->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([4294967295,0,1,2,3],$Y->toArray());

        $dtype = NDArray::uint64;
        $Y = $mo->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());
    }

    public function testreduceMaxSameSizeNormal()
    {
        if($this->checkSkip('reduceMax')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([0,0]);
        [$m,$n,$k,$AA,$offA,$BB,$offB] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $math->reduceMax($m,$n,$k,$AA,$offA,$BB,$offB);
        $this->assertEquals([3,6],$X->toArray());
    }

    public function testreduceArgMaxSameSizeNormal()
    {
        if($this->checkSkip('reduceArgMax')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([0,0],NDArray::float32);
        [$m,$n,$k,$AA,$offA,$BB,$offB] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $math->reduceArgMax($m,$n,$k,$AA,$offA,$BB,$offB);
        $this->assertEquals([2,2],$X->toArray());
    }


    public function testIm2col1dNormal()
    {
        if($this->checkSkip('im2col1d')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $images = $mo->array([1,2,3,4]);
        $cols = $mo->zeros([1,2,3,1],NDArray::float32);

        $images_offset = $images->offset();
        $images_size = $images->size();
        $images_buff = $images->buffer();
        $cols_buff = $cols->buffer();
        $cols_offset = $cols->offset();
        $cols_size = $cols->size();
        $math->im2col1d(
            $reverse=false,
            $images_buff,
            $images_offset,
            $images_size,
            $batches=1,
            $in_w=4,
            $channels=1,
            $filter_w=3,
            $stride_w=1,
            $padding=false,
            $channels_first=false,
            $dilation_w=1,
            $cols_channels_first=false,
            $cols_buff,
            $cols_offset,
            $cols_size
        );
        $this->assertEquals(
            [[[[1],[2],[3]],
              [[2],[3],[4]]]],
            $cols->toArray());
    }

    public function testIm2col2dNormal()
    {
        if($this->checkSkip('im2col2d')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $reverse = false;
        $batches = 1;
        $im_h = 4;
        $im_w = 4;
        $channels = 3;
        $kernel_h = 3;
        $kernel_w = 3;
        $stride_h = 1;
        $stride_w = 1;
        $padding = null;
        $channels_first = null;
        $cols_channels_first=null;
        $cols = null;
        $out_h = 2;
        $out_w = 2;
        $images = $mo->arange(
            $batches*
            $im_h*$im_w*
            $channels,
            null,null,
            NDArray::float32
        )->reshape([
            $batches,
            $im_h,
            $im_w,
            $channels
        ]);
        $cols = $mo->zeros(
            [
                $batches,
                $out_h,$out_w,
                $kernel_h,$kernel_w,
                $channels,
            ]);
        $images_offset = $images->offset();
        $images_size = $images->size();
        $images_buff = $images->buffer();
        $cols_buff = $cols->buffer();
        $cols_offset = $cols->offset();
        $cols_size = $cols->size();
        $math->im2col2d(
            $reverse,
            $images_buff,
            $images_offset,
            $images_size,
            $batches,
            $im_h,
            $im_w,
            $channels,
            $kernel_h,
            $kernel_w,
            $stride_h,
            $stride_w,
            $padding,
            $channels_first,
            $dilation_h=1,
            $dilation_w=1,
            $cols_channels_first,
            $cols_buff,
            $cols_offset,
            $cols_size
        );
        $this->assertEquals(
        [[
          [
           [[[0,1,2],[3,4,5],[6,7,8]],
            [[12,13,14],[15,16,17],[18,19,20]],
            [[24,25,26],[27,28,29],[30,31,32]],],
           [[[3,4,5],[6,7,8],[9,10,11]],
            [[15,16,17],[18,19,20],[21,22,23]],
            [[27,28,29],[30,31,32],[33,34,35]],],
          ],
          [
           [[[12,13,14],[15,16,17],[18,19,20]],
            [[24,25,26],[27,28,29],[30,31,32]],
            [[36,37,38],[39,40,41],[42,43,44]],],
           [[[15,16,17],[18,19,20],[21,22,23]],
            [[27,28,29],[30,31,32],[33,34,35]],
            [[39,40,41],[42,43,44],[45,46,47]],],
          ],
        ]],
        $cols->toArray()
        );
    }

    public function testIm2col3dNormal()
    {
        if($this->checkSkip('im2col3d')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $reverse = false;
        $batches = 1;
        $im_d = 4;
        $im_h = 4;
        $im_w = 4;
        $channels = 3;
        $kernel_d = 3;
        $kernel_h = 3;
        $kernel_w = 3;
        $stride_d = 1;
        $stride_h = 1;
        $stride_w = 1;
        $padding = null;
        $channels_first = null;
        $cols_channels_first=null;
        $cols = null;
        $out_d = 2;
        $out_h = 2;
        $out_w = 2;

        $images = $mo->arange(
            $batches*
            $im_d*$im_h*$im_w*
            $channels,
            null,null,
            NDArray::float32
        )->reshape([
            $batches,
            $im_d,
            $im_h,
            $im_w,
            $channels
        ]);

        $cols = $mo->zeros(
            [
                $batches,
                $out_d,$out_h,$out_w,
                $kernel_d,$kernel_h,$kernel_w,
                $channels,
            ]);
        $images_offset = $images->offset();
        $images_size = $images->size();
        $images_buff = $images->buffer();
        $cols_buff = $cols->buffer();
        $cols_offset = $cols->offset();
        $cols_size = $cols->size();
        $math->im2col3d(
            $reverse,
            $images_buff,
            $images_offset,
            $images_size,
            $batches,
            $im_d,
            $im_h,
            $im_w,
            $channels,
            $kernel_d,
            $kernel_h,
            $kernel_w,
            $stride_d,
            $stride_h,
            $stride_w,
            $padding,
            $channels_first,
            $dilation_d=1,
            $dilation_h=1,
            $dilation_w=1,
            $cols_channels_first,
            $cols_buff,
            $cols_offset,
            $cols_size
        );
        $this->assertNotEquals(
            $cols->toArray(),
            $mo->zerosLike($cols)
            );
    }

}
