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
        float $alpha,
        NDArray $X
        ) : array
    {
        $n = $X->size();
        $XX = $X->buffer();
        $offX = $X->offset();

        return [$n,$alpha,$XX,$offX,1];
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
       NDArray $X=null
       ) : array
   {
       if($axis===null)
           $axis = 0;
       if($axis!==0 && $axis!==1 && $axis!==-1)
           throw new InvalidArgumentException('"axis" must be 0 or 1 or -1.');
       if($axis===-1) {
           $axis = 1;
       }
       $shapeA = $A->shape();
       if($axis==0) {
           $trans = true;
           $rows = array_pop($shapeA);
       } else {
           $trans = false;
           $rows = $shapeA[0];
       }

       if($X==null) {
           $X = $this->alloc([$rows],$A->dtype());
       } else {
           if($X->shape()!=[$rows]) {
               $shapeError = '('.implode(',',$A->shape()).'),('.implode(',',$X->shape()).')';
               throw new InvalidArgumentException("Unmatch shape of dimension: ".$shapeError);
           }
       }

       $m = $A->shape()[0];
       $n = $A->size()/$m;
       $AA = $A->buffer();
       $offA = $A->offset();
       $XX = $X->buffer();
       $offX = $X->offset();

       return [
           $trans,
           $m,
           $n,
           $AA,$offA,$n,
           $XX,$offX,1
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
       $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
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
       $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
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
       $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
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
        $this->expectExceptionMessage('parameter 2 to be float');
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
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
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

        $X = $mo->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,0,0);

        // X := 1 / ( alpha * X + beta )
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Zero divide.');
        $math->reciprocal($N,$alpha,$XX,$offX,$incX,$beta);
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
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
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

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $math->maximum($N,$alpha,$XX,$offX,$incX);
        $this->assertEquals([2,2,3],$X->toArray());
    }

    public function testMaximumMinusN()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $N = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->maximum($N,$alpha,$XX,$offX,$incX);
    }

    public function testMaximumMinusOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $offX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->maximum($N,$alpha,$XX,$offX,$incX);
    }

    public function testMaximumMinusIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $incX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->maximum($N,$alpha,$XX,$offX,$incX);
    }

    public function testMaximumIllegalBufferX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
        $math->maximum($N,$alpha,$XX,$offX,$incX);
    }

    public function testMaximumOverflowBufferXwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $XX = $mo->array([1,2])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->maximum($N,$alpha,$XX,$offX,$incX);
    }

    public function testMaximumOverflowBufferXwithOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $offX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->maximum($N,$alpha,$XX,$offX,$incX);
    }

    public function testMaximumOverflowBufferXwithIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->maximum($N,$alpha,$XX,$offX,$incX);
    }

    public function testMinimumNormal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $math->minimum($N,$alpha,$XX,$offX,$incX);
        $this->assertEquals([1,2,2],$X->toArray());
    }

    public function testMinimumMinusN()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $N = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->minimum($N,$alpha,$XX,$offX,$incX);
    }

    public function testMinimumMinusOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $offX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->minimum($N,$alpha,$XX,$offX,$incX);
    }

    public function testMinimumMinusIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $incX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->minimum($N,$alpha,$XX,$offX,$incX);
    }

    public function testMinimumIllegalBufferX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
        $math->minimum($N,$alpha,$XX,$offX,$incX);
    }

    public function testMinimumOverflowBufferXwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $XX = $mo->array([1,2])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->minimum($N,$alpha,$XX,$offX,$incX);
    }

    public function testMinimumOverflowBufferXwithOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $offX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->minimum($N,$alpha,$XX,$offX,$incX);
    }

    public function testMinimumOverflowBufferXwithIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->minimum($N,$alpha,$XX,$offX,$incX);
    }

    public function testGreaterNormal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $math->greater($N,$alpha,$XX,$offX,$incX);
        $this->assertEquals([0,0,1],$X->toArray());
    }

    public function testGreaterMinusN()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $N = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->greater($N,$alpha,$XX,$offX,$incX);
    }

    public function testGreaterMinusOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $offX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->greater($N,$alpha,$XX,$offX,$incX);
    }

    public function testGreaterMinusIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $incX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->greater($N,$alpha,$XX,$offX,$incX);
    }

    public function testGreaterIllegalBufferX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
        $math->greater($N,$alpha,$XX,$offX,$incX);
    }

    public function testGreaterOverflowBufferXwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $XX = $mo->array([1,2])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->greater($N,$alpha,$XX,$offX,$incX);
    }

    public function testGreaterOverflowBufferXwithOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $offX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->greater($N,$alpha,$XX,$offX,$incX);
    }

    public function testGreaterOverflowBufferXwithIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->greater($N,$alpha,$XX,$offX,$incX);
    }

    public function testLessNormal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $math->less($N,$alpha,$XX,$offX,$incX);
        $this->assertEquals([1,0,0],$X->toArray());
    }

    public function testLessMinusN()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $N = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->less($N,$alpha,$XX,$offX,$incX);
    }

    public function testLessMinusOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $offX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->less($N,$alpha,$XX,$offX,$incX);
    }

    public function testLessMinusIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $incX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->less($N,$alpha,$XX,$offX,$incX);
    }

    public function testLessIllegalBufferX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
        $math->less($N,$alpha,$XX,$offX,$incX);
    }

    public function testLessOverflowBufferXwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $XX = $mo->array([1,2])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->less($N,$alpha,$XX,$offX,$incX);
    }

    public function testLessOverflowBufferXwithOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $offX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->less($N,$alpha,$XX,$offX,$incX);
    }

    public function testLessOverflowBufferXwithIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->less($N,$alpha,$XX,$offX,$incX);
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
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
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
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
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
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
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
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
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
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
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
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
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
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
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

        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Invalid value in sqrt.');
        $math->sqrt($N,$XX,$offX,$incX);
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
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
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

        $X = $mo->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,0,0);

        // X := 1 / (a * sqrt(X) + b)
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Zero divide.');
        $math->rsqrt($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testrsqrtInvalidSqrt()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([3,2,-1]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,1,1);

        // X := 1 / (a * sqrt(X) + b)
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Invalid value in sqrt.');
        $math->rsqrt($N,$alpha,$XX,$offX,$incX,$beta);
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
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
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
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $math->pow($N,$alpha,$XX,$offX,$incX);
        $this->assertEquals([1,4,9],$X->toArray());
    }

    public function testpowMinusN()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $N = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->pow($N,$alpha,$XX,$offX,$incX);
    }

    public function testpowMinusOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $offX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->pow($N,$alpha,$XX,$offX,$incX);
    }

    public function testpowMinusIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $incX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->pow($N,$alpha,$XX,$offX,$incX);
    }

    public function testpowIllegalBufferX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
        $math->pow($N,$alpha,$XX,$offX,$incX);
    }

    public function testpowOverflowBufferXwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $XX = $mo->array([1,2])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->pow($N,$alpha,$XX,$offX,$incX);
    }

    public function testpowOverflowBufferXwithOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $offX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->pow($N,$alpha,$XX,$offX,$incX);
    }

    public function testpowOverflowBufferXwithIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_maximum(2,$X);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->pow($N,$alpha,$XX,$offX,$incX);
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
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
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

        $X = $mo->array([1,2,0]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Invalid value in log.');
        $math->log($N,$XX,$offX,$incX);
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
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
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
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
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

    ##
    public function testSelectAxis0Normal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,2],NDArray::int32);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
        $this->assertEquals([[1,2,3],[7,8,9]],$Y->toArray());

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float64);
        $X = $mo->array([0,2],NDArray::int64);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::float64);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
        $this->assertEquals([[1,2,3],[7,8,9]],$Y->toArray());

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::int64);
        $X = $mo->array([0,2],NDArray::int64);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::int64);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
        $this->assertEquals([[1,2,3],[7,8,9]],$Y->toArray());

        $A = $mo->array([1,2,3,4],NDArray::float32);
        $X = $mo->array([0,2],NDArray::int32);
        $Y = $mo->array([0,0],NDArray::float32);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
        $this->assertEquals([1,3],$Y->toArray());

        $A = $mo->array([1,2,3,4],NDArray::float64);
        $X = $mo->array([0,2],NDArray::int64);
        $Y = $mo->array([0,0],NDArray::float64);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
        $this->assertEquals([1,3],$Y->toArray());

        $A = $mo->array([1,2,3,4],NDArray::int64);
        $X = $mo->array([0,2],NDArray::int64);
        $Y = $mo->array([0,0],NDArray::int64);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
        $this->assertEquals([1,3],$Y->toArray());
    }

    public function testSelectAxis0LabelNumberOutOfBounds1()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,4],NDArray::int32);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Label number is out of bounds.');
        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testSelectAxis0LabelNumberOutOfBounds2()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,-1],NDArray::int32);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Label number is out of bounds.');
        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testselectAxis0MinusM()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $M = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testselectAxis0MinusN()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $N = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testselectAxis0MinusK()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $K = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument k must be greater than 0.');
        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testselectAxis0MinusOffsetA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $offA = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetA must be greater than equals 0.');
        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testselectAxis0MinusLdA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $ldA = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument ldA must be greater than 0.');
        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testselectAxis0IllegalBufferA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $AA = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testselectAxis0OverflowBufferAwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $AA = $mo->array([1,2,3,4,5,6,7,8,9,10,11])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testselectAxis0OverflowBufferAwithOffset()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $offA = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testselectAxis0OverflowBufferAwithLdA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $ldA = 4;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testselectAxis0MinusOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $offX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testselectAxis0MinusIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $incX = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testselectAxis0IllegalBufferX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testselectAxis0OverflowBufferXwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $XX = $mo->array([1])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testselectAxis0OverflowBufferXwithOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $offX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testselectAxis0OverflowBufferXwithIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testselectAxis0MinusOffsetY()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $offY = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetY must be greater than equals 0.');
        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testselectAxis0MinusLdY()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $ldY = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument ldY must be greater than 0.');
        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testselectAxis0IllegalBufferY()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $YY = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testselectAxis0OverflowBufferYwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $YY = $mo->array([0])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferY');
        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testselectAxis0OverflowBufferXwithOffsetY()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $offY = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferY');
        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testselectAxis0OverflowBufferXwithLdY()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $mo->array([0,1],NDArray::int32);
        $Y = $mo->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY] =
            $this->translate_selectAxis0($A,$X,$Y);

        $ldY = 4;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferY');
        $math->selectAxis0($M,$N,$K,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

######################################################################
    public function testSelectAxis1Normal()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,2]);
        $Y = $mo->array([0,0]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_selectAxis1($A,$X,$Y);

        $math->selectAxis1($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([2,6],$Y->toArray());
    }

    public function testSelectAxis1LabelNumberOutOfBounds1()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,3]);
        $Y = $mo->array([0,0]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_selectAxis1($A,$X,$Y);

        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Label number is out of bounds.');
        $math->selectAxis1($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testSelectAxis1LabelNumberOutOfBounds2()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $Y = $mo->array([0,0]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_selectAxis1($A,$X,$Y);

        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Label number is out of bounds.');
        $math->selectAxis1($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testselectAxis1MinusM()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $Y = $mo->array([0,0]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_selectAxis1($A,$X,$Y);

        $M = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->selectAxis1($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testselectAxis1MinusN()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $Y = $mo->array([0,0]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_selectAxis1($A,$X,$Y);

        $N = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->selectAxis1($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testselectAxis1MinusOffsetA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $Y = $mo->array([0,0]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_selectAxis1($A,$X,$Y);

        $offA = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetA must be greater than equals 0.');
        $math->selectAxis1($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testselectAxis1MinusLdA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $Y = $mo->array([0,0]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_selectAxis1($A,$X,$Y);

        $ldA = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument ldA must be greater than 0.');
        $math->selectAxis1($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testselectAxis1IllegalBufferA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $Y = $mo->array([0,0]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_selectAxis1($A,$X,$Y);

        $AA = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
        $math->selectAxis1($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testselectAxis1OverflowBufferAwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $Y = $mo->array([0,0]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_selectAxis1($A,$X,$Y);

        $AA = $mo->array([1,2,3,4,5])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->selectAxis1($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testselectAxis1OverflowBufferAwithOffset()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $Y = $mo->array([0,0]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_selectAxis1($A,$X,$Y);

        $offA = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->selectAxis1($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testselectAxis1OverflowBufferAwithLdA()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $Y = $mo->array([0,0]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_selectAxis1($A,$X,$Y);

        $ldA = 4;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->selectAxis1($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testselectAxis1MinusOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $Y = $mo->array([0,0]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_selectAxis1($A,$X,$Y);

        $offX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->selectAxis1($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testselectAxis1MinusIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $Y = $mo->array([0,0]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_selectAxis1($A,$X,$Y);

        $incX = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->selectAxis1($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testselectAxis1IllegalBufferX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $Y = $mo->array([0,0]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_selectAxis1($A,$X,$Y);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
        $math->selectAxis1($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testselectAxis1OverflowBufferXwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $Y = $mo->array([0,0]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_selectAxis1($A,$X,$Y);

        $XX = $mo->array([1])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->selectAxis1($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testselectAxis1OverflowBufferXwithOffsetX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $Y = $mo->array([0,0]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_selectAxis1($A,$X,$Y);

        $offX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->selectAxis1($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testselectAxis1OverflowBufferXwithIncX()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $Y = $mo->array([0,0]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_selectAxis1($A,$X,$Y);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->selectAxis1($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testselectAxis1MinusOffsetY()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $Y = $mo->array([0,0]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_selectAxis1($A,$X,$Y);

        $offY = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetY must be greater than equals 0.');
        $math->selectAxis1($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testselectAxis1MinusIncY()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $Y = $mo->array([0,0]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_selectAxis1($A,$X,$Y);

        $incY = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incY must be greater than 0.');
        $math->selectAxis1($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testselectAxis1IllegalBufferY()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $Y = $mo->array([0,0]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_selectAxis1($A,$X,$Y);

        $YY = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
        $math->selectAxis1($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testselectAxis1OverflowBufferYwithSize()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $Y = $mo->array([0,0]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_selectAxis1($A,$X,$Y);

        $YY = $mo->array([0])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferY');
        $math->selectAxis1($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testselectAxis1OverflowBufferXwithOffsetY()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $Y = $mo->array([0,0]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_selectAxis1($A,$X,$Y);

        $offY = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferY');
        $math->selectAxis1($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testselectAxis1OverflowBufferXwithIncY()
    {
        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([1,-1]);
        $Y = $mo->array([0,0]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_selectAxis1($A,$X,$Y);

        $incY = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferY');
        $math->selectAxis1($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX,$YY,$offY,$incY);
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
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
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
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
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
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $math->reduceSum($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
        $this->assertEquals([6,15],$X->toArray());
    }

    public function testreduceSumBroadcastTranspose()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $A = $mo->array([[1,2,3],[4,5,6]]);
        $X = $mo->array([0,0,0]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_reduceSum($A,$axis=0,$X);

        $math->reduceSum($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
        $this->assertEquals([5,7,9],$X->toArray());
    }

    public function testreduceSumMinusM()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,0]);
        $A = $mo->array([[1,2,3],[4,5,6]]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $M = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->reduceSum($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testreduceSumMinusN()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,0]);
        $A = $mo->array([[1,2,3],[4,5,6]]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $N = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->reduceSum($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testreduceSumMinusOffsetX()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,0]);
        $A = $mo->array([[1,2,3],[4,5,6]]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $offX = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than equals 0.');
        $math->reduceSum($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testreduceSumMinusIncX()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,0]);
        $A = $mo->array([[1,2,3],[4,5,6]]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $incX = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->reduceSum($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testreduceSumIllegalBufferX()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,0]);
        $A = $mo->array([[1,2,3],[4,5,6]]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
        $math->reduceSum($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testreduceSumOverflowBufferXwithSize()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,0]);
        $A = $mo->array([[1,2,3],[4,5,6]]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $XX = $mo->array([1])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->reduceSum($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testreduceSumOverflowBufferXwithOffsetX()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,0]);
        $A = $mo->array([[1,2,3],[4,5,6]]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $offX = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->reduceSum($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testreduceSumOverflowBufferXwithIncX()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,0]);
        $A = $mo->array([[1,2,3],[4,5,6]]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $incX = 2;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->reduceSum($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testreduceSumMinusOffsetA()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,0]);
        $A = $mo->array([[1,2,3],[4,5,6]]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $offA = -1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument offsetA must be greater than equals 0.');
        $math->reduceSum($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testreduceSumMinusIncA()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,0]);
        $A = $mo->array([[1,2,3],[4,5,6]]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $ldA = 0;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Argument ldA must be greater than 0.');
        $math->reduceSum($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testreduceSumIllegalBufferA()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,0]);
        $A = $mo->array([[1,2,3],[4,5,6]]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $AA = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must be an instance of Rindow\OpenBLAS\Buffer');
        $math->reduceSum($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testreduceSumOverflowBufferAwithSize()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,0]);
        $A = $mo->array([[1,2,3],[4,5,6]]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $AA = $mo->array([1,2,3,4,5])->buffer();
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->reduceSum($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testreduceSumOverflowBufferXwithOffsetA()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,0]);
        $A = $mo->array([[1,2,3],[4,5,6]]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $offA = 1;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->reduceSum($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testreduceSumOverflowBufferXwithLdA()
    {
        if($this->checkSkip('reduceSum')){return;}

        $mo = new MatrixOperator();
        $math = $this->getMath($mo);

        $X = $mo->array([0,0]);
        $A = $mo->array([[10,100,1000],[10,100,1000]]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $ldA = 4;
        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->reduceSum($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
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
}
