<?php
namespace RindowTest\OpenBlas\MathTest;
if(!class_exists('RindowTest\OpenBLAS\Utils')) {
    include_once __DIR__.'/Utils.php';
}
use RindowTest\OpenBLAS\Utils;
use function RindowTest\OpenBLAS\R;

use PHPUnit\Framework\TestCase;
use Interop\Polite\Math\Matrix\NDArray;
use Interop\Polite\Math\Matrix\BLAS;
use Rindow\Math\Matrix\MatrixOperator;
use Rindow\OpenBLAS\Math as Math;
use Rindow\OpenBLAS\BLAS as OpenBLAS;
use InvalidArgumentException;
use RuntimeException;
use TypeError;

/**
 * @requires extension rindow_openblas
 */
class MathTest extends TestCase
{
    use Utils;

    public function getMath()
    {
        $math = new Math();
        return $math;
    }

    public function getBlas()
    {
        $blas = new OpenBLAS();
        return $blas;
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
        NDArray $X
        ) : array
    {
        [$m,$n] = $A->shape();
        $AA = $A->buffer();
        $offA = $A->offset();
        $XX = $X->buffer();
        $offX = $X->offset();

        return [$m,$n,$AA,$offA,$n,$XX,$offX,1];
    }

    public function translate_nan2num(
        NDArray $X,
        float $alpha
        ) : array
    {
        $n = $X->size();
        $XX = $X->buffer();
        $offX = $X->offset();

        return [
            $n,
            $XX,$offX,1,
            $alpha
        ];
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

    public function translate_duplicate(
        NDArray $X, int $n=null, bool $trans=null,NDArray $A=null) : array
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

    public function translate_pow(
       NDArray $A,
       NDArray $alpha,
       bool $trans=null
       ) : array
    {
        if($trans===null) {
            $trans = false;
        }
        $shapeA = $A->shape();
        if(is_numeric($alpha)) {
            $alpha = $this->array($alpha,$A->dtype());
        }
        $shapeX = $alpha->shape();
        if(count($shapeX)==0) {
            $trans = false;
            $shapeA = [(int)array_product($shapeA),1];
            $shapeX = [1];
        }

        if($trans) {
            $shapeA = array_reverse($shapeA);
        }
        while(true) {
            $xd = array_pop($shapeX);
            if($xd===null)
                break;
            $ad = array_pop($shapeA);
            if($xd!==$ad) {
                $shapeA = $trans ? array_reverse($A->shape()) : $A->shape();
                throw new InvalidArgumentException('Unmatch dimension size for broadcast.: '.
                    '['.implode(',',$X->shape()).'] => ['.implode(',',$shapeA).']');
            }
        }
        $n = $alpha->size();
        $XX = $alpha->buffer();
        $offX = $alpha->offset();
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
            $AA,$offA,$n,
            $XX,$offX,1,
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

    public function translate_fill(
        $value,
        NDArray $X
        )
    {
        if(is_scalar($value)) {
            if(is_string($value)) {
                $value = ord($value);
            }
            $V = $this->alloc([1],$X->dtype());
            $V[0] = $value;
        } elseif($value instanceof NDArray) {
            if($value->size()!=1) {
                throw new InvalidArgumentException('Value must be scalar');
            }
            $V = $value;
        } else {
            throw new InvalidArgumentException('Invalid data type');
        }
        $n = $X->size();
        $VV = $V->buffer();
        $offV = $V->offset();
        $XX = $X->buffer();
        $offX = $X->offset();
        return [
            $n,
            $VV, $offV,
            $XX, $offX,1
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
        if($A->ndim()==1) {
            $individual = false;
        } elseif($A->ndim()==2) {
            $individual = true;
        } else {
            throw new InvalidArgumentException('A must be 1D or 2D NDArray.');
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
        if($individual) {
            [$m,$n] = $A->shape();
            if($m!=$X->size()) {
                $shapeError = '('.implode(',',$A->shape()).'),('.implode(',',$X->shape()).')';
                throw new InvalidArgumentException("Unmatch shape of dimension A,X: ".$shapeError);
            }
            $ldA = $n;
        } else {
            $m = $X->size();
            $n = $A->size();
            $ldA = 0;
        }
        $AA = $A->buffer();
        $offA = $A->offset();
        $XX = $X->buffer();
        $offX = $X->offset();
        $YY = $Y->buffer();
        $offY = $Y->offset();
 
        return [
            $m,
            $n,
            $AA,$offA,$ldA,
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

    public function translate_transpose(
        NDArray $A,
        $perm,
        NDArray $B
        ) : array
    {
        $AA = $A->buffer();
        $BB = $B->buffer();
        $offsetA = $A->offset();
        $offsetB = $B->offset();
        $sourceShape = $this->array($A->shape(),NDArray::int32)->buffer();
        if(is_array($perm)) {
            $perm = $this->array($perm,NDArray::int32);
        }
        $permBuf = $perm->buffer();

        return [
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB,
        ];
    }

    public function translate_bandpart(
        NDArray $A,
        int $lower,
        int $upper
    ) : array
    {
        if($A->ndim()<2) {
            throw new InvalidArgumentException('input array must be 2D or upper.');
        }
        $shape = $A->shape();
        $k = array_pop($shape);
        $n = array_pop($shape);
        $m = (int)array_product($shape);
        $buffer = $A->buffer();
        $offset = $A->offset();
        return [
            $m,$n,$k,
            $buffer,$offset,
            $lower,
            $upper,
        ];
    }

    public static function providerDtypesFloats()
    {
        return [
            'float32' => [[
                'dtype' => NDArray::float32,
            ]],
            'float64' => [[
                'dtype' => NDArray::float64,
            ]],
        ];
    }

    public static function providerDtypesFloatsInt3264()
    {
        return [
            'float32' => [[
                'dtype' => NDArray::float32,
            ]],
            'float64' => [[
                'dtype' => NDArray::float64,
            ]],
            'int32' => [[
                'dtype' => NDArray::int32,
            ]],
            'int64' => [[
                'dtype' => NDArray::int64,
            ]],
        ];
    }

    public function testSumNormal()
    {
        $math = $this->getMath();
 
        $X = $this->array([100,-10,-1000],NDArray::float32);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
        $min = $math->sum($N,$XX,$offX,$incX);
        $this->assertEquals(-910,$min);
 
        $X = $this->array([100,-10,-1000],NDArray::float64);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
        $min = $math->sum($N,$XX,$offX,$incX);
        $this->assertEquals(-910,$min);
 
        $X = $this->array([-100,-100,-120],NDArray::int8);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
        $min = $math->sum($N,$XX,$offX,$incX);
        $this->assertEquals(-320,$min);
 
        $X = $this->array([-1,-2,-3],NDArray::uint8);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
        $min = $math->sum($N,$XX,$offX,$incX);
        $this->assertEquals(256*3-1-2-3,$min);
 
        $X = $this->array([-100,-100,-120],NDArray::int16);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
        $min = $math->sum($N,$XX,$offX,$incX);
        $this->assertEquals(-320,$min);
 
        $X = $this->array([-1,-2,-3],NDArray::uint16);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
        $min = $math->sum($N,$XX,$offX,$incX);
        $this->assertEquals(65536*3-1-2-3,$min);
 
        $X = $this->array([-100,-100,-120],NDArray::int32);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
        $min = $math->sum($N,$XX,$offX,$incX);
        $this->assertEquals(-320,$min);
 
        $X = $this->array([-1,-2,-3],NDArray::uint32);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
        $min = $math->sum($N,$XX,$offX,$incX);
        $this->assertEquals((2**32)*3-1-2-3,$min);
 
        $X = $this->array([-100,-100,-120],NDArray::int64);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
        $min = $math->sum($N,$XX,$offX,$incX);
        $this->assertEquals(-320,$min);
 
        $X = $this->array([-1,-2,-3],NDArray::uint64);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
        $min = $math->sum($N,$XX,$offX,$incX);
        //$this->assertEquals((2**64)*3-1-2-3,$min);
        $this->assertEquals(-6,$min);
 
        $X = $this->array([true,false,true],NDArray::bool);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
        $min = $math->sum($N,$XX,$offX,$incX);
        $this->assertEquals(2,$min);
    }
 
    public function testSumMinusN()
    {
        $math = $this->getMath();
 
        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
 
        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $min = $math->sum($N,$XX,$offX,$incX);
    }
 
    public function testSumMinusOffsetX()
    {
        $math = $this->getMath();
 
        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
 
        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $min = $math->sum($N,$XX,$offX,$incX);
    }
 
    public function testSumMinusIncX()
    {
        $math = $this->getMath();
 
        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
 
        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $min = $math->sum($N,$XX,$offX,$incX);
    }
 
    public function testSumIllegalBufferX()
    {
        $math = $this->getMath();
 
        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
 
        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $min = $math->sum($N,$XX,$offX,$incX);
    }
 
    public function testSumOverflowBufferXwithSize()
    {
        $math = $this->getMath();
 
        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
 
        $XX = $this->array([100,-10])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $min = $math->sum($N,$XX,$offX,$incX);
    }
 
    public function testSumOverflowBufferXwithOffsetX()
    {
        $math = $this->getMath();
 
        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
 
        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $min = $math->sum($N,$XX,$offX,$incX);
    }
 
    public function testSumOverflowBufferXwithIncX()
    {
        $math = $this->getMath();
 
        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
 
        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $min = $math->sum($N,$XX,$offX,$incX);
    }
 
    /**
    * @dataProvider providerDtypesFloats
    */
    public function testMinNormal($params)
    {
        extract($params);
        $math = $this->getMath();
 
        $X = $this->array([100,-10,1],$dtype);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
 
        $min = $math->imin($N,$XX,$offX,$incX);
        $this->assertEquals(1,$min);
    }
 
    public function testMinMinusN()
    {
        $math = $this->getMath();
 
        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
 
        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $min = $math->imin($N,$XX,$offX,$incX);
    }
 
    public function testMinMinusOffsetX()
    {
        $math = $this->getMath();
 
        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
 
        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $min = $math->imin($N,$XX,$offX,$incX);
    }
 
    public function testMinMinusIncX()
    {
        $math = $this->getMath();
 
        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
 
        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $min = $math->imin($N,$XX,$offX,$incX);
    }
 
    public function testMinIllegalBufferX()
    {
        $math = $this->getMath();
 
        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
 
        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $min = $math->imin($N,$XX,$offX,$incX);
    }
 
    public function testMinOverflowBufferXwithSize()
    {
        $math = $this->getMath();
 
        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
 
        $XX = $this->array([100,-10])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $min = $math->imin($N,$XX,$offX,$incX);
    }
 
    public function testMinOverflowBufferXwithOffsetX()
    {
        $math = $this->getMath();
 
        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
 
        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $min = $math->imin($N,$XX,$offX,$incX);
    }
 
    public function testMinOverflowBufferXwithIncX()
    {
        $math = $this->getMath();
 
        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
 
        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $min = $math->imin($N,$XX,$offX,$incX);
    }
 
    /**
    * @dataProvider providerDtypesFloats
    */
    public function testMaxNormal($params)
    {
        extract($params);
        $math = $this->getMath();
 
        $X = $this->array([100,-10,-1000],$dtype);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
 
        $min = $math->imax($N,$XX,$offX,$incX);
        $this->assertEquals(0,$min);
    }
 
    public function testMaxMinusN()
    {
        $math = $this->getMath();
 
        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
 
        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $min = $math->imax($N,$XX,$offX,$incX);
    }
 
    public function testMaxMinusOffsetX()
    {
        $math = $this->getMath();
 
        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
 
        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $min = $math->imax($N,$XX,$offX,$incX);
    }
 
    public function testMaxMinusIncX()
    {
        $math = $this->getMath();
 
        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
 
        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $min = $math->imax($N,$XX,$offX,$incX);
    }
 
    public function testMaxIllegalBufferX()
    {
        $math = $this->getMath();
 
        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
 
        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $min = $math->imax($N,$XX,$offX,$incX);
    }
 
    public function testMaxOverflowBufferXwithSize()
    {
        $math = $this->getMath();
 
        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
 
        $XX = $this->array([100,-10])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $min = $math->imax($N,$XX,$offX,$incX);
    }
 
    public function testMaxOverflowBufferXwithOffsetX()
    {
        $math = $this->getMath();
 
        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
 
        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $min = $math->imax($N,$XX,$offX,$incX);
    }
 
    public function testMaxOverflowBufferXwithIncX()
    {
        $math = $this->getMath();
 
        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);
 
        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $min = $math->imax($N,$XX,$offX,$incX);
    }
 
 
    /**
    * @dataProvider providerDtypesFloats
    */
    public function testIncrementNormal($params)
    {
        extract($params);
        $math = $this->getMath();
 
        $X = $this->array([1,2,3],$dtype);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,10,2);

        $math->increment($N,$alpha,$XX,$offX,$incX,$beta);
        $this->assertEquals([12,14,16],$X->toArray());
    }

    public function testIncrementInvalidArgments()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
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
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,10,2);

        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->increment($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testIncrementMinusOffsetX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,10,2);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $math->increment($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testIncrementMinusIncX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,10,2);

        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->increment($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testIncrementIllegalBufferX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,10,2);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->increment($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testIncrementOverflowBufferXwithSize()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,10,2);

        $XX = $this->array([1,2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->increment($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testIncrementOverflowBufferXwithOffsetX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,10,2);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->increment($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testIncrementOverflowBufferXwithIncX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,10,2);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->increment($N,$alpha,$XX,$offX,$incX,$beta);
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testReciprocalNormal($params)
    {
        extract($params);
        $math = $this->getMath();

        $X = $this->array([3,2,0],$dtype);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $math->reciprocal($N,$alpha,$XX,$offX,$incX,$beta);
        // X := 1 / ( alpha * X + beta )
        //    = [1/1, 1/2, 1/4]
        $this->assertEquals([1,0.5,0.25],$X->toArray());
    }

    public function testReciprocalZeroDivide()
    {
        $math = $this->getMath();

        $X = $this->array([4,2,0]);
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
        $math = $this->getMath();

        $X = $this->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->reciprocal($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testReciprocalMinusOffsetX()
    {
        $math = $this->getMath();

        $X = $this->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $math->reciprocal($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testReciprocalMinusIncX()
    {
        $math = $this->getMath();

        $X = $this->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->reciprocal($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testReciprocalIllegalBufferX()
    {
        $math = $this->getMath();

        $X = $this->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->reciprocal($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testReciprocalOverflowBufferXwithSize()
    {
        $math = $this->getMath();

        $X = $this->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $XX = $this->array([3,2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->reciprocal($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testReciprocalOverflowBufferXwithOffsetX()
    {
        $math = $this->getMath();

        $X = $this->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->reciprocal($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testReciprocalOverflowBufferXwithIncX()
    {
        $math = $this->getMath();

        $X = $this->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->reciprocal($N,$alpha,$XX,$offX,$incX,$beta);
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testMaximumNormal($params)
    {
        extract($params);
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]],$dtype);
        $X = $this->array([2,3],$dtype);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
        $this->assertEquals([[2,3],[2,3],[3,4]],$A->toArray());
    }

    public function testMaximumMinusM()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $M = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumMinusN()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumMinusOffsetA()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offA = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetA must be greater than or equals 0.');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumMinusLdA()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $ldA = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument ldA must be greater than 0.');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumIllegalBufferA()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $AA = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumOverflowBufferAwithSize()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $AA = $this->array([1,2,2,3,3])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumOverflowBufferXwithOffsetA()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offA = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumOverflowBufferXwithLdA()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $ldA = 3;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumMinusOffsetX()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumMinusIncX()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $incX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumIllegalBufferX()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumOverflowBufferXwithSize()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $XX = $this->array([1])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumOverflowBufferXwithOffsetX()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMaximumOverflowBufferXwithIncX()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->maximum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testMinimumNormal($params)
    {
        extract($params);
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]],$dtype);
        $X = $this->array([2,3],$dtype);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
        $this->assertEquals([[1,2],[2,3],[2,3]],$A->toArray());
    }

    public function testMinimumMinusM()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $M = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumMinusN()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumMinusOffsetA()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offA = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetA must be greater than or equals 0.');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumMinusLdA()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $ldA = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument ldA must be greater than 0.');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumIllegalBufferA()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $AA = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumOverflowBufferAwithSize()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $AA = $this->array([1,2,2,3,3])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumOverflowBufferXwithOffsetA()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offA = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumOverflowBufferXwithLdA()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $ldA = 3;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumMinusOffsetX()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumMinusIncX()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $incX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumIllegalBufferX()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumOverflowBufferXwithSize()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $XX = $this->array([1])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumOverflowBufferXwithOffsetX()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testMinimumOverflowBufferXwithIncX()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->minimum($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }
    
    /**
    * @dataProvider providerDtypesFloats
    */
    public function testGreaterNormal($params)
    {
        extract($params);
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]],$dtype);
        $X = $this->array([2,3],$dtype);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
        $this->assertEquals([[0,0],[0,0],[1,1]],$A->toArray());
    }

    public function testGreaterMinusM()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $M = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterMinusN()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterMinusOffsetA()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offA = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetA must be greater than or equals 0.');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterMinusLdA()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $ldA = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument ldA must be greater than 0.');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterIllegalBufferA()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $AA = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterOverflowBufferAwithSize()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $AA = $this->array([1,2,2,3,3])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterOverflowBufferXwithOffsetA()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offA = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterOverflowBufferXwithLdA()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $ldA = 3;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterMinusOffsetX()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterMinusIncX()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $incX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterIllegalBufferX()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterOverflowBufferXwithSize()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $XX = $this->array([1])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterOverflowBufferXwithOffsetX()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testGreaterOverflowBufferXwithIncX()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->greater($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testGreaterEqualNormal($params)
    {
        extract($params);
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]],$dtype);
        $X = $this->array([2,3],$dtype);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $math->greaterEqual($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
        $this->assertEquals([[0,0],[1,1],[1,1]],$A->toArray());
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testLessNormal($params)
    {
        extract($params);
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]],$dtype);
        $X = $this->array([2,3],$dtype);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
        $this->assertEquals([[1,1],[0,0],[0,0]],$A->toArray());
    }

    public function testLessMinusM()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $M = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessMinusN()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessMinusOffsetA()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offA = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetA must be greater than or equals 0.');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessMinusLdA()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $ldA = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument ldA must be greater than 0.');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessIllegalBufferA()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $AA = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessOverflowBufferAwithSize()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $AA = $this->array([1,2,2,3,3])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessOverflowBufferXwithOffsetA()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offA = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessOverflowBufferXwithLdA()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $ldA = 3;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessMinusOffsetX()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessMinusIncX()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $incX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessIllegalBufferX()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessOverflowBufferXwithSize()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $XX = $this->array([1])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessOverflowBufferXwithOffsetX()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    public function testLessOverflowBufferXwithIncX()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]]);
        $X = $this->array([2,3]);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->less($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testLessEqualNormal($params)
    {
        extract($params);
        $math = $this->getMath();

        $A = $this->array([[1,2],[2,3],[3,4]],$dtype);
        $X = $this->array([2,3],$dtype);
        [$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_maximum($A,$X);

        $math->lessEqual($M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);
        $this->assertEquals([[1,1],[1,1],[0,0]],$A->toArray());
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testMultiplySameSizeNormal($params)
    {
        extract($params);
        if($this->checkSkip('multiply')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3],$dtype);
        $A = $this->array([10,100,1000],$dtype);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
        $this->assertEquals([10,200,3000],$A->toArray());
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testMultiplyBroadcastNormal($params)
    {
        extract($params);
        if($this->checkSkip('multiply')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3],$dtype);
        $A = $this->array([[10,100,1000],[-1,-1,-1]],$dtype);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
        $this->assertEquals([[10,200,3000],[-1,-2,-3]],$A->toArray());
    }

    public function testMultiplyBroadcastTranspose()
    {
        if($this->checkSkip('multiply')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([[10,100],[1000,10000],[-1,-1]]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A,true);

        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
        $this->assertEquals([[10,100],[2000,20000],[-3,-3]],$A->toArray());
    }

    public function testMultiplyMinusM()
    {
        if($this->checkSkip('multiply')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $M = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testMultiplyMinusN()
    {
        if($this->checkSkip('multiply')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testMultiplyMinusOffsetX()
    {
        if($this->checkSkip('multiply')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testMultiplyMinusIncX()
    {
        if($this->checkSkip('multiply')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testMultiplyIllegalBufferX()
    {
        if($this->checkSkip('multiply')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
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

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $XX = $this->array([1,2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testMultiplyOverflowBufferXwithOffsetX()
    {
        if($this->checkSkip('multiply')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testMultiplyOverflowBufferXwithIncX()
    {
        if($this->checkSkip('multiply')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testMultiplyMinusOffsetA()
    {
        if($this->checkSkip('multiply')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $offA = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetA must be greater than or equals 0.');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testMultiplyMinusIncA()
    {
        if($this->checkSkip('multiply')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $ldA = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument ldA must be greater than 0.');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testMultiplyIllegalBufferA()
    {
        if($this->checkSkip('multiply')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
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

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $AA = $this->array([1,2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testMultiplyOverflowBufferXwithOffsetA()
    {
        if($this->checkSkip('multiply')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $offA = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testMultiplyOverflowBufferXwithLdA()
    {
        if($this->checkSkip('multiply')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([[10,100,1000],[10,100,1000]]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_multiply($X,$A);

        $ldA = 4;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->multiply($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testaddSameSizeNormal($params)
    {
        extract($params);
        if($this->checkSkip('add')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3],$dtype);
        $A = $this->array([10,100,1000],$dtype);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A,-1);

        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
        $this->assertEquals([9,98,997],$A->toArray());
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testaddBroadcastNormal($params)
    {
        extract($params);
        if($this->checkSkip('add')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3],$dtype);
        $A = $this->array([[10,100,1000],[-1,-1,-1]],$dtype);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
        $this->assertEquals([[11,102,1003],[0,1,2]],$A->toArray());
    }

    public function testaddBroadcastTranspose()
    {
        if($this->checkSkip('add')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([[10,100],[1000,10000],[-1,-1]]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A,null,true);

        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
        $this->assertEquals([[11,101],[1002,10002],[2,2]],$A->toArray());
    }

    public function testaddMinusM()
    {
        if($this->checkSkip('add')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $M = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testaddMinusN()
    {
        if($this->checkSkip('add')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testaddMinusOffsetX()
    {
        if($this->checkSkip('add')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testaddMinusIncX()
    {
        if($this->checkSkip('add')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testaddIllegalBufferX()
    {
        if($this->checkSkip('add')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
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

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $XX = $this->array([1,2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testaddOverflowBufferXwithOffsetX()
    {
        if($this->checkSkip('add')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testaddOverflowBufferXwithIncX()
    {
        if($this->checkSkip('add')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testaddMinusOffsetA()
    {
        if($this->checkSkip('add')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $offA = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetA must be greater than or equals 0.');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testaddMinusIncA()
    {
        if($this->checkSkip('add')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $ldA = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument ldA must be greater than 0.');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testaddIllegalBufferA()
    {
        if($this->checkSkip('add')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
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

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $AA = $this->array([1,2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testaddOverflowBufferXwithOffsetA()
    {
        if($this->checkSkip('add')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $offA = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testaddOverflowBufferXwithLdA()
    {
        if($this->checkSkip('add')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([[10,100,1000],[10,100,1000]]);
        [$trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_add($X,$A);

        $ldA = 4;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->add($trans,$M,$N,$alpha,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testDuplicateSameSizeNormal($params)
    {
        extract($params);
        if($this->checkSkip('duplicate')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3],$dtype);
        $A = $this->array([10,100,1000],$dtype);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
        $this->assertEquals([1,2,3],$A->toArray());
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testDuplicateBroadcastNormal($params)
    {
        extract($params);
        if($this->checkSkip('duplicate')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3],$dtype);
        $A = $this->array([[10,100,1000],[-1,-1,-1]],$dtype);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
        $this->assertEquals([[1,2,3],[1,2,3]],$A->toArray());
    }

    public function testDuplicateBroadcastTranspose()
    {
        if($this->checkSkip('duplicate')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2]);
        $A = $this->array([[10,100,1000],[-1,-1,-1]]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,true,$A);
    
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
        $this->assertEquals([[1,1,1],[2,2,2]],$A->toArray());
    }

    public function testDuplicateMinusM()
    {
        if($this->checkSkip('duplicate')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $M = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testDuplicateMinusN()
    {
        if($this->checkSkip('duplicate')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testDuplicateMinusOffsetX()
    {
        if($this->checkSkip('duplicate')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testDuplicateMinusIncX()
    {
        if($this->checkSkip('duplicate')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testDuplicateIllegalBufferX()
    {
        if($this->checkSkip('duplicate')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
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

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $XX = $this->array([2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testDuplicateOverflowBufferXwithOffsetX()
    {
        if($this->checkSkip('duplicate')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testDuplicateOverflowBufferXwithIncX()
    {
        if($this->checkSkip('duplicate')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testDuplicateMinusOffsetA()
    {
        if($this->checkSkip('duplicate')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $offA = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetA must be greater than or equals 0.');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testDuplicateMinusLdA()
    {
        if($this->checkSkip('duplicate')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $ldA = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument ldA must be greater than 0.');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testDuplicateIllegalBufferA()
    {
        if($this->checkSkip('duplicate')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
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

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $AA = $this->array([1,2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testDuplicateOverflowBufferAwithOffsetA()
    {
        if($this->checkSkip('duplicate')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $offA = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    public function testDuplicateOverflowBufferAwithLdA()
    {
        if($this->checkSkip('duplicate')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([[10,100,1000],[10,100,1000]]);
        [$trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA] =
            $this->translate_duplicate($X,null,null,$A);

        $ldA = 4;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->duplicate($trans,$M,$N,$XX,$offX,$incX,$AA,$offA,$ldA);
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testSquareNormal($params)
    {
        extract($params);
        $math = $this->getMath();

        $X = $this->array([1,2,3],$dtype);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $math->square($N,$XX,$offX,$incX);
        $this->assertEquals([1,4,9],$X->toArray());
    }

    public function testsquareMinusN()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->square($N,$XX,$offX,$incX);
    }

    public function testsquareMinusOffsetX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $math->square($N,$XX,$offX,$incX);
    }

    public function testsquareMinusIncX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->square($N,$XX,$offX,$incX);
    }

    public function testsquareIllegalBufferX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->square($N,$XX,$offX,$incX);
    }

    public function testsquareOverflowBufferXwithSize()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $XX = $this->array([1,2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->square($N,$XX,$offX,$incX);
    }

    public function testsquareOverflowBufferXwithOffsetX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->square($N,$XX,$offX,$incX);
    }

    public function testsquareOverflowBufferXwithIncX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->square($N,$XX,$offX,$incX);
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testsqrtNormal($params)
    {
        extract($params);
        $math = $this->getMath();

        $X = $this->array([0,1,4,9],$dtype);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $math->sqrt($N,$XX,$offX,$incX);
        $this->assertEquals([0,1,2,3],$X->toArray());
    }

    public function testsqrtIllegalValue()
    {
        $math = $this->getMath();

        $X = $this->array([1,4,-1]);
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
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->sqrt($N,$XX,$offX,$incX);
    }

    public function testsqrtMinusOffsetX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $math->sqrt($N,$XX,$offX,$incX);
    }

    public function testsqrtMinusIncX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->sqrt($N,$XX,$offX,$incX);
    }

    public function testssqrtIllegalBufferX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->sqrt($N,$XX,$offX,$incX);
    }

    public function testsqrtOverflowBufferXwithSize()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $XX = $this->array([1,2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->sqrt($N,$XX,$offX,$incX);
    }

    public function testsqrtOverflowBufferXwithOffsetX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->sqrt($N,$XX,$offX,$incX);
    }

    public function testsqrtOverflowBufferXwithIncX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->sqrt($N,$XX,$offX,$incX);
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testrsqrtNormal($params)
    {
        extract($params);
        $math = $this->getMath();

        $X = $this->array([1,4,16],$dtype);
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
        $math = $this->getMath();

        $X = $this->array([4,1,0]);
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
        $math = $this->getMath();

        $X = $this->array([4,1,-1]);
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
        $math = $this->getMath();

        $X = $this->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->rsqrt($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testrsqrtMinusOffsetX()
    {
        $math = $this->getMath();

        $X = $this->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $math->rsqrt($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testrsqrtMinusIncX()
    {
        $math = $this->getMath();

        $X = $this->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->rsqrt($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testrsqrtIllegalBufferX()
    {
        $math = $this->getMath();

        $X = $this->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->rsqrt($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testrsqrtOverflowBufferXwithSize()
    {
        $math = $this->getMath();

        $X = $this->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $XX = $this->array([3,2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->rsqrt($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testrsqrtOverflowBufferXwithOffsetX()
    {
        $math = $this->getMath();

        $X = $this->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->rsqrt($N,$alpha,$XX,$offX,$incX,$beta);
    }

    public function testrsqrtOverflowBufferXwithIncX()
    {
        $math = $this->getMath();

        $X = $this->array([3,2,0]);
        [$N,$alpha,$XX,$offX,$incX,$beta] =
            $this->translate_increment($X,4,-1);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->rsqrt($N,$alpha,$XX,$offX,$incX,$beta);
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testPowSameSizeNormal($params)
    {
        extract($params);
        if($this->checkSkip('pow')){return;}

        $math = $this->getMath();

        $A = $this->array([1,2,3],$dtype);
        $X = $this->array([4,3,2],$dtype);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_pow($A,$X);

        $math->pow($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);;
        $this->assertEquals([1,8,9],$A->toArray());
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testPowBroadcastNormal($params)
    {
        extract($params);
        if($this->checkSkip('pow')){return;}

        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6]],$dtype);
        $X = $this->array([4,3,2],$dtype);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_pow($A,$X);

        $math->pow($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);;
        $this->assertEquals([[1,8,9],[256,125,36]],$A->toArray());
    }

    public function testPowBroadcastTranspose()
    {
        if($this->checkSkip('pow')){return;}

        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([3,2]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_pow($A,$X,true);

        $math->pow($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);;
        $this->assertEquals([[1,8,27],[16,25,36]],$A->toArray());
    }

    public function testPowMinusM()
    {
        if($this->checkSkip('pow')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_pow($A,$X);

        $M = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->pow($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);;
    }

    public function testPowMinusN()
    {
        if($this->checkSkip('pow')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_pow($A,$X);

        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->pow($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);;
    }

    public function testPowMinusOffsetX()
    {
        if($this->checkSkip('pow')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_pow($A,$X);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $math->pow($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);;
    }

    public function testPowMinusIncX()
    {
        if($this->checkSkip('pow')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_pow($A,$X);

        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->pow($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);;
    }

    public function testPowIllegalBufferX()
    {
        if($this->checkSkip('pow')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_pow($A,$X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->pow($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);;
    }

    public function testPowOverflowBufferXwithSize()
    {
        if($this->checkSkip('pow')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_pow($A,$X);

        $XX = $this->array([1,2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->pow($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);;
    }

    public function testPowOverflowBufferXwithOffsetX()
    {
        if($this->checkSkip('pow')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_pow($A,$X);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->pow($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);;
    }

    public function testPowOverflowBufferXwithIncX()
    {
        if($this->checkSkip('pow')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_pow($A,$X);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->pow($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);;
    }

    public function testPowMinusOffsetA()
    {
        if($this->checkSkip('pow')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_pow($A,$X);

        $offA = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetA must be greater than or equals 0.');
        $math->pow($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);;
    }

    public function testPowMinusIncA()
    {
        if($this->checkSkip('pow')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_pow($A,$X);

        $ldA = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument ldA must be greater than 0.');
        $math->pow($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);;
    }

    public function testPowIllegalBufferA()
    {
        if($this->checkSkip('pow')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_pow($A,$X);

        $AA = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->pow($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);;
    }

    public function testPowOverflowBufferAwithSize()
    {
        if($this->checkSkip('pow')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_pow($A,$X);

        $AA = $this->array([1,2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->pow($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);;
    }

    public function testPowOverflowBufferXwithOffsetA()
    {
        if($this->checkSkip('pow')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([10,100,1000]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_pow($A,$X);

        $offA = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->pow($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);;
    }

    public function testPowOverflowBufferXwithLdA()
    {
        if($this->checkSkip('pow')){return;}

        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $A = $this->array([[10,100,1000],[10,100,1000]]);
        [$trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX] =
            $this->translate_pow($A,$X);

        $ldA = 4;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->pow($trans,$M,$N,$AA,$offA,$ldA,$XX,$offX,$incX);;
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testexpNormal($params)
    {
        extract($params);
        $math = $this->getMath();

        $X = $this->array([0,2,4,9],$dtype);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $math->exp($N,$XX,$offX,$incX);
        $math->log($N,$XX,$offX,$incX);

        $this->assertTrue($this->isclose($this->array([0,2,4,9],$dtype),$X));
    }

    public function testexpMinusN()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->exp($N,$XX,$offX,$incX);
    }

    public function testexpMinusOffsetX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $math->exp($N,$XX,$offX,$incX);
    }

    public function testexpMinusIncX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->exp($N,$XX,$offX,$incX);
    }

    public function testexpIllegalBufferX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->exp($N,$XX,$offX,$incX);
    }

    public function testexpOverflowBufferXwithSize()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $XX = $this->array([1,2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->exp($N,$XX,$offX,$incX);
    }

    public function testexpOverflowBufferXwithOffsetX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->exp($N,$XX,$offX,$incX);
    }

    public function testexpOverflowBufferXwithIncX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->exp($N,$XX,$offX,$incX);
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testlogNormal($params)
    {
        extract($params);
        $math = $this->getMath();

        $X = $this->array([1,2,4,9],$dtype);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $math->log($N,$XX,$offX,$incX);
        $math->exp($N,$XX,$offX,$incX);

        $this->assertTrue($this->isclose($this->array([1,2,4,9],$dtype),$X));
    }

    public function testlogInvalidValue()
    {
        $math = $this->getMath();

        $X = $this->array([1,0,-1]);
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
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->log($N,$XX,$offX,$incX);
    }

    public function testlogMinusOffsetX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $math->log($N,$XX,$offX,$incX);
    }

    public function testlogMinusIncX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->log($N,$XX,$offX,$incX);
    }

    public function testlogIllegalBufferX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->log($N,$XX,$offX,$incX);
    }

    public function testlogOverflowBufferXwithSize()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $XX = $this->array([1,2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->log($N,$XX,$offX,$incX);
    }

    public function testlogOverflowBufferXwithOffsetX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->log($N,$XX,$offX,$incX);
    }

    public function testlogOverflowBufferXwithIncX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->log($N,$XX,$offX,$incX);
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testfillNormal($params)
    {
        extract($params);
        $math = $this->getMath();

        $X = $this->array([NAN,NAN,NAN,NAN],$dtype);
        [$N, $VV, $offV, $XX, $offX, $incX] =
            $this->translate_fill($this->array(1.0,$X->dtype()),$X);
        $math->fill($N, $VV, $offV, $XX, $offX, $incX);
        $this->assertEquals([1,1,1,1],$X->toArray());
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testnan2numNormal($params)
    {
        extract($params);
        $math = $this->getMath();

        $X = $this->array([NAN,2,4,NAN],$dtype);
        [$N,$XX,$offX,$incX,$alpha] =
            $this->translate_nan2num($X,0.0);
        $math->nan2num($N,$XX,$offX,$incX,$alpha);
        $this->assertEquals([0,2,4,0],$X->toArray());

        $X = $this->array([NAN,2,4,NAN],$dtype);
        [$N,$XX,$offX,$incX,$alpha] =
            $this->translate_nan2num($X,1.0);
        $math->nan2num($N,$XX,$offX,$incX,$alpha);
        $this->assertEquals([1,2,4,1],$X->toArray());
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testisnanNormal($params)
    {
        extract($params);
        $math = $this->getMath();

        $X = $this->array([NAN,2,4,NAN],$dtype);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $math->isnan($N,$XX,$offX,$incX);

        $this->assertEquals([1,0,0,1],$X->toArray());
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testsearchsortedNormal($params)
    {
        extract($params);
        $math = $this->getMath();

        $A = $this->array([0,2,4],$dtype);
        $X = $this->array([-1,1,2,5],$dtype);
        $Y = $this->zeros([4],NDArray::int32);
        [$m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $math->searchsorted($m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
        $this->assertEquals([0,1,1,3],$Y->toArray());

        $math->searchsorted($m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,true,$YY,$offsetY,$incY);
        $this->assertEquals([0,1,2,3],$Y->toArray());
    }

    public function testsearchsortedIndividual()
    {
        $math = $this->getMath();

        $A = $this->array([
            [1,   3,  5,   7,   9],
            [1,   2,  3,   4,   5],
            [0, 100, 20, 300, 400]
        ]);
        $X = $this->array([0, 5, 10]);
        $Y = $this->zeros([3],NDArray::int32);
        [$m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $math->searchsorted($m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
        $this->assertEquals([0, 4, 1],$Y->toArray());

        $math->searchsorted($m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,true,$YY,$offsetY,$incY);
        $this->assertEquals([0, 5, 1],$Y->toArray());
    }

    public function testsearchsortedMinusM()
    {
        $math = $this->getMath();

        $A = $this->array([0,2,4]);
        $X = $this->array([-1,1,2,5]);
        $Y = $this->zeros([4],NDArray::int32);
        [$m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $m = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->searchsorted($m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedMinusOffsetA()
    {
        $math = $this->getMath();

        $A = $this->array([0,2,4]);
        $X = $this->array([-1,1,2,5]);
        $Y = $this->zeros([4],NDArray::int32);
        [$m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $offsetA = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetA must be greater than or equals 0.');
        $math->searchsorted($m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedMinusIncA()
    {
        $math = $this->getMath();

        $A = $this->array([0,2,4]);
        $X = $this->array([-1,1,2,5]);
        $Y = $this->zeros([4],NDArray::int32);
        [$m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $ldA = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument ldA must be greater than or equals 0.');
        $math->searchsorted($m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedIllegalBufferA()
    {
        $math = $this->getMath();

        $A = $this->array([0,2,4]);
        $X = $this->array([-1,1,2,5]);
        $Y = $this->zeros([4],NDArray::int32);
        [$m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $AA = new \stdClass();
        $this->expectException(TypeError::class);
        $this->argExpectExceptionMessage('A must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->searchsorted($m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedOverflowBufferAwithSize()
    {
        $math = $this->getMath();

        $A = $this->array([0,2,4]);
        $X = $this->array([-1,1,2,5]);
        $Y = $this->zeros([4],NDArray::int32);
        [$m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $AA = $this->array([1,2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA.');
        $math->searchsorted($m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedOverflowBufferAwithOffsetA()
    {
        $math = $this->getMath();

        $A = $this->array([0,2,4]);
        $X = $this->array([-1,1,2,5]);
        $Y = $this->zeros([4],NDArray::int32);
        [$m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $offsetA = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA.');
        $math->searchsorted($m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedOverflowBufferAwithIncA()
    {
        $math = $this->getMath();

        $A = $this->array([0,2,4]);
        $X = $this->array([-1,1,2,5]);
        $Y = $this->zeros([4],NDArray::int32);
        [$m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $ldA = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA.');
        $math->searchsorted($m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedMinusN()
    {
        $math = $this->getMath();

        $A = $this->array([0,2,4]);
        $X = $this->array([-1,1,2,5]);
        $Y = $this->zeros([4],NDArray::int32);
        [$m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $n = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->searchsorted($m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedMinusOffsetX()
    {
        $math = $this->getMath();

        $A = $this->array([0,2,4]);
        $X = $this->array([-1,1,2,5]);
        $Y = $this->zeros([4],NDArray::int32);
        [$m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $offsetX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $math->searchsorted($m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedMinusIncX()
    {
        $math = $this->getMath();

        $A = $this->array([0,2,4]);
        $X = $this->array([-1,1,2,5]);
        $Y = $this->zeros([4],NDArray::int32);
        [$m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->searchsorted($m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedIllegalBufferX()
    {
        $math = $this->getMath();

        $A = $this->array([0,2,4]);
        $X = $this->array([-1,1,2,5]);
        $Y = $this->zeros([4],NDArray::int32);
        [$m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->argExpectExceptionMessage('X must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->searchsorted($m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedOverflowBufferXwithSize()
    {
        $math = $this->getMath();

        $A = $this->array([0,2,4]);
        $X = $this->array([-1,1,2,5]);
        $Y = $this->zeros([4],NDArray::int32);
        [$m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $XX = $this->array([1,2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->searchsorted($m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedOverflowBufferXwithOffsetX()
    {
        $math = $this->getMath();

        $A = $this->array([0,2,4]);
        $X = $this->array([-1,1,2,5]);
        $Y = $this->zeros([4],NDArray::int32);
        [$m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $offsetX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->searchsorted($m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedOverflowBufferXwithIncX()
    {
        $math = $this->getMath();

        $A = $this->array([0,2,4]);
        $X = $this->array([-1,1,2,5]);
        $Y = $this->zeros([4],NDArray::int32);
        [$m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->searchsorted($m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedMinusOffsetY()
    {
        $math = $this->getMath();

        $A = $this->array([0,2,4]);
        $X = $this->array([-1,1,2,5]);
        $Y = $this->zeros([4],NDArray::int32);
        [$m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $offsetY = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetY must be greater than or equals 0.');
        $math->searchsorted($m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedMinusIncY()
    {
        $math = $this->getMath();

        $A = $this->array([0,2,4]);
        $X = $this->array([-1,1,2,5]);
        $Y = $this->zeros([4],NDArray::int32);
        [$m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $incY = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incY must be greater than 0.');
        $math->searchsorted($m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedIllegalBufferY()
    {
        $math = $this->getMath();

        $A = $this->array([0,2,4]);
        $X = $this->array([-1,1,2,5]);
        $Y = $this->zeros([4],NDArray::int32);
        [$m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $YY = new \stdClass();
        $this->expectException(TypeError::class);
        $this->argExpectExceptionMessage('Y must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->searchsorted($m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedOverflowBufferYwithSize()
    {
        $math = $this->getMath();

        $A = $this->array([0,2,4]);
        $X = $this->array([-1,1,2,5]);
        $Y = $this->zeros([4],NDArray::int32);
        [$m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $YY = $this->array([1,2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferY');
        $math->searchsorted($m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedOverflowBufferYwithOffsetY()
    {
        $math = $this->getMath();

        $A = $this->array([0,2,4]);
        $X = $this->array([-1,1,2,5]);
        $Y = $this->zeros([4],NDArray::int32);
        [$m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $offsetX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->searchsorted($m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedOverflowBufferYwithIncY()
    {
        $math = $this->getMath();

        $A = $this->array([0,2,4]);
        $X = $this->array([-1,1,2,5]);
        $Y = $this->zeros([4],NDArray::int32);
        [$m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $incY = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferY');
        $math->searchsorted($m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    public function testsearchsortedUnmatchDataType()
    {
        $math = $this->getMath();

        $A = $this->array([0,2,4],NDArray::float32);
        $X = $this->array([-1,1,2,5],NDArray::float64);
        $Y = $this->zeros([4],NDArray::int32);
        [$m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY] =
            $this->translate_searchsorted($A,$X,false,null,$Y);

        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Unmatch data type for A and X');
        $math->searchsorted($m,$n,$AA,$offsetA,$ldA,$XX,$offsetX,$incX,$right,$YY,$offsetY,$incY);
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testcumsumNormal($params)
    {
        extract($params);
        $math = $this->getMath();

        $X = $this->array([1,2,3],$dtype);
        $Y = $this->zeros([3],$dtype);
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

        $X = $this->array([1,NAN,3],$dtype);
        $Y = $this->zeros([3],$dtype);
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
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $Y = $this->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $n = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumMinusOffsetX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $Y = $this->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $offsetX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumMinusIncX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $Y = $this->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumIllegalBufferX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $Y = $this->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->argExpectExceptionMessage('X must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumOverflowBufferXwithSize()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $Y = $this->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $XX = $this->array([1,2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumOverflowBufferXwithOffsetX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $Y = $this->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $offsetX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumOverflowBufferXwithIncX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $Y = $this->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumMinusOffsetY()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $Y = $this->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $offsetY = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetY must be greater than or equals 0.');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumMinusIncY()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $Y = $this->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $incY = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incY must be greater than 0.');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumIllegalBufferY()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $Y = $this->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $YY = new \stdClass();
        $this->expectException(TypeError::class);
        $this->argExpectExceptionMessage('Y must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumOverflowBufferYwithSize()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $Y = $this->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $YY = $this->array([1,2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferY');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumOverflowBufferYwithOffsetY()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $Y = $this->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $offsetX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumOverflowBufferYwithIncY()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $Y = $this->zeros([3],NDArray::float32);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $incY = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferY');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    public function testcumsumUnmatchDataType()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        $Y = $this->zeros([3],NDArray::float64);
        [$n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY] =
            $this->translate_cumsum($X,false,false,$Y);

        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Unmatch data type for X and Y');
        $math->cumsum($n,$XX,$offsetX,$incX,$exclusive,$reverse,$YY,$offsetY,$incY);
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testzerosNormal($params)
    {
        extract($params);
        $math = $this->getMath();

        $X = $this->array([1,2,4,9],$dtype);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $math->zeros($N,$XX,$offX,$incX);

        $this->assertEquals([0,0,0,0],$X->toArray());
    }

    public function testzerosMinusN()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->zeros($N,$XX,$offX,$incX);
    }

    public function testzerosMinusOffsetX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $math->zeros($N,$XX,$offX,$incX);
    }

    public function testzerosMinusIncX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->zeros($N,$XX,$offX,$incX);
    }

    public function testzerosIllegalBufferX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $math->zeros($N,$XX,$offX,$incX);
    }

    public function testzerosOverflowBufferXwithSize()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $XX = $this->array([1,2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->zeros($N,$XX,$offX,$incX);
    }

    public function testzerosOverflowBufferXwithOffsetX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->zeros($N,$XX,$offX,$incX);
    }

    public function testzerosOverflowBufferXwithIncX()
    {
        $math = $this->getMath();

        $X = $this->array([1,2,3]);
        [$N,$XX,$offX,$incX] =
            $this->translate_square($X);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->zeros($N,$XX,$offX,$incX);
    }

    /**
    * @dataProvider providerDtypesFloatsInt3264
    */
    public function testTransposeFloat1DNormal($params)
    {
        extract($params);
        $math = $this->getMath();

        $A = $this->array([1,2,4,9],$dtype);
        $B = $this->zerosLike($A);
        [
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB,
        ] = $this->translate_transpose($A,[0],$B);

        $math->transpose(
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB
        );

        $this->assertEquals([1,2,4,9],$B->toArray());
    }

    /**
    * @dataProvider providerDtypesFloatsInt3264
    */
    public function testTransposeFloat2DNormal($params)
    {
        extract($params);
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6]],$dtype);
        $B = $this->zeros([3,2],$dtype);
        [
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB,
        ] = $this->translate_transpose($A,[1,0],$B);

        $math->transpose(
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB
        );

        $this->assertEquals([[1,4],[2,5],[3,6]],$B->toArray());
    }

    /**
    * @dataProvider providerDtypesFloatsInt3264
    */
    public function testTransposeFloat3DNormal($params)
    {
        extract($params);
        $math = $this->getMath();

        $A = $this->array([
            [[0,1,2,3],
             [4,5,6,7],
             [8,9,10,11]],
            [[12,13,14,15],
             [16,17,18,19],
             [20,21,22,23]],
        ],$dtype);
        $B = $this->zeros([4,3,2],$dtype);
        [
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB,
        ] = $this->translate_transpose($A,[2,1,0],$B);

        $math->transpose(
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB
        );

        $this->assertEquals([
           [[ 0., 12.],
            [ 4., 16.],
            [ 8., 20.]],
    
           [[ 1., 13.],
            [ 5., 17.],
            [ 9., 21.]],
    
           [[ 2., 14.],
            [ 6., 18.],
            [10., 22.]],
    
           [[ 3., 15.],
            [ 7., 19.],
            [11., 23.]]            
        ],$B->toArray());
    }

    /**
    * @dataProvider providerDtypesFloatsInt3264
    */
    public function testTransposeFloat3DWithPerm($params)
    {
        extract($params);
        $math = $this->getMath();

        $A = $this->array([
            [[0,1,2,3],
             [4,5,6,7],
             [8,9,10,11]],
            [[12,13,14,15],
             [16,17,18,19],
             [20,21,22,23]],
        ],$dtype);
        $B = $this->zeros([2,4,3],$dtype);
        [
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB,
        ] = $this->translate_transpose($A,[0,2,1],$B);

        $math->transpose(
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB
        );

        $this->assertEquals([
            [[ 0.,  4.,  8.],
             [ 1.,  5.,  9.],
             [ 2.,  6., 10.],
             [ 3.,  7., 11.]],
    
            [[12., 16., 20.],
             [13., 17., 21.],
             [14., 18., 22.],
             [15., 19., 23.]]
        ],$B->toArray());
    }

    public function testTransposeShapeDtypeError()
    {
        $math = $this->getMath();
    
        $A = $this->array([
            [[0,1,2,3],
             [4,5,6,7],
             [8,9,10,11]],
            [[12,13,14,15],
             [16,17,18,19],
             [20,21,22,23]],
        ],NDArray::int32);
        $B = $this->zeros([2,4,3],NDArray::int32);
        [
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB,
        ] = $this->translate_transpose($A,[0,2,1],$B);
        $sourceShape = $this->array($A->shape(),NDArray::float32)->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('data type of shape buffer must be int32.');
        $math->transpose(
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB
        );
    }

    public function testTransposeShapeValueError()
    {
        $math = $this->getMath();
    
        $A = $this->array([
            [[0,1,2,3],
             [4,5,6,7],
             [8,9,10,11]],
            [[12,13,14,15],
             [16,17,18,19],
             [20,21,22,23]],
        ],NDArray::int32);
        $B = $this->zeros([2,4,3],NDArray::int32);
        [
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB,
        ] = $this->translate_transpose($A,[0,2,1],$B);
        $sourceShape[1] = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('shape values must be greater than 0.');
        $math->transpose(
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB
        );
    }

    public function testTransposePermSizeError()
    {
        $math = $this->getMath();
    
        $A = $this->array([
            [[0,1,2,3],
             [4,5,6,7],
             [8,9,10,11]],
            [[12,13,14,15],
             [16,17,18,19],
             [20,21,22,23]],
        ],NDArray::int32);
        $B = $this->zeros([2,4,3],NDArray::int32);
        [
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB,
        ] = $this->translate_transpose($A,[0,2,1,3],$B);
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('matrix shape and perm must be same size.');
        $math->transpose(
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB
        );
    }

    public function testTransposePermDtypeError()
    {
        $math = $this->getMath();
    
        $A = $this->array([
            [[0,1,2,3],
             [4,5,6,7],
             [8,9,10,11]],
            [[12,13,14,15],
             [16,17,18,19],
             [20,21,22,23]],
        ],NDArray::int32);
        $B = $this->zeros([2,4,3],NDArray::int32);
        [
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB,
        ] = $this->translate_transpose($A,[0,2,1],$B);
        $permBuf = $this->array([0,2,1],NDArray::float32)->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('data type of perm buffer must be int32.');
        $math->transpose(
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB
        );
    }

    public function testTransposeMatrixASizeError()
    {
        $math = $this->getMath();
    
        $A = $this->array([
            [[0,1,2,3],
             [4,5,6,7],
             [8,9,10,11]],
            [[12,13,14,15],
             [16,17,18,19],
             [20,21,22,23]],
        ],NDArray::int32);
        $B = $this->zeros([2,4,3],NDArray::int32);
        [
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB,
        ] = $this->translate_transpose($A,[0,2,1],$B);
        $sourceShape[0] = 3;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferA.');
        $math->transpose(
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB
        );
    }

    public function testTransposeMatrixBSizeError()
    {
        $math = $this->getMath();
    
        $A = $this->array([
            [[0,1,2,3],
             [4,5,6,7],
             [8,9,10,11]],
            [[12,13,14,15],
             [16,17,18,19],
             [20,21,22,23]],
        ],NDArray::int32);
        $B = $this->zeros([2,4,3],NDArray::int32);
        [
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB,
        ] = $this->translate_transpose($A,[0,2,1],$B);
        $BB = $this->zeros([2,4,2],NDArray::int32)->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferB.');
        $math->transpose(
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB
        );
    }
     
    public function testTransposeMatrixABDtypeError()
    {
        $math = $this->getMath();
    
        $A = $this->array([
            [[0,1,2,3],
             [4,5,6,7],
             [8,9,10,11]],
            [[12,13,14,15],
             [16,17,18,19],
             [20,21,22,23]],
        ],NDArray::int32);
        $B = $this->zeros([2,4,3],NDArray::float32);
        [
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB,
        ] = $this->translate_transpose($A,[0,2,1],$B);
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Unmatch data type for A and B.');
        $math->transpose(
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB
        );
    }

    public function testTransposeDuplicatePerm()
    {
        $math = $this->getMath();
    
        $A = $this->array([
            [[0,1,2,3],
             [4,5,6,7],
             [8,9,10,11]],
            [[12,13,14,15],
             [16,17,18,19],
             [20,21,22,23]],
        ],NDArray::float32);
        $B = $this->zeros([2,4,3],NDArray::float32);
        [
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB,
        ] = $this->translate_transpose($A,[0,2,0],$B);
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Perm contained duplicate axis');
        $math->transpose(
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB
        );
    }

    public function testTransposeOutOfAxisPerm()
    {
        $math = $this->getMath();
    
        $A = $this->array([
            [[0,1,2,3],
             [4,5,6,7],
             [8,9,10,11]],
            [[12,13,14,15],
             [16,17,18,19],
             [20,21,22,23]],
        ],NDArray::float32);
        $B = $this->zeros([2,4,3],NDArray::float32);
        [
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB,
        ] = $this->translate_transpose($A,[0,2,3],$B);
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('perm contained an out-of-bounds axis');
        $math->transpose(
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB
        );
    }

    public function testTransposefloatMatrixAOffset()
    {
        $math = $this->getMath();
    
        $A = $this->array([
            [[0,0,0,0],
             [0,0,0,0],
             [0,0,0,0]],
            [[0,1,2,3],
             [4,5,6,7],
             [8,9,10,11]],
            [[12,13,14,15],
             [16,17,18,19],
             [20,21,22,23]],
        ],NDArray::float32);
        $A = $A[R(1,3)];
        $B = $this->zeros([4,3,2],NDArray::float32);
        [
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB,
        ] = $this->translate_transpose($A,[2,1,0],$B);
        $this->assertEquals(3*4,$offsetA);
        $math->transpose(
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB
        );

        $this->assertEquals([
            [[ 0., 12.],
             [ 4., 16.],
             [ 8., 20.]],
         
            [[ 1., 13.],
             [ 5., 17.],
             [ 9., 21.]],
         
            [[ 2., 14.],
             [ 6., 18.],
             [10., 22.]],
         
            [[ 3., 15.],
             [ 7., 19.],
             [11., 23.]]            
         ],$B->toArray());
    }

    public function testTransposefloatMatrixBOffset()
    {
        $math = $this->getMath();
    
        $A = $this->array([
            [[0,1,2,3],
             [4,5,6,7],
             [8,9,10,11]],
            [[12,13,14,15],
             [16,17,18,19],
             [20,21,22,23]],
        ],NDArray::float32);
        $origB = $this->zeros([5,3,2],NDArray::float32);
        $B = $origB[R(1,5)];
        [
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB,
        ] = $this->translate_transpose($A,[2,1,0],$B);
        $math->transpose(
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB
        );

        $this->assertEquals([
            [[ 0.,  0.],
             [ 0.,  0.],
             [ 0.,  0.]],

            [[ 0., 12.],
             [ 4., 16.],
             [ 8., 20.]],
         
            [[ 1., 13.],
             [ 5., 17.],
             [ 9., 21.]],
         
            [[ 2., 14.],
             [ 6., 18.],
             [10., 22.]],
         
            [[ 3., 15.],
             [ 7., 19.],
             [11., 23.]]            
         ],$origB->toArray());
    }

    public function testTransposeDoubleMatrixAOffset()
    {
        $math = $this->getMath();
    
        $A = $this->array([
            [[0,0,0,0],
             [0,0,0,0],
             [0,0,0,0]],
            [[0,1,2,3],
             [4,5,6,7],
             [8,9,10,11]],
            [[12,13,14,15],
             [16,17,18,19],
             [20,21,22,23]],
        ],NDArray::float64);
        $A = $A[R(1,3)];
        $B = $this->zeros([4,3,2],NDArray::float64);
        [
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB,
        ] = $this->translate_transpose($A,[2,1,0],$B);
        $this->assertEquals(3*4,$offsetA);
        $math->transpose(
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB
        );

        $this->assertEquals([
            [[ 0., 12.],
             [ 4., 16.],
             [ 8., 20.]],
         
            [[ 1., 13.],
             [ 5., 17.],
             [ 9., 21.]],
         
            [[ 2., 14.],
             [ 6., 18.],
             [10., 22.]],
         
            [[ 3., 15.],
             [ 7., 19.],
             [11., 23.]]            
         ],$B->toArray());
    }

    public function testTransposeDoubleMatrixBOffset()
    {
        $math = $this->getMath();
    
        $A = $this->array([
            [[0,1,2,3],
             [4,5,6,7],
             [8,9,10,11]],
            [[12,13,14,15],
             [16,17,18,19],
             [20,21,22,23]],
        ],NDArray::float64);
        $origB = $this->zeros([5,3,2],NDArray::float64);
        $B = $origB[R(1,5)];
        [
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB,
        ] = $this->translate_transpose($A,[2,1,0],$B);
        $math->transpose(
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB
        );

        $this->assertEquals([
            [[ 0.,  0.],
             [ 0.,  0.],
             [ 0.,  0.]],

            [[ 0., 12.],
             [ 4., 16.],
             [ 8., 20.]],
         
            [[ 1., 13.],
             [ 5., 17.],
             [ 9., 21.]],
         
            [[ 2., 14.],
             [ 6., 18.],
             [10., 22.]],
         
            [[ 3., 15.],
             [ 7., 19.],
             [11., 23.]]            
         ],$origB->toArray());
    }

    public function testTransposeintMatrixAOffset()
    {
        $math = $this->getMath();
    
        $A = $this->array([
            [[0,0,0,0],
             [0,0,0,0],
             [0,0,0,0]],
            [[0,1,2,3],
             [4,5,6,7],
             [8,9,10,11]],
            [[12,13,14,15],
             [16,17,18,19],
             [20,21,22,23]],
        ],NDArray::int32);
        $A = $A[R(1,3)];
        $B = $this->zeros([4,3,2],NDArray::int32);
        [
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB,
        ] = $this->translate_transpose($A,[2,1,0],$B);
        $this->assertEquals(3*4,$offsetA);
        $math->transpose(
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB
        );

        $this->assertEquals([
            [[ 0., 12.],
             [ 4., 16.],
             [ 8., 20.]],
         
            [[ 1., 13.],
             [ 5., 17.],
             [ 9., 21.]],
         
            [[ 2., 14.],
             [ 6., 18.],
             [10., 22.]],
         
            [[ 3., 15.],
             [ 7., 19.],
             [11., 23.]]            
         ],$B->toArray());
    }

    public function testTransposeintMatrixBOffset()
    {
        $math = $this->getMath();
    
        $A = $this->array([
            [[0,1,2,3],
             [4,5,6,7],
             [8,9,10,11]],
            [[12,13,14,15],
             [16,17,18,19],
             [20,21,22,23]],
        ],NDArray::int32);
        $origB = $this->zeros([5,3,2],NDArray::int32);
        $B = $origB[R(1,5)];
        [
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB,
        ] = $this->translate_transpose($A,[2,1,0],$B);
        $math->transpose(
            $sourceShape,
            $permBuf,
            $AA, $offsetA,
            $BB, $offsetB
        );

        $this->assertEquals([
            [[ 0.,  0.],
             [ 0.,  0.],
             [ 0.,  0.]],

            [[ 0., 12.],
             [ 4., 16.],
             [ 8., 20.]],
         
            [[ 1., 13.],
             [ 5., 17.],
             [ 9., 21.]],
         
            [[ 2., 14.],
             [ 6., 18.],
             [10., 22.]],
         
            [[ 3., 15.],
             [ 7., 19.],
             [11., 23.]]            
         ],$origB->toArray());
    }
    

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testBandpartNormal($params)
    {
        extract($params);
        $math = $this->getMath();

        // under
        $A = $this->ones([2,3,3],$dtype);
        [
            $m,$n,$k,
            $AA, $offsetA,
            $lower,$upper
        ] = $this->translate_bandpart($A,0,-1);
        $math->bandpart(
            $m,$n,$k,
            $AA, $offsetA,
            $lower,$upper
        );
        $this->assertEquals([
            [[1,1,1],
             [0,1,1],
             [0,0,1]],
            [[1,1,1],
             [0,1,1],
             [0,0,1]],
        ],$A->toArray());

        $A = $this->ones([2,3,3],$dtype);
        [
            $m,$n,$k,
            $AA, $offsetA,
            $lower,$upper
        ] = $this->translate_bandpart($A,0,1);
        $math->bandpart(
            $m,$n,$k,
            $AA, $offsetA,
            $lower,$upper
        );
        $this->assertEquals([
            [[1,1,0],
             [0,1,1],
             [0,0,1]],
            [[1,1,0],
             [0,1,1],
             [0,0,1]],
        ],$A->toArray());

        // upper
        $A = $this->ones([2,3,3],$dtype);
        [
            $m,$n,$k,
            $AA, $offsetA,
            $lower,$upper
        ] = $this->translate_bandpart($A,-1,0);
        $math->bandpart(
            $m,$n,$k,
            $AA, $offsetA,
            $lower,$upper
        );
        $this->assertEquals([
            [[1,0,0],
             [1,1,0],
             [1,1,1]],
            [[1,0,0],
             [1,1,0],
             [1,1,1]],
        ],$A->toArray());

        $A = $this->ones([2,3,3],$dtype);
        [
            $m,$n,$k,
            $AA, $offsetA,
            $lower,$upper
        ] = $this->translate_bandpart($A,1,0);
        $math->bandpart(
            $m,$n,$k,
            $AA, $offsetA,
            $lower,$upper
        );
        $this->assertEquals([
            [[1,0,0],
             [1,1,0],
             [0,1,1]],
            [[1,0,0],
             [1,1,0],
             [0,1,1]],
        ],$A->toArray());
    }

    public function testBandpartParallel()
    {
        $math = $this->getMath();

        // m > n
        $A = $this->ones([4,3,3]);
        [
            $m,$n,$k,
            $AA, $offsetA,
            $lower,$upper
        ] = $this->translate_bandpart($A,0,-1);
        $math->bandpart(
            $m,$n,$k,
            $AA, $offsetA,
            $lower,$upper
        );
        $this->assertEquals([
            [[1,1,1],
             [0,1,1],
             [0,0,1]],
            [[1,1,1],
             [0,1,1],
             [0,0,1]],
            [[1,1,1],
             [0,1,1],
             [0,0,1]],
            [[1,1,1],
             [0,1,1],
             [0,0,1]],
        ],$A->toArray());

        // m < n
        $A = $this->ones([2,3,3]);
        [
            $m,$n,$k,
            $AA, $offsetA,
            $lower,$upper
        ] = $this->translate_bandpart($A,0,-1);
        $math->bandpart(
            $m,$n,$k,
            $AA, $offsetA,
            $lower,$upper
        );
        $this->assertEquals([
            [[1,1,1],
             [0,1,1],
             [0,0,1]],
            [[1,1,1],
             [0,1,1],
             [0,0,1]],
        ],$A->toArray());

    }

    public function testBandpartOffset()
    {
        $math = $this->getMath();

        $ORGA = $this->ones([2,3,3]);
        $A = $ORGA[R(1,2)];
        [
            $m,$n,$k,
            $AA, $offsetA,
            $lower,$upper
        ] = $this->translate_bandpart($A,0,-1);
        $math->bandpart(
            $m,$n,$k,
            $AA, $offsetA,
            $lower,$upper
        );
        $this->assertEquals([
            [[1,1,1],
             [1,1,1],
             [1,1,1]],
            [[1,1,1],
             [0,1,1],
             [0,0,1]],
        ],$ORGA->toArray());
    }

    public function testBandpartOverSize()
    {
        $math = $this->getMath();

        $A = $this->ones([2,3,3]);
        [
            $m,$n,$k,
            $AA, $offsetA,
            $lower,$upper
        ] = $this->translate_bandpart($A,0,-1);
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferA');
        $math->bandpart(
            $m,$n,$k+1,
            $AA, $offsetA,
            $lower,$upper
        );
    }

    public function testBandpartUnsupportedDtype()
    {
        $math = $this->getMath();

        $A = $this->ones([2,3,3],NDArray::int32);
        [
            $m,$n,$k,
            $AA, $offsetA,
            $lower,$upper
        ] = $this->translate_bandpart($A,0,-1);
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Unsupported data type.');
        $math->bandpart(
            $m,$n,$k,
            $AA, $offsetA,
            $lower,$upper
        );
    }


    public function testGatherAxisNullNormal()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $this->array([0,2],NDArray::int32);
        $B = $this->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([[1,2,3],[7,8,9]],$B->toArray());

        $A = $this->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float64);
        $X = $this->array([0,2],NDArray::int64);
        $B = $this->array([[0,0,0],[0,0,0]],NDArray::float64);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([[1,2,3],[7,8,9]],$B->toArray());

        $A = $this->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::int64);
        $X = $this->array([0,2],NDArray::int64);
        $B = $this->array([[0,0,0],[0,0,0]],NDArray::int64);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([[1,2,3],[7,8,9]],$B->toArray());

        $A = $this->array([1,2,3,4],NDArray::float32);
        $X = $this->array([0,2],NDArray::int32);
        $B = $this->array([0,0],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([1,3],$B->toArray());

        $A = $this->array([1,2,3,4],NDArray::float64);
        $X = $this->array([0,2],NDArray::int64);
        $B = $this->array([0,0],NDArray::float64);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([1,3],$B->toArray());

        $A = $this->array([1,2,3,4],NDArray::int64);
        $X = $this->array([0,2],NDArray::int64);
        $B = $this->array([0,0],NDArray::int64);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([1,3],$B->toArray());
    }

    public function testGatherAxisNullAddMode()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $this->array([0,2],NDArray::int32);
        $B = $this->array([[1,1,1],[1,1,1]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);
        $addMode = true;

        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([[2,3,4],[8,9,10]],$B->toArray());

        $A = $this->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float64);
        $X = $this->array([0,2],NDArray::int64);
        $B = $this->array([[1,1,1],[1,1,1]],NDArray::float64);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);
        $addMode = true;

        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([[2,3,4],[8,9,10]],$B->toArray());

        $A = $this->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::int64);
        $X = $this->array([0,2],NDArray::int64);
        $B = $this->array([[1,1,1],[1,1,1]],NDArray::int64);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);
        $addMode = true;

        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([[2,3,4],[8,9,10]],$B->toArray());

        $A = $this->array([1,2,3,4],NDArray::float32);
        $X = $this->array([0,2],NDArray::int32);
        $B = $this->array([1,1],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);
        $addMode = true;

        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([2,4],$B->toArray());

        $A = $this->array([1,2,3,4],NDArray::float64);
        $X = $this->array([0,2],NDArray::int64);
        $B = $this->array([1,1],NDArray::float64);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);
        $addMode = true;

        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([2,4],$B->toArray());

        $A = $this->array([1,2,3,4],NDArray::int64);
        $X = $this->array([0,2],NDArray::int64);
        $B = $this->array([1,1],NDArray::int64);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);
        $addMode = true;

        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([2,4],$B->toArray());
    }

    public function testGatherAxisNullReverse()
    {
        $math = $this->getMath();

        $A = $this->array([[0,0,0],[0,0,0],[0,0,0],[0,0,0]],NDArray::float32);
        $X = $this->array([0,2],NDArray::int32);
        $B = $this->array([[1,2,3],[7,8,9]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);
        $reverse = true;

        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([[1,2,3],[0,0,0],[7,8,9],[0,0,0]],$A->toArray());

        $A = $this->array([[0,0,0],[0,0,0],[0,0,0],[0,0,0]],NDArray::float64);
        $X = $this->array([0,2],NDArray::int64);
        $B = $this->array([[1,2,3],[7,8,9]],NDArray::float64);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);
        $reverse = true;

        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([[1,2,3],[0,0,0],[7,8,9],[0,0,0]],$A->toArray());

        $A = $this->array([[0,0,0],[0,0,0],[0,0,0],[0,0,0]],NDArray::int64);
        $X = $this->array([0,2],NDArray::int64);
        $B = $this->array([[1,2,3],[7,8,9]],NDArray::int64);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);
        $reverse = true;

        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([[1,2,3],[0,0,0],[7,8,9],[0,0,0]],$A->toArray());

        $A = $this->array([0,0,0,0],NDArray::float32);
        $X = $this->array([0,2],NDArray::int32);
        $B = $this->array([1,3],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);
        $reverse = true;

        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([1,0,3,0],$A->toArray());

        $A = $this->array([0,0,0,0],NDArray::float64);
        $X = $this->array([0,2],NDArray::int64);
        $B = $this->array([1,3],NDArray::float64);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);
        $reverse = true;

        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([1,0,3,0],$A->toArray());

        $A = $this->array([0,0,0,0],NDArray::int64);
        $X = $this->array([0,2],NDArray::int64);
        $B = $this->array([1,3],NDArray::int64);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);
        $reverse = true;

        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([1,0,3,0],$A->toArray());
    }

    public function testGatherAxisNullLabelNumberOutOfBounds1()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $this->array([0,4],NDArray::int32);
        $B = $this->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Label number is out of bounds.');
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testGatherAxisNullLabelNumberOutOfBounds2()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $this->array([0,-1],NDArray::int32);
        $B = $this->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Label number is out of bounds.');
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testGatherAxisNullMinusN()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $this->array([0,1],NDArray::int32);
        $B = $this->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $n = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testGatherAxisNullMinusK()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $this->array([0,1],NDArray::int32);
        $B = $this->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $k = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument k must be greater than 0.');
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testGatherAxisNullMinusNumClass()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $this->array([0,1],NDArray::int32);
        $B = $this->array([[0,0,0],[0,0,0]],NDArray::float32);
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
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $this->array([0,1],NDArray::int32);
        $B = $this->array([[0,0,0],[0,0,0]],NDArray::float32);
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
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $this->array([0,1],NDArray::int32);
        $B = $this->array([[0,0,0],[0,0,0]],NDArray::float32);
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
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $this->array([0,1],NDArray::int32);
        $B = $this->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $AA = $this->array([1,2,3,4,5,6,7,8,9,10,11])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix A specification too large for buffer');
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testGatherAxisNullOverflowBufferAwithOffset()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $this->array([0,1],NDArray::int32);
        $B = $this->array([[0,0,0],[0,0,0]],NDArray::float32);
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
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $this->array([0,1],NDArray::int32);
        $B = $this->array([[0,0,0],[0,0,0]],NDArray::float32);
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
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $this->array([0,1],NDArray::int32);
        $B = $this->array([[0,0,0],[0,0,0]],NDArray::float32);
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
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $this->array([0,1],NDArray::int32);
        $B = $this->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $XX = $this->array([1])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix X specification too large for buffer.');
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testGatherAxisNullOverflowBufferXwithOffsetX()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $this->array([0,1],NDArray::int32);
        $B = $this->array([[0,0,0],[0,0,0]],NDArray::float32);
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
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $this->array([0,1],NDArray::int32);
        $B = $this->array([[0,0,0],[0,0,0]],NDArray::float32);
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
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $this->array([0,1],NDArray::int32);
        $B = $this->array([[0,0,0],[0,0,0]],NDArray::float32);
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
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $this->array([0,1],NDArray::int32);
        $B = $this->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $BB = $this->array([0])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix B specification too large for buffer.');
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testGatherAxisNullOverflowBufferBwithOffsetB()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],NDArray::float32);
        $X = $this->array([0,1],NDArray::int32);
        $B = $this->array([[0,0,0],[0,0,0]],NDArray::float32);
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=null,$B,$A->dtype());
        $this->assertFalse($reduce);

        $offB = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix B specification too large for buffer.');
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testReduceGatherAxis1Normal($params)
    {
        extract($params);
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6]],$dtype);
        $X = $this->array([1,2],NDArray::int32);
        $B = $this->array([0,0],$dtype);
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=1,$B,$A->dtype());
        $this->assertTrue($reduce);

        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
        $this->assertEquals([2,6],$B->toArray());
    }

    public function testReduceGatherAxis1LabelNumberOutOfBounds1()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([1,3]);
        $B = $this->array([0,0]);
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=1,$B,$A->dtype());
        $this->assertTrue($reduce);

        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Label number is out of bounds.');
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testReduceGatherAxis1LabelNumberOutOfBounds2()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([1,-1]);
        $B = $this->array([0,0]);
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=1,$B,$A->dtype());
        $this->assertTrue($reduce);

        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Label number is out of bounds.');
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testReduceGatherAxis1MinusM()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([1,-1]);
        $B = $this->array([0,0]);
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=1,$B,$A->dtype());
        $this->assertTrue($reduce);

        $m = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testReduceGatherAxis1MinusN()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([1,-1]);
        $B = $this->array([0,0]);
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=1,$B,$A->dtype());
        $this->assertTrue($reduce);

        $n = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testReduceGatherAxis1MinusOffsetA()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([1,-1]);
        $B = $this->array([0,0]);
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
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([1,-1]);
        $B = $this->array([0,0]);
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
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([1,-1]);
        $B = $this->array([0,0]);
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=1,$B,$A->dtype());
        $this->assertTrue($reduce);

        $AA = $this->array([1,2,3,4,5])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix A specification too large for buffer.');
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testReduceGatherAxis1OverflowBufferAwithOffset()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([1,-1]);
        $B = $this->array([0,0]);
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
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([1,-1]);
        $B = $this->array([0,0]);
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
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([1,-1]);
        $B = $this->array([0,0]);
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
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([1,-1]);
        $B = $this->array([0,0]);
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=1,$B,$A->dtype());
        $this->assertTrue($reduce);
        $XX = $this->array([1])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix X specification too large for buffer.');
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testReduceGatherAxis1OverflowBufferXwithOffsetX()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([1,-1]);
        $B = $this->array([0,0]);
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
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([1,-1]);
        $B = $this->array([0,0]);
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
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([1,-1]);
        $B = $this->array([0,0]);
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
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([1,-1]);
        $B = $this->array([0,0]);
        [$reduce,$reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB]
            = $this->translate_gather($scatterAdd=false,$A,$X,$axis=1,$B,$A->dtype());
        $this->assertTrue($reduce);

        $BB = $this->array([0])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix B specification too large for buffer.');
        $math->reduceGather($reverse,$addMode,$m,$n,$numClass,$XX,$offX,$AA,$offA,$BB,$offB);
    }

    public function testReduceGatherAxis1OverflowBufferXwithOffsetB()
    {
        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([1,-1]);
        $B = $this->array([0,0]);
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
        $math = $this->getMath();
        // float32
        $numClass = 4;
        $X = $this->array([0,2],NDArray::int64);
        $A = $this->array([[1,2,3],[7,8,9]],NDArray::float32);
        $B = $this->zeros([4,3],$A->dtype());
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
        $X = $this->array([0,2],NDArray::int64);
        $A = $this->array([[1,2,3],[7,8,9]],NDArray::float64);
        $B = $this->zeros([4,3],$A->dtype());
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
        $X = $this->array([0,2],NDArray::int64);
        $A = $this->array([[1,2,3],[7,8,9]],NDArray::int64);
        $B = $this->zeros([4,3],$A->dtype());
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
        $X = $this->array([0,2],NDArray::int64);
        $A = $this->array([[1,2,3],[7,8,9]],NDArray::uint8);
        $B = $this->zeros([4,3],$A->dtype());
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
        $X = $this->array([0,2],NDArray::int64);
        $A = $this->array([1,3],NDArray::float32);
        $B = $this->zeros([4],$A->dtype());
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA]
            = $this->translate_scatter($X,$A,$numClass,$axis=null,$B);
        $this->assertFalse($reduce);
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA);
        $this->assertEquals(
           [1,0,3,0],
            $B->toArray()
        );
        // int32
        $X = $this->array([0,2],NDArray::int64);
        $A = $this->array([1,3],NDArray::int32);
        $B = $this->zeros([4],$A->dtype());
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA]
            = $this->translate_scatter($X,$A,$numClass,$axis=null,$B);
        $this->assertFalse($reduce);
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA);
        $this->assertEquals(
           [1,0,3,0],
            $B->toArray()
        );
        // float64
        $X = $this->array([0,2],NDArray::int64);
        $A = $this->array([1,3],NDArray::float64);
        $B = $this->zeros([4],$A->dtype());
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA]
            = $this->translate_scatter($X,$A,$numClass,$axis=null,$B);
        $this->assertFalse($reduce);
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA);
        $this->assertEquals(
           [1,0,3,0],
            $B->toArray()
        );
        // int64
        $X = $this->array([0,2],NDArray::int64);
        $A = $this->array([1,3],NDArray::int64);
        $B = $this->zeros([4],$A->dtype());
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA]
            = $this->translate_scatter($X,$A,$numClass,$axis=null,$B);
        $this->assertFalse($reduce);
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA);
        $this->assertEquals(
           [1,0,3,0],
            $B->toArray()
        );
        // uint8
        $X = $this->array([0,2],NDArray::int64);
        $A = $this->array([252,254],NDArray::uint8);
        $B = $this->zeros([4],$A->dtype());
        [$reduce,$reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA]
            = $this->translate_scatter($X,$A,$numClass,$axis=null,$B);
        $this->assertFalse($reduce);
        $math->gather($reverse,$addMode,$n,$k,$numClass,$XX,$offX,$BB,$offB,$AA,$offA);
        $this->assertEquals(
           [252,0,254,0],
            $B->toArray()
        );
        // x=uint8
        $X = $this->array([0,255],NDArray::uint8);
        $A = $this->array([252,254],NDArray::uint8);
        $B = $this->zeros([256],$A->dtype());
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
        $math = $this->getMath();
        $numClass = 3;
        $X = $this->array([0,1,2,0],NDArray::int32);
        $A = $this->array([1,5,9,10],NDArray::float32);
        $B = $this->zeros([4,3],$A->dtype());
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

        $X = $this->array([0,1,2,0],NDArray::int64);
        $B = $this->zeros([4,3],$A->dtype());
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

        $X = $this->array([0,1,2,0],NDArray::float32);
        $B = $this->zeros([4,3],$A->dtype());
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

        $X = $this->array([0,1,2,0],NDArray::float64);
        $B = $this->zeros([4,3],$A->dtype());
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

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testupdateAddOnehotNormal($params)
    {
        extract($params);
        $math = $this->getMath();
        $X = $this->array([1, 2],NDArray::int32);
        $Y = $this->array([[10,10,10],[10,10,10]],$dtype);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
        $this->assertEquals([[10,9,10],[10,10,9]],$Y->toArray());
    }

    public function testupdateAddOnehotOutOfboundsLabelNumber1()
    {
        $math = $this->getMath();
        $X = $this->array([1, 3]);
        $Y = $this->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Label number is out of bounds.');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotOutOfboundsLabelNumber2()
    {
        $math = $this->getMath();
        $X = $this->array([1, -1]);
        $Y = $this->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $this->expectException(RuntimeException::class);
        $this->expectExceptionMessage('Label number is out of bounds.');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotMinusM()
    {
        $math = $this->getMath();
        $X = $this->array([1, -1]);
        $Y = $this->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $m = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotMinusN()
    {
        $math = $this->getMath();
        $X = $this->array([1, -1]);
        $Y = $this->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $n = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotMinusOffsetX()
    {
        $math = $this->getMath();
        $X = $this->array([1, -1]);
        $Y = $this->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotMinusIncX()
    {
        $math = $this->getMath();
        $X = $this->array([1, -1]);
        $Y = $this->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotIllegalBufferX()
    {
        $math = $this->getMath();
        $X = $this->array([1, -1]);
        $Y = $this->array([[10,10,10],[10,10,10]]);
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
        $math = $this->getMath();
        $X = $this->array([1, -1]);
        $Y = $this->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $XX = $this->array([1])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotOverflowBufferXwithOffsetX()
    {
        $math = $this->getMath();
        $X = $this->array([1, -1]);
        $Y = $this->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotOverflowBufferXwithIncX()
    {
        $math = $this->getMath();
        $X = $this->array([1, -1]);
        $Y = $this->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotMinusOffsetY()
    {
        $math = $this->getMath();
        $X = $this->array([1, -1]);
        $Y = $this->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $offY = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetY must be greater than or equals 0.');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotMinusLdY()
    {
        $math = $this->getMath();
        $X = $this->array([1, -1]);
        $Y = $this->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $ldY = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument ldY must be greater than 0.');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotIllegalBufferY()
    {
        $math = $this->getMath();
        $X = $this->array([1, -1]);
        $Y = $this->array([[10,10,10],[10,10,10]]);
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
        $math = $this->getMath();
        $X = $this->array([1, -1]);
        $Y = $this->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $YY = $this->array([1])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferY');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotOverflowBufferYwithOffsetY()
    {
        $math = $this->getMath();
        $X = $this->array([1, -1]);
        $Y = $this->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $offY = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferY');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    public function testupdateAddOnehotOverflowBufferYwithIncY()
    {
        $math = $this->getMath();
        $X = $this->array([1, -1]);
        $Y = $this->array([[10,10,10],[10,10,10]]);
        $numClass = 3;
        [$m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY] = $this->translate_onehot(
            $X,$numClass,-1,$Y);

        $ldY = 4;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferY');
        $math->updateAddOnehot($m,$n,$a,$XX,$offX,$incX,$YY,$offY,$ldY);
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testreduceSumSameSizeNormal($params)
    {
        extract($params);
        if($this->checkSkip('reduceSum')){return;}

        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6]],$dtype);
        $X = $this->array([0,0],$dtype);
        [$m,$n,$k,$AA,$offA,$BB,$offB] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $math->reduceSum($m,$n,$k,$AA,$offA,$BB,$offB);
        $this->assertEquals([6,15],$X->toArray());
    }

    public function testreduceSumBroadcastTranspose()
    {
        if($this->checkSkip('reduceSum')){return;}

        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([0,0,0]);
        [$m,$n,$k,$AA,$offA,$BB,$offB] =
            $this->translate_reduceSum($A,$axis=0,$X);

        $math->reduceSum($m,$n,$k,$AA,$offA,$BB,$offB);
        $this->assertEquals([5,7,9],$X->toArray());
    }

    public function testreduceSumMinusM()
    {
        if($this->checkSkip('reduceSum')){return;}

        $math = $this->getMath();

        $X = $this->array([0,0]);
        $A = $this->array([[1,2,3],[4,5,6]]);
        [$m,$n,$k,$AA,$offA,$BB,$offB] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $m = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $math->reduceSum($m,$n,$k,$AA,$offA,$BB,$offB);
    }

    public function testreduceSumMinusK()
    {
        if($this->checkSkip('reduceSum')){return;}

        $math = $this->getMath();

        $X = $this->array([0,0]);
        $A = $this->array([[1,2,3],[4,5,6]]);
        [$m,$n,$k,$AA,$offA,$BB,$offB] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $k = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument k must be greater than 0.');
        $math->reduceSum($m,$n,$k,$AA,$offA,$BB,$offB);
    }

    public function testreduceSumMinusN()
    {
        if($this->checkSkip('reduceSum')){return;}

        $math = $this->getMath();

        $X = $this->array([0,0]);
        $A = $this->array([[1,2,3],[4,5,6]]);
        [$m,$n,$k,$AA,$offA,$BB,$offB] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $n = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $math->reduceSum($m,$n,$k,$AA,$offA,$BB,$offB);
    }

    public function testreduceSumMinusOffsetB()
    {
        if($this->checkSkip('reduceSum')){return;}

        $math = $this->getMath();

        $X = $this->array([0,0]);
        $A = $this->array([[1,2,3],[4,5,6]]);
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

        $math = $this->getMath();

        $X = $this->array([0,0]);
        $A = $this->array([[1,2,3],[4,5,6]]);
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

        $math = $this->getMath();

        $X = $this->array([0,0]);
        $A = $this->array([[1,2,3],[4,5,6]]);
        [$m,$n,$k,$AA,$offA,$BB,$offB] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $BB = $this->array([1])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferB');
        $math->reduceSum($m,$n,$k,$AA,$offA,$BB,$offB);
    }

    public function testreduceSumOverflowBufferBwithOffsetB()
    {
        if($this->checkSkip('reduceSum')){return;}

        $math = $this->getMath();

        $X = $this->array([0,0]);
        $A = $this->array([[1,2,3],[4,5,6]]);
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

        $math = $this->getMath();

        $X = $this->array([0,0]);
        $A = $this->array([[1,2,3],[4,5,6]]);
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

        $math = $this->getMath();

        $X = $this->array([0,0]);
        $A = $this->array([[1,2,3],[4,5,6]]);
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

        $math = $this->getMath();

        $X = $this->array([0,0]);
        $A = $this->array([[1,2,3],[4,5,6]]);
        [$m,$n,$k,$AA,$offA,$BB,$offB] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $AA = $this->array([1,2,3,4,5])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->reduceSum($m,$n,$k,$AA,$offA,$BB,$offB);
    }

    public function testreduceSumOverflowBufferXwithOffsetA()
    {
        if($this->checkSkip('reduceSum')){return;}

        $math = $this->getMath();

        $X = $this->array([0,0]);
        $A = $this->array([[1,2,3],[4,5,6]]);
        [$m,$n,$k,$AA,$offA,$BB,$offB] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $offA = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $math->reduceSum($m,$n,$k,$AA,$offA,$BB,$offB);
    }

    public function testsoftmax()
    {
        $math = $this->getMath();
        $X = $this->array([-1.0,-0.5,0.0,0.5,1.0]);
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
        $y = $this->array([
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
        $math = $this->getMath();
        $n = 5;
        $offX = 0;
        $incX = 1;
        $offY = 0;
        $incY = 1;

        $X = $this->array([-1.0,-0.5,0.0,0.5,-1.0],NDArray::float32);
        $Y = $this->array([1.0,-0.5,0.0,0.5,1.0],NDArray::float32);
        $XX = $X->buffer();
        $YY = $Y->buffer();
        $math->equal($n,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([0,1,1,1,0],$Y->toArray());

        $X = $this->array([-1.0,-0.5,0.0,0.5,-1.0],NDArray::float64);
        $Y = $this->array([1.0,-0.5,0.0,0.5,1.0],NDArray::float64);
        $XX = $X->buffer();
        $YY = $Y->buffer();
        $math->equal($n,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([0,1,1,1,0],$Y->toArray());

        $X = $this->array([-2,-1,0,1,-1],NDArray::int8);
        $Y = $this->array([ 1,-1,0,1, 1],NDArray::int8);
        $XX = $X->buffer();
        $YY = $Y->buffer();
        $math->equal($n,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([0,1,1,1,0],$Y->toArray());

        $X = $this->array([-2,-1,0,1,-1],NDArray::int16);
        $Y = $this->array([ 1,-1,0,1, 1],NDArray::int16);
        $XX = $X->buffer();
        $YY = $Y->buffer();
        $math->equal($n,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([0,1,1,1,0],$Y->toArray());

        $X = $this->array([-2,-1,0,1,-1],NDArray::int32);
        $Y = $this->array([ 1,-1,0,1, 1],NDArray::int32);
        $XX = $X->buffer();
        $YY = $Y->buffer();
        $math->equal($n,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([0,1,1,1,0],$Y->toArray());

        $X = $this->array([-2,-1,0,1,-1],NDArray::int64);
        $Y = $this->array([ 1,-1,0,1, 1],NDArray::int64);
        $XX = $X->buffer();
        $YY = $Y->buffer();
        $math->equal($n,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([0,1,1,1,0],$Y->toArray());

        $X = $this->array([false,false,true ,true,true ],NDArray::bool);
        $Y = $this->array([true ,false,true ,true,false],NDArray::bool);
        $XX = $X->buffer();
        $YY = $Y->buffer();
        $math->equal($n,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([false,true,true,true,false],$Y->toArray());
    }


    public function testastype()
    {
        $math = $this->getMath();

        #### int to any
        $X = $this->array([-1,0,1,2,3],NDArray::int32);
        $dtype = NDArray::float32;
        $Y = $this->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());

        $dtype = NDArray::float64;
        $Y = $this->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());

        $dtype = NDArray::int8;
        $Y = $this->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());

        $dtype = NDArray::int16;
        $Y = $this->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());

        $dtype = NDArray::int32;
        $Y = $this->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());

        $dtype = NDArray::int64;
        $Y = $this->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());

        $dtype = NDArray::bool;
        $Y = $this->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([true,false,true,true,true],$Y->toArray());

        #### float to any ######
        $X = $this->array([-1,0,1,2,3],NDArray::float32);
        $dtype = NDArray::float32;
        $Y = $this->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());

        $dtype = NDArray::float64;
        $Y = $this->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());

        $dtype = NDArray::int8;
        $Y = $this->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());

        $dtype = NDArray::int16;
        $Y = $this->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());

        $dtype = NDArray::int32;
        $Y = $this->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());

        $dtype = NDArray::int64;
        $Y = $this->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());

        $dtype = NDArray::bool;
        $Y = $this->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([true,false,true,true,true],$Y->toArray());

        #### bool to any ######
        $X = $this->array([true,false,true,true,true],NDArray::bool);
        $dtype = NDArray::float32;
        $Y = $this->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([1,0,1,1,1],$Y->toArray());

        $dtype = NDArray::float64;
        $Y = $this->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([1,0,1,1,1],$Y->toArray());

        $dtype = NDArray::int8;
        $Y = $this->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([1,0,1,1,1],$Y->toArray());

        $dtype = NDArray::int16;
        $Y = $this->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([1,0,1,1,1],$Y->toArray());

        $dtype = NDArray::int32;
        $Y = $this->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([1,0,1,1,1],$Y->toArray());

        $dtype = NDArray::int64;
        $Y = $this->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([1,0,1,1,1],$Y->toArray());

        $dtype = NDArray::bool;
        $Y = $this->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([true,false,true,true,true],$Y->toArray());

        #### float to unsigned ######
        $X = $this->array([-1,0,1,2,3],NDArray::float32);
        $dtype = NDArray::uint8;
        $Y = $this->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([255,0,1,2,3],$Y->toArray());

        $dtype = NDArray::uint16;
        $Y = $this->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([65535,0,1,2,3],$Y->toArray());

        $dtype = NDArray::uint32;
        $Y = $this->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([4294967295,0,1,2,3],$Y->toArray());

        $dtype = NDArray::uint64;
        $Y = $this->zeros($X->shape(),$dtype);
        [$n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY] = $this->translate_astype($X, $dtype, $Y);
        $math->astype($n,$dtype,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([-1,0,1,2,3],$Y->toArray());
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testreduceMaxSameSizeNormal($params)
    {
        extract($params);
        if($this->checkSkip('reduceMax')){return;}

        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6]],$dtype);
        $X = $this->array([0,0],$dtype);
        [$m,$n,$k,$AA,$offA,$BB,$offB] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $math->reduceMax($m,$n,$k,$AA,$offA,$BB,$offB);
        $this->assertEquals([3,6],$X->toArray());
    }

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testreduceArgMaxSameSizeNormal($params)
    {
        extract($params);
        if($this->checkSkip('reduceArgMax')){return;}

        $math = $this->getMath();

        $A = $this->array([[1,2,3],[4,5,6]],$dtype);
        $X = $this->array([0,0],NDArray::int32);
        [$m,$n,$k,$AA,$offA,$BB,$offB] =
            $this->translate_reduceSum($A,$axis=1,$X);

        $math->reduceArgMax($m,$n,$k,$AA,$offA,$BB,$offB);
        $this->assertEquals([2,2],$X->toArray());
    }


    /**
    * @dataProvider providerDtypesFloats
    */
    public function testIm2col1dNormal($params)
    {
        extract($params);
        if($this->checkSkip('im2col1d')){return;}

        $math = $this->getMath();

        $images = $this->array([1,2,3,4],$dtype);
        $cols = $this->zeros([1,2,3,1],$dtype);

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

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testIm2col2dNormal($params)
    {
        extract($params);
        if($this->checkSkip('im2col2d')){return;}

        $math = $this->getMath();

        $reverse = false;
        $batches = 1;
        $im_h = 4;
        $im_w = 4;
        $channels = 3;
        $kernel_h = 3;
        $kernel_w = 3;
        $stride_h = 1;
        $stride_w = 1;
        $padding = false;
        $channels_first = false;
        $cols_channels_first=false;
        $cols = null;
        $out_h = 2;
        $out_w = 2;
        $images = $this->arange(
            $batches*
            $im_h*$im_w*
            $channels,
            null,null,
            $dtype
        )->reshape([
            $batches,
            $im_h,
            $im_w,
            $channels
        ]);
        $cols = $this->zeros(
            [
                $batches,
                $out_h,$out_w,
                $kernel_h,$kernel_w,
                $channels,
            ],$dtype);
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

    /**
    * @dataProvider providerDtypesFloats
    */
    public function testIm2col3dNormal($params)
    {
        extract($params);
        if($this->checkSkip('im2col3d')){return;}

        $math = $this->getMath();

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
        $padding = false;
        $channels_first = false;
        $cols_channels_first=false;
        $cols = null;
        $out_d = 2;
        $out_h = 2;
        $out_w = 2;

        $images = $this->arange(
            $batches*
            $im_d*$im_h*$im_w*
            $channels,
            null,null,
            $dtype
        )->reshape([
            $batches,
            $im_d,
            $im_h,
            $im_w,
            $channels
        ]);

        $cols = $this->zeros(
            [
                $batches,
                $out_d,$out_h,$out_w,
                $kernel_d,$kernel_h,$kernel_w,
                $channels,
            ],$dtype);
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
            $this->zerosLike($cols)
            );
    }

}
