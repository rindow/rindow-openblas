<?php
namespace RindowTest\OpenBLAS\BlasTest;
if(!class_exists('RindowTest\OpenBLAS\Utils')) {
    include_once __DIR__.'/Utils.php';
}
use RindowTest\OpenBLAS\Utils;

use PHPUnit\Framework\TestCase;
use Interop\Polite\Math\Matrix\NDArray;
use Interop\Polite\Math\Matrix\BLAS;
use Rindow\Math\Matrix\MatrixOperator;
use Rindow\OpenBLAS\BLAS as OpenBLAS;
use InvalidArgumentException;
use RuntimeException;
use TypeError;

/**
 * @requires extension rindow_openblas
 */
class BlasTest extends TestCase
{
    use Utils;

    public function getBlas()
    {
        $blas = new OpenBLAS();
        return $blas;
    }

    public function getOpenBLASVersion($blas)
    {
        $config = $blas->getConfig();
        if(strpos($config,'OpenBLAS')===0) {
            $config = explode(' ',$config);
            return $config[1];
        } else {
            return '0.0.0';
        }
    }

    public function skipiamin()
    {
        $blas = $this->getBlas();
        if(version_compare($this->getOpenBLASVersion($blas),'0.3.6','>=')) {
            return false;
        }
        $this->markTestSkipped("openblas has no iamin");
        return true;
    }

    public function translate_scal(
        float $a,NDArray $X) : array
    {
        $N = $X->size();
        $XX = $X->buffer();
        $offX = $X->offset();
        return [$N,$a,$XX,$offX,1];
    }

    public function translate_axpy(
        NDArray $X,NDArray $Y, float $alpha=null) : array
    {
        if($X->shape()!=$Y->shape()) {
            $shapeError = '('.implode(',',$X->shape()).'),('.implode(',',$Y->shape()).')';
            throw new InvalidArgumentException("Unmatch shape of dimension: ".$shapeError);
        }
        $N = $X->size();
        $XX = $X->buffer();
        $offX = $X->offset();
        $N = $X->size();
        $YY = $Y->buffer();
        $offY = $Y->offset();
        return [$N,$alpha,$XX,$offX,1,$YY,$offY,1];
    }

    public function translate_dot(
        NDArray $X,NDArray $Y) : array
    {
        if($X->shape()!=$Y->shape()) {
            $shapeError = '('.implode(',',$X->shape()).'),('.implode(',',$Y->shape()).')';
            throw new InvalidArgumentException("Unmatch shape of dimension: ".$shapeError);
        }
        $N = $X->size();
        $XX = $X->buffer();
        $offX = $X->offset();
        $N = $X->size();
        $YY = $Y->buffer();
        $offY = $Y->offset();
        return [$N,$XX,$offX,1,$YY,$offY,1];
    }

    public function translate_amin(
        NDArray $X) : array
    {
        $N = $X->size();
        $XX = $X->buffer();
        $offX = $X->offset();
        return [$N,$XX,$offX,1];
    }

    public function translate_copy(
        NDArray $X,
        NDArray $Y ) : array
    {
        if($X->shape()!=$Y->shape()) {
            $shapeError = '('.implode(',',$X->shape()).'),('.implode(',',$Y->shape()).')';
            throw new InvalidArgumentException("Unmatch shape of dimension: ".$shapeError);
        }
        $N = $X->size();
        $XX = $X->buffer();
        $offX = $X->offset();
        $YY = $Y->buffer();
        $offY = $Y->offset();
        return [$N,$XX,$offX,1,$YY,$offY,1];
    }

    public function translate_gemv(
        NDArray $A,
        NDArray $X,
        float $alpha=null,
        float $beta=null,
        NDArray $Y=null,
        bool $trans=null)
    {
        if($A->ndim()!=2 || $X->ndim()!=1) {
            throw new InvalidArgumentException('"A" must be 2D-NDArray and "X" must 1D-NDArray.');
        }
        $shapeA = $A->shape();
        $shapeX = $X->shape();
        $rows = (!$trans) ? $shapeA[0] : $shapeA[1];
        $cols = (!$trans) ? $shapeA[1] : $shapeA[0];
        if($cols!=$shapeX[0]) {
            throw new InvalidArgumentException('The number of columns in "A" and The number of item in "X" must be the same');
        }
        $AA = $A->buffer();
        $XX = $X->buffer();
        $offA = $A->offset();
        $offX = $X->offset();
        $m = $shapeA[0];
        $n = $shapeA[1];
        if($alpha===null) {
            $alpha = 1.0;
        }
        if($beta===null) {
            $beta = 0.0;
        }
        if($Y!=null) {
            if($Y->ndim()!=1) {
                throw new InvalidArgumentException('"Y" must 1D-NDArray.');
            }
            $shapeY = $Y->shape();
            if($rows!=$shapeY[0]) {
                throw new InvalidArgumentException('The number of rows in "A" and The number of item in "Y" must be the same');
            }
        } else {
            $Y = $this->mo->zeros([$rows]);
        }
        $YY = $Y->buffer();
        $offY = $Y->offset();
        $trans = (!$trans) ? BLAS::NoTrans : BLAS::Trans;

        return [
            $trans,
            $m,$n,
            $alpha,
            $AA,$offA,$n,
            $XX,$offX,1,
            $beta,
            $YY,$offY,1,
        ];
    }

    public function translate_gemm(
        NDArray $A,
        NDArray $B,
        float $alpha=null,
        float $beta=null,
        NDArray $C=null,
        bool $transA=null,
        bool $transB=null)
    {
        $shapeA = $A->shape();
        if($transA) {
            $shapeA = [$shapeA[1],$shapeA[0]];
        }
        $shapeB = $B->shape();
        if($transB) {
            $shapeB = [$shapeB[1],$shapeB[0]];
        }
        if($shapeA[1]!=$shapeB[0]) {
            throw new InvalidArgumentException('The number of columns in "A" and the number of rows in "B" must be the same');
        }
        $AA = $A->buffer();
        $BB = $B->buffer();
        $offA = $A->offset();
        $offB = $B->offset();
        $M = $shapeA[0];
        $N = $shapeB[1];
        $K = $shapeA[1];

        if($alpha===null) {
            $alpha = 1.0;
        }
        if($beta===null) {
            $beta = 0.0;
        }
        if($C!=null) {
            $shapeC = $C->shape();
            if($M!=$shapeC[0] || $N!=$shapeC[1]) {
                throw new InvalidArgumentException('"A" and "C" must have the same number of rows."B" and "C" must have the same number of columns');
            }
        } else {
            $C = $this->mo->zeros([$M,$N]);
        }
        $CC = $C->buffer();
        $offC = $C->offset();

        $lda = ($transA) ? $M : $K;
        $ldb = ($transB) ? $K : $N;
        $ldc = $N;
        $transA = ($transA) ? BLAS::Trans : BLAS::NoTrans;
        $transB = ($transB) ? BLAS::Trans : BLAS::NoTrans;

        return [
            $transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc,
        ];
    }

    public function testGetNumThreads()
    {
        $blas = $this->getBlas();
        $n = $blas->getNumThreads();
        $this->assertGreaterThan(0,$n);
    }

    public function testGetNumProcs()
    {
        $blas = $this->getBlas();
        $n = $blas->getNumProcs();
        $this->assertGreaterThan(0,$n);
    }

    public function testGetConfig()
    {
        $blas = $this->getBlas();
        $s = $blas->getConfig();

        $this->assertTrue(
            strpos($s,'OpenBLAS')===0 ||
            strpos($s,'NO_LAPACKE')===0);
    }

    public function testGetCorename()
    {
        $blas = $this->getBlas();
        $s = $blas->getCorename();
        $this->assertTrue(is_string($s));
    }

    public function testScalNormal()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_scal(2,$X);

        $blas->scal($N,$alpha,$XX,$offX,$incX);
        $this->assertEquals([2,4,6],$X->toArray());
    }

    public function testScalMinusN()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_scal(2,$X);

        $N = 0;

        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $blas->scal($N,$alpha,$XX,$offX,$incX);
    }

    public function testScalMinusOffsetX()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_scal(2,$X);

        $offX = -1;

        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $blas->scal($N,$alpha,$XX,$offX,$incX);
    }

    public function testScalMinusIncX()
    {

        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_scal(2,$X);

        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $blas->scal($N,$alpha,$XX,$offX,$incX);
    }

    public function testScalIllegalBufferX()
    {

        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_scal(2,$X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $blas->scal($N,$alpha,$XX,$offX,$incX);
    }

    public function testScalOverflowBufferXwithSize()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_scal(2,$X);

        $XX = $this->array([1,2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $blas->scal($N,$alpha,$XX,$offX,$incX);
    }

    public function testScalOverflowBufferXwithOffsetX()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_scal(2,$X);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $blas->scal($N,$alpha,$XX,$offX,$incX);
    }

    public function testScalOverflowBufferXwithIncX()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        [$N,$alpha,$XX,$offX,$incX] =
            $this->translate_scal(2,$X);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $blas->scal($N,$alpha,$XX,$offX,$incX);
    }

    public function testAxpyNormal()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([10,20,30]);
        [$N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_axpy($X,$Y,2);

        $blas->axpy($N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([12,24,36],$Y->toArray());
    }

    public function testAxpyMinusN()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([10,20,30]);
        [$N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_axpy($X,$Y,2);

        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $blas->axpy($N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testAxpyMinusOffsetX()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([10,20,30]);
        [$N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_axpy($X,$Y,2);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $blas->axpy($N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testAxpyMinusIncX()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([10,20,30]);
        [$N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_axpy($X,$Y,2);

        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $blas->axpy($N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testAxpyIllegalBufferX()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([10,20,30]);
        [$N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_axpy($X,$Y,2);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $blas->axpy($N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testAxpyOverflowBufferXwithSize()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([10,20,30]);
        [$N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_axpy($X,$Y,2);

        $XX = $this->array([1,2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for buffer');
        $blas->axpy($N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testAxpyOverflowBufferXwithOffsetX()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([10,20,30]);
        [$N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_axpy($X,$Y,2);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for buffer');
        $blas->axpy($N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testAxpyOverflowBufferXwithIncX()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([10,20,30]);
        [$N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_axpy($X,$Y,2);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for buffer');
        $blas->axpy($N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testAxpyMinusOffsetY()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([10,20,30]);
        [$N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_axpy($X,$Y,2);

        $offY = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetY must be greater than or equals 0.');
        $blas->axpy($N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testAxpyMinusIncY()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([10,20,30]);
        [$N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_axpy($X,$Y,2);

        $incY = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incY must be greater than 0.');
        $blas->axpy($N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testAxpyIllegalBufferY()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([10,20,30]);
        [$N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_axpy($X,$Y,2);

        $YY = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $blas->axpy($N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testAxpyOverflowBufferYwithSize()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([10,20,30]);
        [$N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_axpy($X,$Y,2);

        $YY = $this->array([1,2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for buffer');
        $blas->axpy($N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testAxpyOverflowBufferXwithOffsetY()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([10,20,30]);
        [$N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_axpy($X,$Y,2);

        $offY = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for buffer');
        $blas->axpy($N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testAxpyOverflowBufferYwithIncY()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([10,20,30]);
        [$N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_axpy($X,$Y,2);

        $incY = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for buffer');
        $blas->axpy($N,$alpha,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testDotNormal()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([4,5,6]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_dot($X,$Y);

        $dot = $blas->dot($N,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals(32,$dot);
    }

    public function testDotMinusN()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([4,5,6]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_dot($X,$Y);

        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $dot = $blas->dot($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testDotMinusOffsetX()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([4,5,6]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_dot($X,$Y);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $dot = $blas->dot($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testDotMinusIncX()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([4,5,6]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_dot($X,$Y);

        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $dot = $blas->dot($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testDotIllegalBufferX()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([4,5,6]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_dot($X,$Y);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $dot = $blas->dot($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testDotOverflowBufferXwithSize()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([4,5,6]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_dot($X,$Y);

        $XX = $this->array([1,2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for buffer');
        $dot = $blas->dot($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testDotOverflowBufferXwithOffsetX()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([4,5,6]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_dot($X,$Y);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for buffer');
        $dot = $blas->dot($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testDotOverflowBufferXwithIncX()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([4,5,6]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_dot($X,$Y);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for buffer');
        $dot = $blas->dot($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testDotMinusOffsetY()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([4,5,6]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_dot($X,$Y);

        $offY = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetY must be greater than or equals 0.');
        $dot = $blas->dot($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testDotMinusIncY()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([4,5,6]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_dot($X,$Y);

        $incY = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incY must be greater than 0.');
        $dot = $blas->dot($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testDotIllegalBufferY()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([4,5,6]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_dot($X,$Y);

        $YY = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $dot = $blas->dot($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testDotOverflowBufferYwithSize()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([4,5,6]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_dot($X,$Y);

        $YY = $this->array([1,2])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for buffer');
        $dot = $blas->dot($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testDotOverflowBufferXwithOffsetY()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([4,5,6]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_dot($X,$Y);

        $offY = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for buffer');
        $dot = $blas->dot($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testDotOverflowBufferYwithIncY()
    {
        $blas = $this->getBlas();

        $X = $this->array([1,2,3]);
        $Y = $this->array([4,5,6]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_dot($X,$Y);

        $incY = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for buffer');
        $dot = $blas->dot($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testAsumNormal()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,-10,-1000]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);

        $min = $blas->asum($N,$XX,$offX,$incX);
        $this->assertEquals(1110,$min);
    }

    public function testAsumMinusN()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);

        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $min = $blas->asum($N,$XX,$offX,$incX);
    }

    public function testAsumMinusOffsetX()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $min = $blas->asum($N,$XX,$offX,$incX);
    }

    public function testAsumMinusIncX()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);

        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $min = $blas->asum($N,$XX,$offX,$incX);
    }

    public function testAsumIllegalBufferX()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $min = $blas->asum($N,$XX,$offX,$incX);
    }

    public function testAsumOverflowBufferXwithSize()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);

        $XX = $this->array([100,-10])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for buffer');
        $min = $blas->asum($N,$XX,$offX,$incX);
    }

    public function testAsumOverflowBufferXwithOffsetX()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for buffer');
        $min = $blas->asum($N,$XX,$offX,$incX);
    }

    public function testAsumOverflowBufferXwithIncX()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for buffer');
        $min = $blas->asum($N,$XX,$offX,$incX);
    }

    public function testAMaxNormal()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);

        $min = $blas->iamax($N,$XX,$offX,$incX);
        $this->assertEquals(0,$min);
    }

    public function testAMaxMinusN()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);

        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $min = $blas->iamax($N,$XX,$offX,$incX);
    }

    public function testAMaxMinusOffsetX()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $min = $blas->iamax($N,$XX,$offX,$incX);
    }

    public function testAMaxMinusIncX()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);

        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $min = $blas->iamax($N,$XX,$offX,$incX);
    }

    public function testAMaxIllegalBufferX()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $min = $blas->iamax($N,$XX,$offX,$incX);
    }

    public function testAMaxOverflowBufferXwithSize()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);

        $XX = $this->array([100,-10])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for buffer');
        $min = $blas->iamax($N,$XX,$offX,$incX);
    }

    public function testAMaxOverflowBufferXwithOffsetX()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for buffer');
        $min = $blas->iamax($N,$XX,$offX,$incX);
    }

    public function testAMaxOverflowBufferXwithIncX()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for buffer');
        $min = $blas->iamax($N,$XX,$offX,$incX);
    }

    public function testAminNormal()
    {
        if($this->skipiamin()) return;
        $blas = $this->getBlas();

        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);

        $min = $blas->iamin($N,$XX,$offX,$incX);
        $this->assertEquals(2,$min);
    }

    public function testAminMinusN()
    {
        if($this->skipiamin()) return;
        $blas = $this->getBlas();

        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);

        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $min = $blas->iamin($N,$XX,$offX,$incX);
    }

    public function testAminMinusOffsetX()
    {
        if($this->skipiamin()) return;
        $blas = $this->getBlas();

        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $min = $blas->iamin($N,$XX,$offX,$incX);
    }

    public function testAminMinusIncX()
    {
        if($this->skipiamin()) return;
        $blas = $this->getBlas();

        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);

        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $min = $blas->iamin($N,$XX,$offX,$incX);
    }

    public function testAminIllegalBufferX()
    {
        if($this->skipiamin()) return;
        $blas = $this->getBlas();

        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $min = $blas->iamin($N,$XX,$offX,$incX);
    }

    public function testAminOverflowBufferXwithSize()
    {
        if($this->skipiamin()) return;
        $blas = $this->getBlas();

        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);

        $XX = $this->array([100,-10])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for buffer');
        $min = $blas->iamin($N,$XX,$offX,$incX);
    }

    public function testAminOverflowBufferXwithOffsetX()
    {
        if($this->skipiamin()) return;
        $blas = $this->getBlas();

        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for buffer');
        $min = $blas->iamin($N,$XX,$offX,$incX);
    }

    public function testAminOverflowBufferXwithIncX()
    {
        if($this->skipiamin()) return;
        $blas = $this->getBlas();

        $X = $this->array([100,-10,1]);
        [$N,$XX,$offX,$incX] =
            $this->translate_amin($X);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for buffer');
        $min = $blas->iamin($N,$XX,$offX,$incX);
    }

    public function testCopyNormal()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,10,1]);
        $Y = $this->array([0,0,0]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_copy($X,$Y);

        $blas->copy($N,$XX,$offX,$incX,$YY,$offY,$incY);
        $this->assertEquals([100,10,1],$Y->toArray());
    }

    public function testCopyMinusN()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,10,1]);
        $Y = $this->array([0,0,0]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_copy($X,$Y);

        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $blas->copy($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testCopyMinusOffsetX()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,10,1]);
        $Y = $this->array([0,0,0]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_copy($X,$Y);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $blas->copy($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testCopyMinusIncX()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,10,1]);
        $Y = $this->array([0,0,0]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_copy($X,$Y);

        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $blas->copy($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testCopyIllegalBufferX()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,10,1]);
        $Y = $this->array([0,0,0]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_copy($X,$Y);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $blas->copy($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testCopyMinusOffsetY()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,10,1]);
        $Y = $this->array([0,0,0]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_copy($X,$Y);

        $offY = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetY must be greater than or equals 0.');
        $blas->copy($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testCopyMinusIncY()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,10,1]);
        $Y = $this->array([0,0,0]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_copy($X,$Y);

        $incY = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incY must be greater than 0.');
        $blas->copy($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testCopyIllegalBufferY()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,10,1]);
        $Y = $this->array([0,0,0]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_copy($X,$Y);

        $YY = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $blas->copy($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testCopyOverflowBufferXWithSize()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,10,1]);
        $Y = $this->array([0,0,0]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_copy($X,$Y);

        $XX = $this->array([100,10])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX.');
        $blas->copy($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testCopyOverflowBufferXWithOffsetX()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,10,1]);
        $Y = $this->array([0,0,0]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_copy($X,$Y);

        $offX = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX.');
        $blas->copy($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testCopyOverflowBufferXWithIncX()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,10,1]);
        $Y = $this->array([0,0,0]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_copy($X,$Y);

        $incX = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX.');
        $blas->copy($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testCopyOverflowBufferYWithSize()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,10,1]);
        $Y = $this->array([0,0,0]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_copy($X,$Y);

        $YY = $this->array([100,10])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferY.');
        $blas->copy($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testCopyOverflowBufferXWithOffsetY()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,10,1]);
        $Y = $this->array([0,0,0]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_copy($X,$Y);

        $offY = 1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferY.');
        $blas->copy($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testCopyOverflowBufferXWithIncY()
    {
        $blas = $this->getBlas();

        $X = $this->array([100,10,1]);
        $Y = $this->array([0,0,0]);
        [$N,$XX,$offX,$incX,$YY,$offY,$incY] =
            $this->translate_copy($X,$Y);

        $incY = 2;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferY.');
        $blas->copy($N,$XX,$offX,$incX,$YY,$offY,$incY);
    }

    public function testGemvNormal()
    {
        $blas = $this->getBlas();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([100,10,1]);
        $Y = $this->zeros([2]);

        [ $trans,$m,$n,$alpha,$AA,$offA,$ldA,
          $XX,$offX,$incX,$beta,$YY,$offY,$incY] =
            $this->translate_gemv($A,$X,null,null,$Y);

        $blas->gemv(
            BLAS::RowMajor,$trans,
            $m,$n,
            $alpha,
            $AA,$offA,$ldA,
            $XX,$offX,$incX,
            $beta,
            $YY,$offY,$incY);

        $this->assertEquals(
            [123,456]
        ,$Y->toArray());
    }

    public function testGemvTranspose()
    {
        $blas = $this->getBlas();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([10,1]);
        $Y = $this->zeros([3]);

        [ $trans,$m,$n,$alpha,$AA,$offA,$ldA,
          $XX,$offX,$incX,$beta,$YY,$offY,$incY] =
            $this->translate_gemv($A,$X,null,null,$Y,true);

        $blas->gemv(
            BLAS::RowMajor,$trans,
            $m,$n,
            $alpha,
            $AA,$offA,$ldA,
            $XX,$offX,$incX,
            $beta,
            $YY,$offY,$incY);

        $this->assertEquals(
            [14,25,36]
        ,$Y->toArray());
    }

    public function testGemvMinusM()
    {
        $blas = $this->getBlas();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([100,10,1]);
        $Y = $this->zeros([2]);

        [ $trans,$m,$n,$alpha,$AA,$offA,$ldA,
          $XX,$offX,$incX,$beta,$YY,$offY,$incY] =
            $this->translate_gemv($A,$X,null,null,$Y);

        $m = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $blas->gemv(
            BLAS::RowMajor,$trans,
            $m,$n,
            $alpha,
            $AA,$offA,$ldA,
            $XX,$offX,$incX,
            $beta,
            $YY,$offY,$incY);
    }

    public function testGemvMinusN()
    {
        $blas = $this->getBlas();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([100,10,1]);
        $Y = $this->zeros([2]);

        [ $trans,$m,$n,$alpha,$AA,$offA,$ldA,
          $XX,$offX,$incX,$beta,$YY,$offY,$incY] =
            $this->translate_gemv($A,$X,null,null,$Y);

        $n = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $blas->gemv(
            BLAS::RowMajor,$trans,
            $m,$n,
            $alpha,
            $AA,$offA,$ldA,
            $XX,$offX,$incX,
            $beta,
            $YY,$offY,$incY);
    }

    public function testGemvMinusOffsetA()
    {
        $blas = $this->getBlas();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([100,10,1]);
        $Y = $this->zeros([2]);

        [ $trans,$m,$n,$alpha,$AA,$offA,$ldA,
          $XX,$offX,$incX,$beta,$YY,$offY,$incY] =
            $this->translate_gemv($A,$X,null,null,$Y);

        $offA = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetA must be greater than or equals 0.');
        $blas->gemv(
            BLAS::RowMajor,$trans,
            $m,$n,
            $alpha,
            $AA,$offA,$ldA,
            $XX,$offX,$incX,
            $beta,
            $YY,$offY,$incY);
    }

    public function testGemvMinusLdA()
    {
        $blas = $this->getBlas();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([100,10,1]);
        $Y = $this->zeros([2]);

        [ $trans,$m,$n,$alpha,$AA,$offA,$ldA,
          $XX,$offX,$incX,$beta,$YY,$offY,$incY] =
            $this->translate_gemv($A,$X,null,null,$Y);

        $ldA = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument ldA must be greater than 0.');
        $blas->gemv(
            BLAS::RowMajor,$trans,
            $m,$n,
            $alpha,
            $AA,$offA,$ldA,
            $XX,$offX,$incX,
            $beta,
            $YY,$offY,$incY);
    }

    public function testGemvIllegalBufferA()
    {
        $blas = $this->getBlas();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([100,10,1]);
        $Y = $this->zeros([2]);

        [ $trans,$m,$n,$alpha,$AA,$offA,$ldA,
          $XX,$offX,$incX,$beta,$YY,$offY,$incY] =
            $this->translate_gemv($A,$X,null,null,$Y);

        $AA = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $blas->gemv(
            BLAS::RowMajor,$trans,
            $m,$n,
            $alpha,
            $AA,$offA,$ldA,
            $XX,$offX,$incX,
            $beta,
            $YY,$offY,$incY);
    }

    public function testGemvMinusOffsetX()
    {
        $blas = $this->getBlas();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([100,10,1]);
        $Y = $this->zeros([2]);

        [ $trans,$m,$n,$alpha,$AA,$offA,$ldA,
          $XX,$offX,$incX,$beta,$YY,$offY,$incY] =
            $this->translate_gemv($A,$X,null,null,$Y);

        $offX = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetX must be greater than or equals 0.');
        $blas->gemv(
            BLAS::RowMajor,$trans,
            $m,$n,
            $alpha,
            $AA,$offA,$ldA,
            $XX,$offX,$incX,
            $beta,
            $YY,$offY,$incY);
    }

    public function testGemvMinusIncX()
    {
        $blas = $this->getBlas();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([100,10,1]);
        $Y = $this->zeros([2]);

        [ $trans,$m,$n,$alpha,$AA,$offA,$ldA,
          $XX,$offX,$incX,$beta,$YY,$offY,$incY] =
            $this->translate_gemv($A,$X,null,null,$Y);

        $incX = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incX must be greater than 0.');
        $blas->gemv(
            BLAS::RowMajor,$trans,
            $m,$n,
            $alpha,
            $AA,$offA,$ldA,
            $XX,$offX,$incX,
            $beta,
            $YY,$offY,$incY);
    }

    public function testGemvIllegalBufferX()
    {
        $blas = $this->getBlas();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([100,10,1]);
        $Y = $this->zeros([2]);

        [ $trans,$m,$n,$alpha,$AA,$offA,$ldA,
          $XX,$offX,$incX,$beta,$YY,$offY,$incY] =
            $this->translate_gemv($A,$X,null,null,$Y);

        $XX = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $blas->gemv(
            BLAS::RowMajor,$trans,
            $m,$n,
            $alpha,
            $AA,$offA,$ldA,
            $XX,$offX,$incX,
            $beta,
            $YY,$offY,$incY);
    }

    public function testGemvMinusOffsetY()
    {
        $blas = $this->getBlas();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([100,10,1]);
        $Y = $this->zeros([2]);

        [ $trans,$m,$n,$alpha,$AA,$offA,$ldA,
          $XX,$offX,$incX,$beta,$YY,$offY,$incY] =
            $this->translate_gemv($A,$X,null,null,$Y);

        $offY = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetY must be greater than or equals 0.');
        $blas->gemv(
            BLAS::RowMajor,$trans,
            $m,$n,
            $alpha,
            $AA,$offA,$ldA,
            $XX,$offX,$incX,
            $beta,
            $YY,$offY,$incY);
    }

    public function testGemvMinusIncY()
    {
        $blas = $this->getBlas();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([100,10,1]);
        $Y = $this->zeros([2]);

        [ $trans,$m,$n,$alpha,$AA,$offA,$ldA,
          $XX,$offX,$incX,$beta,$YY,$offY,$incY] =
            $this->translate_gemv($A,$X,null,null,$Y);

        $incY = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument incY must be greater than 0.');
        $blas->gemv(
            BLAS::RowMajor,$trans,
            $m,$n,
            $alpha,
            $AA,$offA,$ldA,
            $XX,$offX,$incX,
            $beta,
            $YY,$offY,$incY);
    }

    public function testGemvIllegalBufferY()
    {
        $blas = $this->getBlas();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([100,10,1]);
        $Y = $this->zeros([2]);

        [ $trans,$m,$n,$alpha,$AA,$offA,$ldA,
          $XX,$offX,$incX,$beta,$YY,$offY,$incY] =
            $this->translate_gemv($A,$X,null,null,$Y);

        $YY = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $blas->gemv(
            BLAS::RowMajor,$trans,
            $m,$n,
            $alpha,
            $AA,$offA,$ldA,
            $XX,$offX,$incX,
            $beta,
            $YY,$offY,$incY);
    }

    public function testGemvMatrixOverFlowNormal()
    {
        $blas = $this->getBlas();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([100,10,1]);
        $Y = $this->zeros([2]);

        [ $trans,$m,$n,$alpha,$AA,$offA,$ldA,
          $XX,$offX,$incX,$beta,$YY,$offY,$incY] =
            $this->translate_gemv($A,$X,null,null,$Y);

        $AA = $this->array([1,2,3,4,5])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $blas->gemv(
            BLAS::RowMajor,$trans,
            $m,$n,
            $alpha,
            $AA,$offA,$ldA,
            $XX,$offX,$incX,
            $beta,
            $YY,$offY,$incY);
    }

    public function testGemvVectorXOverFlowNormal()
    {
        $blas = $this->getBlas();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([100,10,1]);
        $Y = $this->zeros([2]);

        [ $trans,$m,$n,$alpha,$AA,$offA,$ldA,
          $XX,$offX,$incX,$beta,$YY,$offY,$incY] =
            $this->translate_gemv($A,$X,null,null,$Y);

        $XX = $this->array([10,1])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferX');
        $blas->gemv(
            BLAS::RowMajor,$trans,
            $m,$n,
            $alpha,
            $AA,$offA,$ldA,
            $XX,$offX,$incX,
            $beta,
            $YY,$offY,$incY);
    }

    public function testGemvVectorYOverFlowNormal()
    {
        $blas = $this->getBlas();

        $A = $this->array([[1,2,3],[4,5,6]]);
        $X = $this->array([100,10,1]);
        $Y = $this->zeros([2]);

        [ $trans,$m,$n,$alpha,$AA,$offA,$ldA,
          $XX,$offX,$incX,$beta,$YY,$offY,$incY] =
            $this->translate_gemv($A,$X,null,null,$Y);

        $YY = $this->array([0])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Vector specification too large for bufferY');
        $blas->gemv(
            BLAS::RowMajor,$trans,
            $m,$n,
            $alpha,
            $AA,$offA,$ldA,
            $XX,$offX,$incX,
            $beta,
            $YY,$offY,$incY);
    }

    public function testGemmNormal()
    {
        $blas = $this->getBlas();
        $A = $this->array([[1,2,3],[4,5,6],[7,8,9]]);
        $B = $this->array([[1,0,0],[0,1,0],[0,0,1]]);
        $alpha = 1.0;
        $beta  = 0.0;
        $C = $this->zeros([3,3]);
        $transA = false;
        $transB = false;

        [ $transA,$transB,$M,$N,$K,$alpha,$AA,$offA,$lda,
          $BB,$offB,$ldb,$beta,$CC,$offC,$ldc] =
            $this->translate_gemm($A,$B,$alpha,$beta,$C,$transA,$transB);

        $blas->gemm(
            BLAS::RowMajor,$transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc);

        $this->assertEquals([
            [1,2,3],
            [4,5,6],
            [7,8,9]
        ],$C->toArray());
    }

    public function testGemmTransposeSquareA()
    {
        $blas = $this->getBlas();
        $A = $this->array([[1,2,3],[4,5,6],[7,8,9]]);
        $B = $this->array([[1,0,0],[0,1,0],[0,0,1]]);
        $alpha = 1.0;
        $beta  = 0.0;
        $C = $this->zeros([3,3]);
        $transA = true;
        $transB = false;

        [ $transA,$transB,$M,$N,$K,$alpha,$AA,$offA,$lda,
          $BB,$offB,$ldb,$beta,$CC,$offC,$ldc] =
            $this->translate_gemm($A,$B,$alpha,$beta,$C,$transA,$transB);

        $blas->gemm(
            BLAS::RowMajor,$transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc);

        $this->assertEquals([
            [1,4,7],
            [2,5,8],
            [3,6,9]
        ],$C->toArray());
    }

    public function testGemmTransposeSquareB()
    {
        $blas = $this->getBlas();
        $A = $this->array([[1,0,0],[0,1,0],[0,0,1]]);
        $B = $this->array([[1,2,3],[4,5,6],[7,8,9]]);
        $alpha = 1.0;
        $beta  = 0.0;
        $C = $this->zeros([3,3]);
        $transA = false;
        $transB = true;

        [ $transA,$transB,$M,$N,$K,$alpha,$AA,$offA,$lda,
          $BB,$offB,$ldb,$beta,$CC,$offC,$ldc] =
            $this->translate_gemm($A,$B,$alpha,$beta,$C,$transA,$transB);

        $blas->gemm(
            BLAS::RowMajor,$transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc);

        $this->assertEquals([
            [1,4,7],
            [2,5,8],
            [3,6,9]
        ],$C->toArray());
    }

    public function testGemmNoTransRectangleA23()
    {
        $blas = $this->getBlas();
        $A = $this->array([[1,2,3],[4,5,6]]);
        $B = $this->array([[1,0,0],[0,1,0],[0,0,1]]);
        $alpha = 1.0;
        $beta  = 0.0;
        $C = $this->zeros([2,3]);
        $transA = false;
        $transB = false;

        [ $transA,$transB,$M,$N,$K,$alpha,$AA,$offA,$lda,
          $BB,$offB,$ldb,$beta,$CC,$offC,$ldc] =
            $this->translate_gemm($A,$B,$alpha,$beta,$C,$transA,$transB);

        $blas->gemm(
            BLAS::RowMajor,$transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc);

        $this->assertEquals([
            [1,2,3],
            [4,5,6],
        ],$C->toArray());
    }

    public function testGemmTransposeRectangleA32()
    {
        $blas = $this->getBlas();
        $A = $this->array([[1,2],[3,4],[5,6]]);
        $B = $this->array([[1,0,0],[0,1,0],[0,0,1]]);
        $alpha = 1.0;
        $beta  = 0.0;
        $C = $this->zeros([2,3]);
        $transA = true;
        $transB = false;

        [ $transA,$transB,$M,$N,$K,$alpha,$AA,$offA,$lda,
          $BB,$offB,$ldb,$beta,$CC,$offC,$ldc] =
            $this->translate_gemm($A,$B,$alpha,$beta,$C,$transA,$transB);

        $blas->gemm(
            BLAS::RowMajor,$transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc);

        $this->assertEquals([
            [1,3,5],
            [2,4,6],
        ],$C->toArray());
    }

    public function testGemmNoTransRectangleB32()
    {
        $blas = $this->getBlas();
        $A = $this->array([[1,0,0],[0,1,0],[0,0,1]]);
        $B = $this->array([[1,2],[3,4],[5,6]]);
        $alpha = 1.0;
        $beta  = 0.0;
        $C = $this->zeros([3,2]);
        $transA = false;
        $transB = false;

        [ $transA,$transB,$M,$N,$K,$alpha,$AA,$offA,$lda,
          $BB,$offB,$ldb,$beta,$CC,$offC,$ldc] =
            $this->translate_gemm($A,$B,$alpha,$beta,$C,$transA,$transB);

        $blas->gemm(
            BLAS::RowMajor,$transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc);

        $this->assertEquals([
            [1,2],
            [3,4],
            [5,6],
        ],$C->toArray());
    }

    public function testGemmTransposeRectangleB23()
    {
        $blas = $this->getBlas();
        $A = $this->array([[1,0,0],[0,1,0],[0,0,1]]);
        $B = $this->array([[1,2,3],[4,5,6]]);
        $alpha = 1.0;
        $beta  = 0.0;
        $C = $this->zeros([3,2]);
        $transA = false;
        $transB = true;

        [ $transA,$transB,$M,$N,$K,$alpha,$AA,$offA,$lda,
          $BB,$offB,$ldb,$beta,$CC,$offC,$ldc] =
            $this->translate_gemm($A,$B,$alpha,$beta,$C,$transA,$transB);

        $blas->gemm(
            BLAS::RowMajor,$transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc);

        $this->assertEquals([
            [1,4],
            [2,5],
            [3,6],
        ],$C->toArray());
    }

    public function testGemmMinusM()
    {
        $blas = $this->getBlas();
        $A = $this->array([[1,2],[3,4]]);
        $B = $this->array([[1,2],[3,4]]);
        $alpha = 1.0;
        $beta  = 0.0;
        $C = $this->zeros([2,2]);
        $transA = true;
        $transB = false;

        [ $transA,$transB,$M,$N,$K,$alpha,$AA,$offA,$lda,
          $BB,$offB,$ldb,$beta,$CC,$offC,$ldc] =
            $this->translate_gemm($A,$B,$alpha,$beta,$C,$transA,$transB);

        $M = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument m must be greater than 0.');
        $blas->gemm(
            BLAS::RowMajor,$transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc);
    }

    public function testGemmMinusN()
    {
        $blas = $this->getBlas();
        $A = $this->array([[1,2],[3,4]]);
        $B = $this->array([[1,2],[3,4]]);
        $alpha = 1.0;
        $beta  = 0.0;
        $C = $this->zeros([2,2]);
        $transA = true;
        $transB = false;

        [ $transA,$transB,$M,$N,$K,$alpha,$AA,$offA,$lda,
          $BB,$offB,$ldb,$beta,$CC,$offC,$ldc] =
            $this->translate_gemm($A,$B,$alpha,$beta,$C,$transA,$transB);

        $N = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument n must be greater than 0.');
        $blas->gemm(
            BLAS::RowMajor,$transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc);
    }

    public function testGemmMinusK()
    {
        $blas = $this->getBlas();
        $A = $this->array([[1,2],[3,4]]);
        $B = $this->array([[1,2],[3,4]]);
        $alpha = 1.0;
        $beta  = 0.0;
        $C = $this->zeros([2,2]);
        $transA = true;
        $transB = false;

        [ $transA,$transB,$M,$N,$K,$alpha,$AA,$offA,$lda,
          $BB,$offB,$ldb,$beta,$CC,$offC,$ldc] =
            $this->translate_gemm($A,$B,$alpha,$beta,$C,$transA,$transB);

        $K = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument k must be greater than 0.');
        $blas->gemm(
            BLAS::RowMajor,$transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc);
    }

    public function testGemmMinusOffsetA()
    {
        $blas = $this->getBlas();
        $A = $this->array([[1,2],[3,4]]);
        $B = $this->array([[1,2],[3,4]]);
        $alpha = 1.0;
        $beta  = 0.0;
        $C = $this->zeros([2,2]);
        $transA = true;
        $transB = false;

        [ $transA,$transB,$M,$N,$K,$alpha,$AA,$offA,$lda,
          $BB,$offB,$ldb,$beta,$CC,$offC,$ldc] =
            $this->translate_gemm($A,$B,$alpha,$beta,$C,$transA,$transB);

        $offA = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetA must be greater than or equals 0.');
        $blas->gemm(
            BLAS::RowMajor,$transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc);
    }

    public function testGemmMinusLdA()
    {
        $blas = $this->getBlas();
        $A = $this->array([[1,2],[3,4]]);
        $B = $this->array([[1,2],[3,4]]);
        $alpha = 1.0;
        $beta  = 0.0;
        $C = $this->zeros([2,2]);
        $transA = true;
        $transB = false;

        [ $transA,$transB,$M,$N,$K,$alpha,$AA,$offA,$lda,
          $BB,$offB,$ldb,$beta,$CC,$offC,$ldc] =
            $this->translate_gemm($A,$B,$alpha,$beta,$C,$transA,$transB);

        $lda = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument ldA must be greater than 0.');
        $blas->gemm(
            BLAS::RowMajor,$transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc);
    }

    public function testGemmIllegalBufferA()
    {
        $blas = $this->getBlas();
        $A = $this->array([[1,2],[3,4]]);
        $B = $this->array([[1,2],[3,4]]);
        $alpha = 1.0;
        $beta  = 0.0;
        $C = $this->zeros([2,2]);
        $transA = true;
        $transB = false;

        [ $transA,$transB,$M,$N,$K,$alpha,$AA,$offA,$lda,
          $BB,$offB,$ldb,$beta,$CC,$offC,$ldc] =
            $this->translate_gemm($A,$B,$alpha,$beta,$C,$transA,$transB);

        $AA = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $blas->gemm(
            BLAS::RowMajor,$transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc);
    }

    public function testGemmMinusOffsetB()
    {
        $blas = $this->getBlas();
        $A = $this->array([[1,2],[3,4]]);
        $B = $this->array([[1,2],[3,4]]);
        $alpha = 1.0;
        $beta  = 0.0;
        $C = $this->zeros([2,2]);
        $transA = true;
        $transB = false;

        [ $transA,$transB,$M,$N,$K,$alpha,$AA,$offA,$lda,
          $BB,$offB,$ldb,$beta,$CC,$offC,$ldc] =
            $this->translate_gemm($A,$B,$alpha,$beta,$C,$transA,$transB);

        $offB = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetB must be greater than or equals 0.');
        $blas->gemm(
            BLAS::RowMajor,$transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc);
    }

    public function testGemmMinusLdB()
    {
        $blas = $this->getBlas();
        $A = $this->array([[1,2],[3,4]]);
        $B = $this->array([[1,2],[3,4]]);
        $alpha = 1.0;
        $beta  = 0.0;
        $C = $this->zeros([2,2]);
        $transA = true;
        $transB = false;

        [ $transA,$transB,$M,$N,$K,$alpha,$AA,$offA,$lda,
          $BB,$offB,$ldb,$beta,$CC,$offC,$ldc] =
            $this->translate_gemm($A,$B,$alpha,$beta,$C,$transA,$transB);

        $ldb = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument ldB must be greater than 0.');
        $blas->gemm(
            BLAS::RowMajor,$transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc);
    }

    public function testGemmIllegalBufferB()
    {
        $blas = $this->getBlas();
        $A = $this->array([[1,2],[3,4]]);
        $B = $this->array([[1,2],[3,4]]);
        $alpha = 1.0;
        $beta  = 0.0;
        $C = $this->zeros([2,2]);
        $transA = true;
        $transB = false;

        [ $transA,$transB,$M,$N,$K,$alpha,$AA,$offA,$lda,
          $BB,$offB,$ldb,$beta,$CC,$offC,$ldc] =
            $this->translate_gemm($A,$B,$alpha,$beta,$C,$transA,$transB);

        $BB = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $blas->gemm(
            BLAS::RowMajor,$transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc);
    }

    public function testGemmMinusOffsetC()
    {
        $blas = $this->getBlas();
        $A = $this->array([[1,2],[3,4]]);
        $B = $this->array([[1,2],[3,4]]);
        $alpha = 1.0;
        $beta  = 0.0;
        $C = $this->zeros([2,2]);
        $transA = true;
        $transB = false;

        [ $transA,$transB,$M,$N,$K,$alpha,$AA,$offA,$lda,
          $BB,$offB,$ldb,$beta,$CC,$offC,$ldc] =
            $this->translate_gemm($A,$B,$alpha,$beta,$C,$transA,$transB);

        $offC = -1;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument offsetC must be greater than or equals 0.');
        $blas->gemm(
            BLAS::RowMajor,$transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc);
    }

    public function testGemmMinusLdC()
    {
        $blas = $this->getBlas();
        $A = $this->array([[1,2],[3,4]]);
        $B = $this->array([[1,2],[3,4]]);
        $alpha = 1.0;
        $beta  = 0.0;
        $C = $this->zeros([2,2]);
        $transA = true;
        $transB = false;

        [ $transA,$transB,$M,$N,$K,$alpha,$AA,$offA,$lda,
          $BB,$offB,$ldb,$beta,$CC,$offC,$ldc] =
            $this->translate_gemm($A,$B,$alpha,$beta,$C,$transA,$transB);

        $ldc = 0;
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Argument ldC must be greater than 0.');
        $blas->gemm(
            BLAS::RowMajor,$transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc);
    }

    public function testGemmIllegalBufferC()
    {
        $blas = $this->getBlas();
        $A = $this->array([[1,2],[3,4]]);
        $B = $this->array([[1,2],[3,4]]);
        $alpha = 1.0;
        $beta  = 0.0;
        $C = $this->zeros([2,2]);
        $transA = true;
        $transB = false;

        [ $transA,$transB,$M,$N,$K,$alpha,$AA,$offA,$lda,
          $BB,$offB,$ldb,$beta,$CC,$offC,$ldc] =
            $this->translate_gemm($A,$B,$alpha,$beta,$C,$transA,$transB);

        $CC = new \stdClass();
        $this->expectException(TypeError::class);
        $this->expectExceptionMessage('must implement interface Interop\Polite\Math\Matrix\LinearBuffer');
        $blas->gemm(
            BLAS::RowMajor,$transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc);
    }

    public function testGemmMatrixAOverFlowTransposeRectangleA32()
    {
        $blas = $this->getBlas();
        $A = $this->array([[1,2],[3,4],[5,6]]);
        $B = $this->array([[1,0,0],[0,1,0],[0,0,1]]);
        $alpha = 1.0;
        $beta  = 0.0;
        $C = $this->zeros([2,3]);
        $transA = true;
        $transB = false;

        [ $transA,$transB,$M,$N,$K,$alpha,$AA,$offA,$lda,
          $BB,$offB,$ldb,$beta,$CC,$offC,$ldc] =
            $this->translate_gemm($A,$B,$alpha,$beta,$C,$transA,$transB);

        $AA = $this->array([1,2,3,4,5])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $blas->gemm(
            BLAS::RowMajor,$transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc);
    }

    public function testGemmMatrixBOverFlowTransposeRectangleA32()
    {
        $blas = $this->getBlas();
        $A = $this->array([[1,2],[3,4],[5,6]]);
        $B = $this->array([[1,0,0],[0,1,0],[0,0,1]]);
        $alpha = 1.0;
        $beta  = 0.0;
        $C = $this->zeros([2,3]);
        $transA = true;
        $transB = false;

        [ $transA,$transB,$M,$N,$K,$alpha,$AA,$offA,$lda,
          $BB,$offB,$ldb,$beta,$CC,$offC,$ldc] =
            $this->translate_gemm($A,$B,$alpha,$beta,$C,$transA,$transB);

        $BB = $this->array([1,0,0, 0,1,0, 0,0])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferB');
        $blas->gemm(
            BLAS::RowMajor,$transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc);
    }

    public function testGemmOutputOverFlowTransposeRectangleA32()
    {
        $blas = $this->getBlas();
        $A = $this->array([[1,2],[3,4],[5,6]]);
        $B = $this->array([[1,0,0],[0,1,0],[0,0,1]]);
        $alpha = 1.0;
        $beta  = 0.0;
        $C = $this->zeros([2,3]);
        $transA = true;
        $transB = false;

        [ $transA,$transB,$M,$N,$K,$alpha,$AA,$offA,$lda,
          $BB,$offB,$ldb,$beta,$CC,$offC,$ldc] =
            $this->translate_gemm($A,$B,$alpha,$beta,$C,$transA,$transB);


        $CC = $this->zeros([5])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferC');
        $blas->gemm(
            BLAS::RowMajor,$transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc);
    }

    public function testGemmMatrixAOverFlowTransposeRectangleB23()
    {
        $blas = $this->getBlas();
        $A = $this->array([[1,0,0],[0,1,0],[0,0,1]]);
        $B = $this->array([[1,2,3],[4,5,6]]);
        $alpha = 1.0;
        $beta  = 0.0;
        $C = $this->zeros([3,2]);
        $transA = false;
        $transB = true;

        [ $transA,$transB,$M,$N,$K,$alpha,$AA,$offA,$lda,
          $BB,$offB,$ldb,$beta,$CC,$offC,$ldc] =
            $this->translate_gemm($A,$B,$alpha,$beta,$C,$transA,$transB);

        $AA = $this->array([1,0,0, 0,1,0, 0,0])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferA');
        $blas->gemm(
            BLAS::RowMajor,$transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc);
    }

    public function testGemmMatrixBOverFlowTransposeRectangleB23()
    {
        $blas = $this->getBlas();
        $A = $this->array([[1,0,0],[0,1,0],[0,0,1]]);
        $B = $this->array([[1,2,3],[4,5,6]]);
        $alpha = 1.0;
        $beta  = 0.0;
        $C = $this->zeros([3,2]);
        $transA = false;
        $transB = true;

        [ $transA,$transB,$M,$N,$K,$alpha,$AA,$offA,$lda,
          $BB,$offB,$ldb,$beta,$CC,$offC,$ldc] =
            $this->translate_gemm($A,$B,$alpha,$beta,$C,$transA,$transB);

        $BB = $this->array([1,2,3,4,5])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferB');
        $blas->gemm(
            BLAS::RowMajor,$transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc);
    }

    public function testGemmOutputOverFlowTransposeRectangleB23()
    {
        $blas = $this->getBlas();
        $A = $this->array([[1,0,0],[0,1,0],[0,0,1]]);
        $B = $this->array([[1,2,3],[4,5,6]]);
        $alpha = 1.0;
        $beta  = 0.0;
        $C = $this->zeros([3,2]);
        $transA = false;
        $transB = true;

        [ $transA,$transB,$M,$N,$K,$alpha,$AA,$offA,$lda,
          $BB,$offB,$ldb,$beta,$CC,$offC,$ldc] =
            $this->translate_gemm($A,$B,$alpha,$beta,$C,$transA,$transB);

        $CC = $this->zeros([5])->buffer();
        $this->expectException(InvalidArgumentException::class);
        $this->expectExceptionMessage('Matrix specification too large for bufferC');
        $blas->gemm(
            BLAS::RowMajor,$transA,$transB,
            $M,$N,$K,
            $alpha,
            $AA,$offA,$lda,
            $BB,$offB,$ldb,
            $beta,
            $CC,$offC,$ldc);
    }
}
