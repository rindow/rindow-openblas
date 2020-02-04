<?php
require __DIR__.'/../phpunit/vendor/autoload.php';
if(!extension_loaded('rindow_openblas')) {
    echo "rindow_openblas extention is not loaded.\n";
    exit;
}
use Interop\Polite\Math\Matrix\NDArray;
use Interop\Polite\Math\Matrix\BLAS;
$mo = new Rindow\Math\Matrix\MatrixOperator();
$blas = new Rindow\OpenBlas\Blas();


$a = $mo->array([1,2],NDArray::float32);

$blas->scal(2,10,$a->buffer(),0,1);
var_dump($a[0]);
var_dump($a[1]);

$a = $mo->array([1,2],NDArray::float32);
$b = $mo->array([50,60],NDArray::float32);

$blas->scal($a->size(),2.0,$a->buffer(),0,1);
var_dump($a->toArray());

$blas->axpy($a->size(),2.0,$a->buffer(),0,1,$b->buffer(),0,1);
var_dump($b->toArray());

$a = $mo->array([1,2,3,4,5,6],NDArray::float32);
$b = $mo->array([3,4,5,6,7,8],NDArray::float32);
var_dump($blas->dot($a->size(),$a->buffer(),0,1,$b->buffer(),0,1));

$a = $mo->array([-1,2,3,4,5,6],NDArray::float32);
var_dump($blas->asum($a->size(),$a->buffer(),0,1));

$a = $mo->array([-1,2,3,4,5,-6],NDArray::float32);
var_dump($blas->iamax($a->size(),$a->buffer(),0,1));

$a = $mo->array([-1,2,3,4,-5,6],NDArray::float32);
var_dump($blas->iamin($a->size(),$a->buffer(),0,1));

$a = $mo->array([-1,2,3,4,-5,6],NDArray::float32);
$b = $mo->zerosLike($a);
$blas->copy($a->size(),$a->buffer(),0,1,$b->buffer(),0,1);
var_dump($b->toArray());

// Use BLAS
// Use BLAS
// 3x3 * 3
//    y := alpha * Ax + beta * y
$A = $mo->array([[1,2,3],[4,5,6],[7,8,9]],NDArray::float32);
$X = $mo->array([2,3,4],NDArray::float32);
$Y = $mo->zeros([3],NDArray::float32);
// [ 20,  47,  74],
$blas->gemv(BLAS::RowMajor,BLAS::NoTrans,$Y->size(),$X->size(),
    1.0,$A->buffer(),0,$Y->size(),$X->buffer(),0,1,0.0,$Y->buffer(),0,1);
echo 'gemv=';
var_dump($Y->toArray());

// Use BLAS
// 3x3 * 3x3
//    C := alpha * AB + beta * C
$A = $mo->array([[1,2,3],[4,5,6],[7,8,9]],NDArray::float32);
$B = $mo->array([[2,3,4],[5,6,7],[8,9,1]],NDArray::float32);
$C = $mo->zeros([3,3],NDArray::float32);
//    [[ 36,  42,  21],
//     [ 81,  96,  57],
//     [126, 150,  93]],
//var_dump($A->shape());
//var_dump($B->shape());
//var_dump($C->shape());
$blas->gemm(
    BLAS::RowMajor,BLAS::NoTrans,BLAS::NoTrans,
    3,3,3,
    1.0,
    $A->buffer(),0,3,
    $B->buffer(),0,3,
    0.0,
    $C->buffer(),0,3);
echo 'gemm=';
var_dump($C->toArray());
