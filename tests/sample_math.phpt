--TEST--
Sample BLAS
--SKIPIF--
<?php
if (!extension_loaded('rindow_openblas')) {
	echo 'skip';
}
?>
--FILE--
<?php
require __DIR__.'/../vendor/autoload.php';
if(!extension_loaded('rindow_openblas')) {
    echo "rindow_openblas extention is not loaded.\n";
    exit;
}
use Interop\Polite\Math\Matrix\NDArray;
use Interop\Polite\Math\Matrix\BLAS;
$mo = new Rindow\Math\Matrix\MatrixOperator();
$math = new Rindow\OpenBlas\Math();

$a = $mo->array([1,2],NDArray::float32);
$math->increment($a->size(),$alpha=2.0,$a->buffer(),0,1,$beta=1.0);
var_dump($a[0]);
var_dump($a[1]);
?>
--EXPECT--
float(3)
float(5)