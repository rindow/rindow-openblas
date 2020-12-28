--TEST--
Blas
--SKIPIF--
<?php
if (!extension_loaded('rindow_openblas')) {
	echo 'skip';
}
?>
--FILE--
<?php
require __DIR__.'/../vendor/autoload.php';
include __DIR__.'/../testPHP/HostBuffer.php';
use Interop\Polite\Math\Matrix\NDArray;
use Interop\Polite\Math\Matrix\BLAS;
use RindowTest\OpenBlas\HostBuffer as Buffer;
$buf = new Buffer(4,NDArray::float32);
$buf[0] = 1;
$buf[1] = 2;
$buf[2] = 3;
$buf[3] = 4;
$blas = new Rindow\OpenBLAS\Blas();
echo $blas->asum(4,$buf,0,1);
?>
--EXPECT--
10
