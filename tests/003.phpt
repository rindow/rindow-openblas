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
define('float32',12);
$buf = new Rindow\OpenBLAS\Buffer(4,float32);
$buf[0] = 1;
$buf[1] = 2;
$buf[2] = 3;
$buf[3] = 4;
$blas = new Rindow\OpenBLAS\Blas();
echo $blas->asum(4,$buf,0,1);
?>
--EXPECT--
10
