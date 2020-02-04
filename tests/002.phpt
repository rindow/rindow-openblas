--TEST--
Buffer
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
echo count($buf);
for($i=0;$i<4;$i++) {
    echo $buf[$i];
}
?>
--EXPECT--
41234
