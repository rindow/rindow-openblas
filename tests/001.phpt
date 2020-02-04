--TEST--
Check if rindow_openblas is loaded
--SKIPIF--
<?php
if (!extension_loaded('rindow_openblas')) {
	echo 'skip';
}
?>
--FILE--
<?php
echo 'The extension "rindow_openblas" is available';
?>
--EXPECT--
The extension "rindow_openblas" is available
