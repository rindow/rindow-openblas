<?php
require __DIR__.'/../phpunit/vendor/autoload.php';
if(!extension_loaded('rindow_openblas')) {
    echo "rindow_openblas extention is not loaded.\n";
    exit;
}
use Interop\Polite\Math\Matrix\NDArray;

$x = new Rindow\OpenBlas\Buffer(2,NDArray::float32);

$x[0] = 10;
$x[1] = 20;

assert(is_float($x[0]),'float type');
assert($x[0]==10.0,'set and get value');
assert($x[1]==20.0,'set and get value');
assert(count($x)==2,'Countable buffer size');

$fn=function($x,$i){$t=false;try{
    $x[$i]=1;
}catch(RuntimeException $e){$t=true;}return $t;};
assert($fn($x,2), 'catch buffer overflow');
assert($fn($x,-1),'catch buffer overflow');
$fn=function($x,$i){$t=false;try{
    echo $x[$i];
}catch(RuntimeException $e){$t=true;}return $t;};
assert($fn($x,2), 'catch buffer overflow');
assert($fn($x,-1),'catch buffer overflow');

assert(isset($x[0])==true,'exists');
assert(isset($x[1])==true,'exists');
assert(isset($x[2])==false,'not exists');
assert(isset($x[-1])==false,'not exists');
unset($x[0]);
assert($x[0]==0.0,'unset means set zero');
$x[0] = 10;

assert($x->dtype()==NDArray::float32,'dtype');

$binary = $x->dump();
assert(strlen($binary)==32/8*2,'dump binary');
$y = new Rindow\OpenBlas\Buffer(2,NDArray::float32);
$y->load($binary);
assert($y[0]==10.0,'load binary');
assert($y[1]==20.0,'load binary');
