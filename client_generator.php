<?php

class OpenBLASClientGenerator
{
    public array $excludeFuncs = [
        'openblas_setaffinity',
        'openblas_getaffinity',
        'cblas_csrot',
        'cblas_zdrot',
        'cblas_crotg',
        'cblas_zrotg',
        'cblas_xerbla',
        'cblas_sbstobf16',
        'cblas_sbdtobf16',
        'cblas_sbf16tos',
        'cblas_dbf16tod',
        'cblas_sbdot',
        'cblas_sbgemv',
        'cblas_sbgemm',
    ];

    public array $funcDecrWord = [
        'int','char*','void','float','double','CBLAS_INDEX',
        //'openblas_complex_float','openblas_complex_double',
        'lapack_int',
        //'lapack_complex_float','lapack_complex_double',
    ];

    public string $preg_types;

    public function __construct()
    {
        $this->preg_types = implode('|',array_map(fn($x) => str_replace('*','\\*',$x),$this->funcDecrWord));
    }

    public function generator($argv)
    {
        $dummy = array_shift($argv);
        $outputfile = array_shift($argv);
        $fpi = null;
        $fpo = null;
        $funcs = [];
        try {
            $fpo = fopen($outputfile,'w');
            if($fpo==null) {
                throw new RuntimeException("Error opening input file: $outputfile");
            }
            $code = $this->beginTemplate();
            if(!fwrite($fpo,$code)) {
                throw new RuntimeException("fwrite error");
            }
            while($inputfile=array_shift($argv)) {
                $fpi = fopen($inputfile,'r');
                if($fpi==null) {
                    throw new RuntimeException("Error opening input file: $inputfile");
                }
                $funcs = array_merge($funcs,$this->generate($fpo,$fpi));
            }
            $code = $this->endTemplate($funcs);
            if(!fwrite($fpo,$code)) {
                throw new RuntimeException("fwrite error");
            }
        } finally {
            if($fpi) {
                fclose($fpi);
            }
            if($fpo) {
                fclose($fpo);
            }
        }
    }
    
    public function generate($fpo,$fpi) : array
    {
        $excludeFuncs = $this->excludeFuncs;
        $eof = false;
        $funcs = [];
        while(!$eof) {
            $line='';
            while(true) {
                $next=fgets($fpi);
                if(!$next) {
                    $eof = true;
                    break;
                }
                $next = trim($next);
                if(substr($next,0,1)=='#') {
                    break;
                }
                if(substr($next,0,2)=='/*') {
                    break;
                }
                if(substr($next,0,2)=='//') {
                    break;
                }
                $line .= trim($next);
                if(substr($next,-1,1)==';') {
                    break;
                }
            }
            $declare = $this->parser($line);
            //var_dump($line);
            //var_dump($declare);
            //echo "PAUSE>";
            //fgets(STDIN);
            if($declare==null) {
                continue;
            }
            if(in_array($declare['func'],$excludeFuncs)) {
                continue;
            }
            $code = $this->funcTemplate($declare);
            if(!fwrite($fpo,$code)) {
                throw new RuntimeException("fwrite error");
            }
            $funcs[] = $declare;
        }
        return $funcs;
    }
    
    public function parser(string $line) : ?array
    {
        $funcDecrWord = $this->funcDecrWord;
        $tmp = explode(' ',$line);
        $head = $tmp[0];
        if(!in_array($head,$funcDecrWord)) {
            return null;
        }
        $pattern = "/^(".$this->preg_types.") *([A-Za-z0-9_]+) *\\(([^)]+)/";
        //var_dump($pattern);
        preg_match($pattern,$line,$match);
        if(array_key_exists(1,$match)) {
            $return = $match[1];
        } else {
            var_dump("unmatch return $line");
            $return = "unknown";
        }
        if(array_key_exists(2,$match)) {
            $func = $match[2];
        } else {
            var_dump("unmatch func $line");
            $func = "unknown";
        }
        if(array_key_exists(3,$match)) {
            $args = $match[3];
        } else {
            var_dump("unmatch args $line");
            $args = ["unknown","unknown"];
        }
        $args = explode(',',$args);
        $strargs = array_map('trim',$args);
        $args = [];
        foreach($strargs as $arg) {
            $var = null;
            $arg = str_replace('*',' *',$arg);
            $types2 = explode(' ',$arg);
            $types = [];
            foreach($types2 as $t) {
                if($t) {
                    $types[] = $t;
                }
            }
            $types = array_map('trim',$types);
            if(count($types)>1) {
                $var = array_pop($types);
                $var = trim($var);
                if(substr($var,0,1)=='*') {
                    $var = substr($var,1);
                    array_push($types,'*');
                }
            }
            $type = implode(' ',$types);
            $args[] = ['type'=>$type,'var'=>$var];
        }
        //var_dump("$return $func");
        //echo "args:";
        //var_dump($args);
        //if($func=='cblas_xerbla') {
        //    var_dump(['return' => $return, 'func'=>$func, 'args'=>$args]);
        //}
        
        return ['return' => $return, 'func'=>$func, 'args'=>$args];
    }
    public function funcTemplate($declare)
    {
        $return = $declare['return'];
        $funcname = $declare['func'];
        // typedef 
        $code = "typedef {$return} (CALLBACK* PFN{$funcname})( /* {$funcname} */";
        $isNext = false;
        foreach($declare['args'] as $arg) {
            $type = $arg['type'];
            $var = $arg['var'];
            if($isNext) {
                $code .= ",";
            }
            $code .= "\n";
            $code .= "    {$type}            /* {$var} */";
            $isNext = true;
        }
        $code .= "\n";
        $code .= ");\n";
    
        // static function pointer
        $code .= "static PFN{$funcname} _g_{$funcname} = NULL;\n";
    
        // proxy function
        $code .= "{$return} {$funcname}(";
        $isNext = false;
        foreach($declare['args'] as $arg) {
            $type = $arg['type'];
            $var = $arg['var'];
            if($isNext) {
                $code .= ",";
            }
            $code .= "\n";
            $code .= "    {$type}            {$var}";
            $isNext = true;
        }
        $code .= "\n";
        $code .= ")\n";
        $code .= "{\n";
        $code .= "    if(_h_openblas==NULL || _g_{$funcname}==NULL) {\n";
        if($return=='void') {
            $code .= "        return;\n";
        } else {
            $code .= "        return 0;\n";
        }
        $code .= "    }\n";
        if($return=='void') {
            $code .= "    _g_{$funcname}(";
        } else {
            $code .= "    return _g_{$funcname}(";
        }
        $isNext = false;
        foreach($declare['args'] as $arg) {
            $type = $arg['type'];
            $var = $arg['var'];
            if($isNext) {
                $code .= ",";
            }
            $code .= "\n";
            if($var) {
                $code .= "        {$var}";
            }
            $isNext = true;
        }
        $code .= "    \n";
        $code .= "    );\n";
        $code .= "}\n";
        return $code;
    }
    
    function beginTemplate()
    {
        $code  = "#include <stdio.h>\n";
        $code .= "#include <Windows.h>\n";
        $code .= "#include <cblas.h>\n";
        $code .= "\n";
        $code .= "#define LOADFUNC(funcname) \\\n";
        $code .= "_g_##funcname = (PFN##funcname)GetProcAddress( _h_openblas, #funcname ); \\\n";
        $code .= "if(_g_##funcname==NULL) { \\\n";
        $code .= "    printf(\"load error: %s\\n\",  #funcname); \\\n";
        $code .= "    return -1; \\\n";
        $code .= "} \\\n";
        $code .= "\n";
        $code .= "static HMODULE _h_openblas = NULL;\n";
        return $code;
    }
    
    public function endTemplate($funcs)
    {
        $code  = "int rindow_load_openblas_dll()\n";
        $code .= "{\n";
        $code .= "    if(_h_openblas!=NULL) {\n";
        $code .= "        return 0;\n";
        $code .= "    }\n";
        $code .= "    _h_openblas = LoadLibraryA( \"libopenblas.dll\" );\n";
        $code .= "    if(_h_openblas==NULL) {\n";
        $code .= "        printf(\"load error: libopenblas\\n\");\n";
        $code .= "        return -1;\n";
        $code .= "    }\n";
        foreach($funcs as $declare) {
            $funcname = $declare['func'];
            $code .= "    LOADFUNC({$funcname})\n";
        }
        $code .= "    return 0;\n";
        $code .= "}\n";
        $code .= "void rindow_unload_openblas_dll()\n";
        $code .= "{\n";
        $code .= "    FreeLibrary( _h_openblas );\n";
        $code .= "    _h_openblas = NULL;\n";
        $code .= "}\n";
        return $code;
    }
}

$generator = new OpenBLASClientGenerator();
$generator->generator($argv);
