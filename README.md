Rindow OpenBLAS PHP extension
=============================
The Rindow OpenBLAS PHP extension is universal Buffer for N-dimension and OpenBLAS and Mathematical library.

- Provides Universal Buffer for 1-dimension for data exchange between C,C+ language and PHP.
- The OpenBLAS library available to PHP. Only the commonly used functions in OpenBLAS are provided.
- Provides commonly used Mathematical libraries not included in OpenBLAS.

You can do very fast N-dimensional array operations in conjunction with the [Rindow Math Matrix](https://github.com/rindow/rindow-math-matrix).

Very useful when you want to do deep learning with PHP!

Requirements
============

- PHP7.2 or PHP7.3 or PHP7.4
- Linux or Windows 10

Download binaries
=================
You can download and use pre-built Windows binaries.
Download the binary for your version of PHP.

- https://github.com/rindow/rindow-openblas-binaries/

Please download the following two binaries.

- The PHP extension of rindow-openblas that matches the php version.
- The flang runtime DLLs.

If you are already using miniconda3, please use the flang runtime DLLs in it.

How to setup
============
Copy the shared library to the PHP extension directory and set it in php.ini.

And then Set the execution path to the unzipped flang runtime DLLs.

How to build from source code on Linux
======================================
You can also build and use from source code.

Install build tools and OpenBLAS libray
---------------------------------------
Install gcc development environment and openblas library.
Then install the php development environment according to the target php version.

```shell
$ sudo apt install build-essential autoconf automake libtool bison re2c
$ sudo apt install pkg-config
$ sudo apt install libopenblas-dev
$ sudo apt install liblapacke-dev
$ sudo apt install php7.4-dev
```
 If you want to use the latest version of openblas, download the source code from [the site](https://github.com/xianyi/OpenBLAS/releases), build it, and set the installation location of openblas in PKG_CONFIG_PATH

### Build
Run the target php version of phpize and build.

```shell
$ git clone https://github.com/rindow/rindow-openblas
$ cd rindow_openblas
$ phpize7.4
$ ./configure --enable-rindow_openblas --with-php-config=php-config7.4
$ make clean
$ make
$ make test
```

### Install from built directory

```shell
$ sudo make install
```
Add the "extension=rindow_openblas.so" entry to php.ini



How to build from source code on Windows
========================================
You can also build and use from source code.


Build OpenBLAS for MSVC 15 on Windows
-------------------------------------

### Download binaries for the OpenBLAS libray
You need to build OpenBLAS libray for MSVC or download built binaries of libray for MSVC.
You can download the binaries from us.
This binary can be used in pure msvc only environment unlike the binary on openblas site.

- [Rindow OpenBLAS Development Kit for Windows](https://github.com/rindow/rindow-openblas-binaries/tree/master/devel/windows)

### Install VC15
Developing PHP extensions from php7.2 to php7.4 requires VC15 instead of the latest VC.

- Install Microsoft Visual Studio 2019 or later installer
- Run Installer with vs2017 build tools option.

### Build OpenBLAS for pure MSVC
If you want to build the OpenBLAS on MSVC with static library instead you use pre-build binary on our site, you can build it yourself.

https://github.com/xianyi/OpenBLAS/wiki/How-to-use-OpenBLAS-in-Microsoft-Visual-Studio
> 1. Native (MSVC) ABI
> Install Miniconda3 for 64 bits. And then follow the procedure described on the above page.

You want to build a DLL of OpenBLAS, you can run cmake with shared libray option -DBUILD_SHARED_LIBS=ON

```shell
Anaconda3>vcvars64 -vcvars_ver=14.16

...... omission

Anaconda3>cmake .. -G "Ninja" -DCMAKE_CXX_COMPILER=clang-cl -DCMAKE_C_COMPILER=clang-cl -DCMAKE_Fortran_COMPILER=flang -DBUILD_WITHOUT_LAPACK=no -DNOFORTRAN=0 -DDYNAMIC_ARCH=ON -DCMAKE_BUILD_TYPE=Release
Anaconda3>cmake --build . --config Release
```


Build the extension for Windows
-------------------------------

Please refer to the following URL for details.

https://wiki.php.net/internals/windows/stepbystepbuild_sdk_2

### php sdk and devel-pack binaries for windows

- You must know that PHP 7.2,7.3 and 7.4 needs environment for the MSVC version vc15 (that means Visual Studio 2017). php-sdk releases 2.1.9 supports vc15.
- Download https://github.com/microsoft/php-sdk-binary-tools/releases/tag/php-sdk-2.1.9
- Extract to c:\php-sdk
- Download target dev-pack from https://windows.php.net/downloads/releases/
- Extract to /path/to/php-devel-pack-7.x.x-Win32-VC15-x64/

### start php-sdk for target PHP version

Open Visual Studio Command Prompt for VS for the target PHP version(see stepbystepbuild.)
Note that you must explicitly specify the version of vc15 for which php.exe was built.
The -vcvars_ver=14.16 means vc15.

```shell
C:\visual\studio\path>vcvars64 -vcvars_ver=14.16

C:\tmp>cd c:\php-sdk
C:\php-sdk>phpsdk-vc15-x64.bat
```

### Build

```shell
$ cd /path/to/here
$ /path/to/php-devel-pack-7.x.x-Win32-VC15-x64/phpize.bat
$ configure --enable-rindow_openblas --with-prefix=/path/to/php-installation-path --with-openblas=/path/to/OpenBLAS-libray-built-directory --with-flang=/path/to/miniconda3-directory/Library
$ nmake clean
$ nmake
$ nmake test
```

### Install from built directory

- Copy the php extension binary(.dll) to the php/ext directory from here/arch/Releases_TS/php_rindow_openblas.dll
- Add the "extension=" entry to php.ini
