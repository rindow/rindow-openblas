IMPORTANT
=========
Development of this program has ended.
Please migrate to the program below instead.

- https://github.com/rindow/rindow-openblas-ffi

We stopped using PHP extensions because it was too difficult to prepare binary files for each PHP version and Linux version.

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

- PHP7.2 or PHP7.3 or PHP7.4 or PHP8.0 or PHP8.1 or PHP8.2 or PHP8.3
- Linux or Windows 10
- [OpenBLAS](https://github.com/OpenMathLib/OpenBLAS)
- [Rindow-Matlib](https://github.com/rindow/rindow-matlib)

Recommend environment
=====================

- PHP8.1 or PHP8.2 or PHP8.3
- OpenCL binding for PHP - rindow-opencl 0.2.0. [sources](https://github.com/rindow/rindow-opencl), [binaries](https://github.com/rindow-opencl/releases)
- BLAS libray for OpenCL implements - rindow-clblast 0.2.0. [sources](https://github.com/rindow/rindow-clblast), [binaries](https://github.com/rindow-clblast/releases)
- Matrix PHP library - rindow-math-matrix 2.0.0 or later. [sources](https://github.com/rindow/rindow-math-matrix)
- Driver Pack - rindow-math-matrix-matlibext 1.0.0 or later. [sources](https://github.com/rindow/rindow-math-matrix-matlibext)

How to build from source code on Linux
======================================
You can also build and use from source code.
(Please change the notation of the php version number according to your environment.)

Install build tools and OpenBLAS libray
---------------------------------------
Install gcc development environment and openblas library.
Then install the php development environment according to the target php version.

```shell
$ sudo apt install build-essential autoconf automake libtool bison re2c
$ sudo apt install pkg-config
$ sudo apt install libopenblas-dev
$ sudo apt install liblapacke-dev
$ sudo apt install phpX.X-dev (ex. php8.1-dev)
```
 If you want to use the latest version of openblas, download the source code from [the site](https://github.com/xianyi/OpenBLAS/releases), build it, and set the installation location of openblas in PKG_CONFIG_PATH

### Build
Run the target php version of phpize and build.

```shell
$ cd /path/to/build/directory
$ git clone https://github.com/rindow/rindow-matlib
$ cd rindow-matlib
$ cmake -S . -B build
$ cd ..
$ git clone https://github.com/rindow/rindow-openblas
$ cd rindow-openblas
$ composer update
$ phpizeX.X (ex. phpize8.1)
$ ./configure --enable-rindow_openblas --with-rindow_matlib=/path/to/rindow-matlib --with-php-config=php-configX.X (ex. php-config8.1)
$ make clean
$ make
$ make test
```

### Install from built directory

```shell
$ sudo make install
```
Add the "extension=rindow_openblas" entry to php.ini or Make the file rindow_openblas.ini.

If you want an easier install, use the following spell instead of "make install" and creating an ini file.

```shell
$ sh ./packaging.sh X.X   (ex. sh ./packaging.sh 8.1)
$ sudo apt install ./rindow-openblas-phpX.X_X.X.X_amd64.deb (ex. ./rindow-openblas-php8.1_...)
```

How to build from source code on Windows
========================================
You can also build and use from source code.


Download or Build OpenBLAS for MSVC on Windows
----------------------------------------------
### Download binaries for the OpenBLAS libray
You need to build OpenBLAS libray for MSVC or download built binaries of libray for MSVC.

If you want to use the pre-built OpenBLAS libray, you need OpenBLAS release 0.3.10 or later.

- https://github.com/xianyi/OpenBLAS/releases

### Install VC15 or VC16
Developing PHP extensions for php7.x requires VC15 instead of the latest VC.

- Install Microsoft Visual Studio 2019 or later installer
- Run Installer with vs2017 build tools option.

Developing PHP extensions from php8.x, requires VS16. You can use Visual Studio 2019.

### Build OpenBLAS for pure MSVC from sources
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

- You must know that PHP 7.x needs environment for the MSVC version vc15 (that means Visual Studio 2017). php-sdk releases 2.1.9 supports vc15.
- For PHP 7.x, Download the php-sdk from https://github.com/microsoft/php-sdk-binary-tools/releases/tag/php-sdk-2.1.9
- If you want to build extensions for PHP 8.x, You have to use php-sdk release 2.2.0. It supports vs16.
- For PHP 8.x, Download the php-sdk from https://github.com/microsoft/php-sdk-binary-tools/releases/tag/php-sdk-2.2.0
- Extract to c:\php-sdk
- Download target Development package from https://windows.php.net/download
- Extract to /path/to/php-devel-pack-x.x.x-Win32-Vxxx-x64/

### start php-sdk for target PHP version

Open Visual Studio Command Prompt for VS for the target PHP version(see stepbystepbuild.)
Note that you must explicitly specify the version of vc15 for which php.exe was built.
The -vcvars_ver=14.16 means vc15.

If you want to build for PHP 8.x, No options required.

```shell
C:\visual\studio\path>vcvars64 -vcvars_ver=14.16
or
C:\visual\studio\path>vcvars64

C:\tmp>cd c:\php-sdk-x.x.x

C:\php-sdk-2.1.9>phpsdk-vc15-x64.bat
or
C:\php-sdk-2.2.0>phpsdk-vs16-x64.bat

```

### Build

```shell
$ cd /path/to/Rindow-Matlib-project-directory
$ cmake -S . -B build
$ PATH %PATH%;/path/to/OpenBLAS/bin
$ cd /path/to/here
$ composer update
$ /path/to/php-devel-pack-x.x.x-Win32-VXXX-x64/phpize.bat
$ configure --enable-rindow_openblas --with-prefix=/path/to/php-installation-path --with-openblas=/path/to/OpenBLAS-libray-built-directory --with-matlib=/path/to/Rindow-Matlib-project-directory
```
Edit "#define PHP_BUILD_SYSTEM" line in the "/php-devel-pack-xxx/include/main/config.w32.h"
Change the PHP_BUILD_SYSTEM definition to the same value as in "/php-devel-pack-xxx/include/main/config.pickle.h". If the values are not the same, a warning error will occur during build.

And then Build.

```shell
$ nmake clean
$ nmake
$ nmake test
```

### Install from built directory

- Copy the php extension binary(.dll) to the php/ext directory from here/arch/Releases_TS/php_rindow_openblas.dll
- Add the "extension=php_rindow_openblas" entry to php.ini
