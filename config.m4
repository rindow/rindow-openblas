dnl config.m4 for extension rindow_openblas

PHP_ARG_ENABLE(rindow_openblas, whether to enable rindow_openblas support,
dnl Make sure that the comment is aligned:
[  --enable-rindow_openblas          Enable rindow_openblas support], no)
PHP_ARG_WITH(rindow_matlib, rindow_matlib support,
[  --with-rindow_matlib=DIR          Specify rindow_matlib path])

if test "$PHP_RINDOW_OPENBLAS" != "no"; then

  dnl # get library OpenBLAS build options from pkg-config output
  AC_PATH_PROG(PKG_CONFIG, pkg-config, no)
  AC_MSG_CHECKING(for openblas)
  if test -x "$PKG_CONFIG"; then
    if $PKG_CONFIG --exists openblas; then
      if $PKG_CONFIG --exists lapacke; then
        if $PKG_CONFIG openblas --atleast-version 0.2.0; then
          LIBBLAS_CFLAGS="`$PKG_CONFIG openblas --cflags` `$PKG_CONFIG lapacke --libs`"
          LIBBLAS_LIBDIR="`$PKG_CONFIG openblas --libs` `$PKG_CONFIG lapacke --libs`"
          LIBBLAS_VERSON=`$PKG_CONFIG openblas --modversion`
          AC_MSG_RESULT(from pkgconfig: version $LIBBLAS_VERSON)
          if $PKG_CONFIG openblas --atleast-version 0.3.6; then
              AC_DEFINE(OPENBLAS_HAVE_IAMIN, 1, [openblas have iamin])
          fi
        else
          AC_MSG_ERROR(system openblas is too old: version 0.2.0 required)
        fi
      else
        AC_MSG_ERROR(lapacke not found)
      fi
    else
      AC_MSG_ERROR(openblas not found)
    fi
  else
    AC_MSG_ERROR(pkg-config not found)
  fi
  LIBMATLIB_LIBDIR=PHP_EXT_SRCDIR(rindow_openblas)[/lib]
  dnl # LIBBLAS_LIBDIR="-L$LIBMATLIB_LIBDIR -lrindowmatlib $LIBBLAS_LIBDIR"
  PHP_EVAL_LIBLINE($LIBBLAS_LIBDIR, RINDOW_OPENBLAS_SHARED_LIBADD)
  PHP_EVAL_INCLINE($LIBBLAS_CFLAGS)
  dnl # RINDOW_OPENBLAS_SHARED_LIBADD="-fopenmp $RINDOW_OPENBLAS_SHARED_LIBADD"
  AC_MSG_CHECKING([LIBBLAS_LIBDIR])
  AC_MSG_RESULT($LIBBLAS_LIBDIR)

  AC_DEFINE(CL_TARGET_OPENCL_VERSION, 120, [ Target OpenCL version 1.2 ])
  PHP_EVAL_LIBLINE($LIBOPENCL_LIBDIR" "$LIBCLBLAST_LIBDIR, RINDOW_CLBLAST_SHARED_LIBADD)
  PHP_EVAL_INCLINE($LIBOPENCL_CFLAGS" "$LIBCLBLAST_CFLAGS)

  dnl # PHP_ADD_INCLUDE($RINDOW_OPENBLAS_DIR/include)
  AC_MSG_CHECKING(for Interop/Polite/Math/Matrix.h)
  if test -f "PHP_EXT_SRCDIR(rindow_openblas)/vendor/interop-phpobjects/polite-math/include/Interop/Polite/Math/Matrix.h" ; then
    AC_MSG_RESULT(ok)
    PHP_ADD_INCLUDE(PHP_EXT_SRCDIR(rindow_openblas)[/vendor/interop-phpobjects/polite-math/include])
  else
    AC_MSG_RESULT(no)
    AC_MSG_ERROR(Interop/Polite/Math/Matrix.h not found. Please type "composer update")
  fi

  if test "$PHP_RINDOW_MATLIB" != "no"; then
    if test "$PHP_RINDOW_MATLIB" == "yes"; then
      AC_MSG_ERROR([You must specify a path when using --with-rindow_matlib])
    fi
    AC_MSG_CHECKING(for rindow/matlib.h)
    if test -f "$PHP_RINDOW_MATLIB/include/rindow/matlib.h" ; then
      AC_MSG_RESULT(ok)
      PHP_ADD_INCLUDE($PHP_RINDOW_MATLIB/include)
      PHP_ADD_INCLUDE($PHP_RINDOW_MATLIB/src)
      PHP_ADD_INCLUDE($PHP_RINDOW_MATLIB/build/src)
    else
      AC_MSG_RESULT(no)
      AC_MSG_ERROR(rindow/matlib.h not found. Please specify directory by --with-rindow_matlib option)
    fi
  fi

  PHP_SUBST(RINDOW_OPENBLAS_SHARED_LIBADD)
  

  dnl # In case of no dependencies
  AC_DEFINE(HAVE_RINDOW_OPENBLAS, 1, [ Have rindow_openblas support ])

  RINDOW_OPENBLAS_SOURCES="\
     rindow_openblas.c \
     src/Rindow/OpenBLAS/Buffer.c \
     src/Rindow/OpenBLAS/Blas.c \
     src/Rindow/OpenBLAS/Lapack.c \
     src/Rindow/OpenBLAS/Math.c \
  "

  dnl # PHP_NEW_EXTENSION(rindow_openblas, $RINDOW_OPENBLAS_SOURCES, $ext_shared,, -fopenmp -msse2)
  PHP_NEW_EXTENSION(rindow_openblas, $RINDOW_OPENBLAS_SOURCES, $ext_shared)
fi
