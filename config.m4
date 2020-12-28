dnl config.m4 for extension rindow_openblas

PHP_ARG_ENABLE(rindow_openblas, whether to enable rindow_openblas support,
dnl Make sure that the comment is aligned:
[  --enable-rindow_openblas          Enable rindow_openblas support], no)

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
  PHP_EVAL_LIBLINE($LIBBLAS_LIBDIR, RINDOW_OPENBLAS_SHARED_LIBADD)
  PHP_EVAL_INCLINE($LIBBLAS_CFLAGS)

  dnl # PHP_ADD_INCLUDE($RINDOW_OPENBLAS_DIR/include)
  AC_MSG_CHECKING(for Interop/Polite/Math/Matrix.h)
  if test -f "PHP_EXT_SRCDIR(rindow_openblas)/vendor/interop-phpobjects/polite-math/include/Interop/Polite/Math/Matrix.h" ; then
    AC_MSG_RESULT(ok)
    PHP_ADD_INCLUDE(PHP_EXT_SRCDIR(rindow_openblas)[/vendor/interop-phpobjects/polite-math/include])
  else
    AC_MSG_RESULT(no)
    AC_MSG_ERROR(Interop/Polite/Math/Matrix.h not found. Please type "composer update")
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

  PHP_NEW_EXTENSION(rindow_openblas, $RINDOW_OPENBLAS_SOURCES, $ext_shared)
fi
