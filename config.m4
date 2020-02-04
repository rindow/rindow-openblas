dnl config.m4 for extension rindow_openblas

dnl Comments in this file start with the string 'dnl'.
dnl Remove where necessary.

dnl If your extension references something external, use with:

dnl PHP_ARG_WITH(rindow_openblas, for rindow_openblas support,
dnl Make sure that the comment is aligned:
dnl [  --with-rindow_openblas             Include rindow_openblas support])

dnl Otherwise use enable:

PHP_ARG_ENABLE(rindow_openblas, whether to enable rindow_openblas support,
dnl Make sure that the comment is aligned:
[  --enable-rindow_openblas          Enable rindow_openblas support], no)

if test "$PHP_RINDOW_OPENBLAS" != "no"; then
  dnl Write more examples of tests here...

  dnl # get library FOO build options from pkg-config output
  dnl AC_PATH_PROG(PKG_CONFIG, pkg-config, no)
  dnl AC_MSG_CHECKING(for libfoo)
  dnl if test -x "$PKG_CONFIG" && $PKG_CONFIG --exists foo; then
  dnl   if $PKG_CONFIG foo --atleast-version 1.2.3; then
  dnl     LIBFOO_CFLAGS=\`$PKG_CONFIG foo --cflags\`
  dnl     LIBFOO_LIBDIR=\`$PKG_CONFIG foo --libs\`
  dnl     LIBFOO_VERSON=\`$PKG_CONFIG foo --modversion\`
  dnl     AC_MSG_RESULT(from pkgconfig: version $LIBFOO_VERSON)
  dnl   else
  dnl     AC_MSG_ERROR(system libfoo is too old: version 1.2.3 required)
  dnl   fi
  dnl else
  dnl   AC_MSG_ERROR(pkg-config not found)
  dnl fi
  dnl PHP_EVAL_LIBLINE($LIBFOO_LIBDIR, RINDOW_OPENBLAS_SHARED_LIBADD)
  dnl PHP_EVAL_INCLINE($LIBFOO_CFLAGS)

  dnl # --with-rindow_openblas -> check with-path
  dnl SEARCH_PATH="/usr/local /usr"     # you might want to change this
  dnl SEARCH_FOR="/include/rindow_openblas.h"  # you most likely want to change this
  dnl if test -r $PHP_RINDOW_OPENBLAS/$SEARCH_FOR; then # path given as parameter
  dnl   RINDOW_OPENBLAS_DIR=$PHP_RINDOW_OPENBLAS
  dnl else # search default path list
  dnl   AC_MSG_CHECKING([for rindow_openblas files in default path])
  dnl   for i in $SEARCH_PATH ; do
  dnl     if test -r $i/$SEARCH_FOR; then
  dnl       RINDOW_OPENBLAS_DIR=$i
  dnl       AC_MSG_RESULT(found in $i)
  dnl     fi
  dnl   done
  dnl fi
  dnl
  dnl if test -z "$RINDOW_OPENBLAS_DIR"; then
  dnl   AC_MSG_RESULT([not found])
  dnl   AC_MSG_ERROR([Please reinstall the rindow_openblas distribution])
  dnl fi

  dnl # --with-rindow_openblas -> add include path
  dnl PHP_ADD_INCLUDE($RINDOW_OPENBLAS_DIR/include)

  dnl # --with-rindow_openblas -> check for lib and symbol presence
  dnl LIBNAME=RINDOW_OPENBLAS # you may want to change this
  dnl LIBSYMBOL=RINDOW_OPENBLAS # you most likely want to change this

  dnl PHP_CHECK_LIBRARY($LIBNAME,$LIBSYMBOL,
  dnl [
  dnl   PHP_ADD_LIBRARY_WITH_PATH($LIBNAME, $RINDOW_OPENBLAS_DIR/$PHP_LIBDIR, RINDOW_OPENBLAS_SHARED_LIBADD)
  dnl   AC_DEFINE(HAVE_RINDOW_OPENBLASLIB,1,[ ])
  dnl ],[
  dnl   AC_MSG_ERROR([wrong rindow_openblas lib version or lib not found])
  dnl ],[
  dnl   -L$RINDOW_OPENBLAS_DIR/$PHP_LIBDIR -lm
  dnl ])
  dnl
  dnl PHP_SUBST(RINDOW_OPENBLAS_SHARED_LIBADD)

  dnl # In case of no dependencies
  AC_DEFINE(HAVE_RINDOW_OPENBLAS, 1, [ Have rindow_openblas support ])

  PHP_NEW_EXTENSION(rindow_openblas, rindow_openblas.c, $ext_shared)
fi
