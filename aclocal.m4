dnl
dnl   This file is part of CONCORDE
dnl
dnl   (c) Copyright 1995--1999 by David Applegate, Robert Bixby,
dnl   Vasek Chvatal, and William Cook
dnl
dnl   Permission is granted for academic research use.  For other uses,
dnl   contact the authors for licensing options.
dnl
dnl   Use at your own risk.  We make no guarantees about the
dnl   correctness or usefulness of this code.
dnl

AC_DEFUN(concorde_CCDEFAULTS,
[AC_REQUIRE([AC_CANONICAL_HOST])
AC_PROVIDE([concorde_CCDEFAULTS])
dnl some canonical hosts
dnl binky     mips-sgi-irix5.3
dnl trick     powerpc-ibm-aix4.2.1.0
dnl badenia   rs6000-ibm-aix4.2.1.0
dnl r53aix01  rs6000-ibm-aix3.2.5
dnl macdoon   i386-pc-solaris2.6
dnl cycle     i686-pc-linux-gnu
dnl sp01      i386-unknown-freebsd2.2.6
dnl truck     sparc-sun-solaris2.6
dnl north     mips-sgi-irix6.2
dnl fry       mips-sgi-irix6.2
dnl pickle    alphaev56-dec-osf4.0a
dnl malarkey  alphaev6-dec-osf4.0f
dnl clayton   alpha-unknown-linux-gnu

AC_MSG_CHECKING([for prespecified compiler options])
found_ccdefault=no
case "$host" in
  alphaev6-*-osf4*)   : ${CC=cc.alt} ;;
  *-*-osf4*)          : ${CC=cc}  ;;
  alpha*-*-linux* )   : ${CC=gcc} ;;
  i*-*-freebsd*)      : ${CC=gcc} ;;
  i*-*-linux*)        : ${CC=gcc} ;;
  i*-*-solaris*)      : ${CC=gcc} ;;
  i*-pc-cygwin32*)    : ${CC=gcc} ;;
  *-*-aix*)           : ${CC=cc}  ;;
  *-*-irix5*)         : ${CC=cc}  ;;
  *-*-irix6*)         : ${CC=cc}  ;;
  sparc-*-solaris2.6) : ${CC=cc}  ;;
  sparc-*-solaris*)
      if test "x$with_purify" = xyes ; then
          : ${CC=cc}
      else
          : ${CC=gcc}
      fi
      ;;
  sparc-*-sunos*)   : ${CC=acc} ;;
esac

if test "x$enable_debugger" = xyes ; then
   : ${CFLAGS="-g"}
fi

if test "x$CC" = xcc ; then
   case "$host" in
     *-*-osf4*)        : ${CFLAGS="-arch host -O4 -g3"}
                       : ${CPPFLAGS="-std1 -warnprotos -portable"}
                       found_ccdefault=yes ;;
     *-*-solaris*)     : ${CFLAGS="-xO2"}
                       : ${CPPFLAGS="-v -Xc"}
                       found_ccdefault=yes ;;
     *-*-aix*)         : ${CFLAGS="-O"}
                       : ${CPPFLAGS="-qlanglvl=ansi"}
                       : ${LDFLAGS="-bmaxdata:0x80000000"} 
                       found_ccdefault=yes ;;
     *-*-irix5*)       : ${CFLAGS="-O2"}
                       : ${CPPFLAGS="-ansi -fullwarn"} 
                       found_ccdefault=yes ;;
     *-*-irix6*)       : ${CFLAGS="-O2"}
                       : ${CPPFLAGS="-64 -mips4 -xansi -fullwarn -woff 1174,1209"}
                       : ${LDFLAGS="-64 -mips4 -Wl,-woff,84"}
                       found_ccdefault=yes ;;
     sparc-*-sunos*)   : ${CFLAGS="-O2"}
                       : ${CPPFLAGS=""}
                       found_ccdefault=yes ;;
   esac
fi
if test "x$CC" = xgcc ; then
   case "$host" in
     dnl -pedantic -Wtraditional -Wmissing-prototypes -Wmissing-declarations -Wshadow
     dnl cause errors in header files on alpha-*-linux*
     alpha*-*-linux* ) : ${CFLAGS="-O3 -g"}
                       : ${CPPFLAGS="-ansi -Wall -W -Wstrict-prototypes -Wpointer-arith -Wnested-externs"}
                       found_ccdefault=yes ;;
     dnl -ansi disrupts gethostname, sigaction
     dnl -pedantic generates header errors in socket.h
     i*-*-linux* )     : ${CFLAGS="-O3 -g"}
                       : ${CPPFLAGS="-Wall -Wshadow -W -Wtraditional -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wpointer-arith -Wnested-externs -Wundef -Wcast-qual -Wcast-align -Wwrite-strings"}
                       found_ccdefault=yes ;;
     dnl -ansi causes errors in header files in *-*-irix6*
     *-*-irix6*)       : ${CFLAGS="-O3 -g"}
                       : ${CPPFLAGS="-pedantic -Wall -Wshadow -W -Wtraditional -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wpointer-arith -Wnested-externs -Wundef -Wcast-qual -Wcast-align -Wwrite-strings"}
                       : ${LDFLAGS="-Wl,-woff,84"}
                       found_ccdefault=yes ;;
     *)                : ${CFLAGS="-O3 -g"}
                       : ${CPPFLAGS="-ansi -pedantic -Wall -W -Wtraditional -Wundef -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wnested-externs"}
                       found_ccdefault=yes ;;
     esac
fi
if test "x$CC" = xacc ; then
   : ${CFLAGS="-O2"}
   : ${CPPFLAGS=""}
   found_ccdefault=yes
fi
if test "x$CC" = xcc.alt ; then
   case "$host" in
     *-*-osf4*)        : ${CFLAGS="-arch host -O4 -g3"}
                       : ${CPPFLAGS="-std1 -warnprotos -portable"}
                       found_ccdefault=yes ;;
   esac
fi

if test "x$found_ccdefault" = "xyes" ; then
   AC_MSG_RESULT(yes)
   echo "    Using default CC:       $CC"
   echo "    Using default CPPFLAGS: $CPPFLAGS"
   echo "    Using default CFLAGS:   $CFLAGS"
   echo "    Using default LDFLAGS:  $LDFLAGS"
else
   AC_MSG_RESULT(no)
fi
])

AC_DEFUN(concorde_C_VOID,
[AC_REQUIRE([AC_PROG_CC])
AC_PROVIDE([concorde_C_VOID])
AC_CACHE_CHECK([for void], cc_cv_c_void,
[AC_TRY_COMPILE(,[void x();],cc_cv_c_void=yes, cc_cv_c_void=no)])
if test $cc_cv_c_void = no; then
  AC_DEFINE(void,int)
fi
])

dnl concorde_CHECK_LIB(LIBRARY, FUNCTION [, ACTION-IF-FOUND
dnl         [, ACTION-IF-NOT-FOUND [, OTHER-LIBRARIES]]])
dnl derived from AC_CHECK_LIB, but only uses the library if it is also
dnl needed (ie, you can't get FUNCTION without LIBRARY)
AC_DEFUN(concorde_CHECK_LIB,
[AC_REQUIRE([AC_PROG_CC])
AC_PROVIDE([concorde_CHECK_LIB])
AC_MSG_CHECKING([for -l$1])
dnl Use a cache variable name containing both the library and function name,
dnl because the test really is for library $1 defining function $2, not
dnl just for library $1.  Separate tests with the same $1 and different $2s
dnl may have different results.
cc_lib_var=`echo $1['_']$2 | tr './+\055' '__p_'`
AC_CACHE_VAL(cc_cv_lib_$cc_lib_var,
[AC_TRY_LINK(dnl
ifelse([$2], [main], , dnl Avoid conflicting decl of main.
[/* Override any gcc2 internal prototype to avoid an error.  */
]ifelse(AC_LANG, CPLUSPLUS, [#ifdef __cplusplus
extern "C"
#endif
])dnl
[/* We use char because int might match the return type of a gcc2
    builtin and then its argument prototype would still apply.  */
char $2();
]),
	    [$2()],
	    eval "cc_cv_lib_$cc_lib_var=no",
[cc_save_LIBS="$LIBS"
LIBS="-l$1 $5 $LIBS"
AC_TRY_LINK(dnl
ifelse([$2], [main], , dnl Avoid conflicting decl of main.
[/* Override any gcc2 internal prototype to avoid an error.  */
]ifelse(AC_LANG, CPLUSPLUS, [#ifdef __cplusplus
extern "C"
#endif
])dnl
[/* We use char because int might match the return type of a gcc2
    builtin and then its argument prototype would still apply.  */
char $2();
]),
	    [$2()],
	    eval "cc_cv_lib_$cc_lib_var=yes",
	    eval "cc_cv_lib_$cc_lib_var=no")
LIBS="$cc_save_LIBS"
])dnl
])dnl
if eval "test \"`echo '$cc_cv_lib_'$cc_lib_var`\" = yes"; then
  AC_MSG_RESULT(yes)
  ifelse([$3], ,
[changequote(, )dnl
  cc_tr_lib=HAVE_LIB`echo $1 | tr 'abcdefghijklmnopqrstuvwxyz' 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'`
changequote([, ])dnl
  AC_DEFINE_UNQUOTED($cc_tr_lib)
  LIBS="-l$1 $LIBS"
], [$3])
else
  AC_MSG_RESULT(no)
ifelse([$4], , , [$4
])dnl
fi
])

dnl concorde_TRY_CPLEX tries to link to cplex in the specified library
AC_DEFUN(concorde_TRY_CPLEX,
[AC_REQUIRE([AC_PROG_CC])
cc_save_LIBS="$LIBS"
LIBS="$1 $LIBS"
AC_TRY_LINK(dnl
ifelse(AC_LANG, CPLUSPLUS, [#ifdef __cplusplus
extern "C"
#endif
])dnl
[
char CPXloadprob();
], [CPXloadprob()], cc_try_cplex=yes, cc_try_cplex=no)
LIBS=$cc_save_LIBS
])

dnl concorde_CHECK_CPLEX tests whether cplex is available
AC_DEFUN(concorde_CHECK_CPLEX,
[AC_REQUIRE([AC_PROG_CC])
AC_PROVIDE([concorde_CHECK_CPLEX])
cc_cplexloc=$1
cc_loc_var=`echo $cc_cplexloc | tr './+\055' '__p_'`
AC_MSG_CHECKING([cplex])
AC_CACHE_VAL(cc_cv_cplex_$cc_loc_var,
[if test "$cc_cplexloc" = "yes" ; then
   cc_cplex_lib="-lcplex"
else
   cc_cplex_lib="$cc_cplexloc/${CPLEX_LIBNAME-libcplex.a}"
fi
concorde_TRY_CPLEX($cc_cplex_lib)
if test "$cc_try_cplex" = "yes" ; then
   eval "cc_cv_cplex_$cc_loc_var=yes"
else
   concorde_TRY_CPLEX($cc_cplex_lib -lpthread)
   if test "$cc_try_cplex" = "yes" ; then
      eval "cc_cv_cplex_$cc_loc_var=pthread"
   else
      concorde_TRY_CPLEX($cc_cplex_lib -lpthreads)
      if test "$cc_try_cplex" = "yes" ; then
         eval "cc_cv_cplex_$cc_loc_var=pthreads"
      else
         eval "cc_cv_cplex_$cc_loc_var=no"
      fi
   fi
fi
])dnl
if eval "test \"`echo '$cc_cv_cplex_'$cc_loc_var`\" = no"; then
   AC_MSG_RESULT(no)
   LPSOLVER_LIB=""
   LPSOLVER_INCFLAG=""
   cc_cplex_found="no"
else
   if test "$cc_cplexloc" = "yes" ; then
      LPSOLVER_LIB="-lcplex"
      LPSOLVER_INCFLAG=""
   else
      LPSOLVER_LIB="$cc_cplexloc/${CPLEX_LIBNAME-libcplex.a}"
      LPSOLVER_INCFLAG="-I$cc_cplexloc"
   fi
   cc_cplex_found="yes"

   if eval "test \"`echo '$cc_cv_cplex_'$cc_loc_var`\" = pthread"; then
      AC_MSG_RESULT(yes, with -lpthread)
      LPSOLVER_LIB="$LPSOLVER_LIB -lpthread"
   else
      if eval "test \"`echo '$cc_cv_cplex_'$cc_loc_var`\" = pthreads"; then
         AC_MSG_RESULT(yes, with -lpthreads)
         LPSOLVER_LIB="$LPSOLVER_LIB -lpthreads"
      else
         AC_MSG_RESULT(yes)
      fi
   fi
fi
])

dnl concorde_CPLEX_VERSION tests which version of cplex we are using.
AC_DEFUN(concorde_CPLEX_VERSION,
[AC_REQUIRE([AC_PROG_CC])
AC_PROVIDE([concorde_CPLEX_VERSION])
AC_MSG_CHECKING([cplex version])
cc_cplexloc=$1
cc_loc_var=`echo $cc_cplexloc['_']version | tr './+\055' '__p_'`
AC_CACHE_VAL(cc_cv_cplex_$cc_loc_var,
[cc_save_LIBS="$LIBS"
LIBS="$LPSOLVER_LIB $LIBS"
AC_TRY_LINK(dnl
ifelse(AC_LANG, CPLUSPLUS, [#ifdef __cplusplus
extern "C"
#endif
])dnl
[
char CPXcopybase();
], [CPXcopybase()], eval "cc_cv_cplex_$cc_loc_var=6", [AC_TRY_LINK(dnl
ifelse(AC_LANG, CPLUSPLUS, [#ifdef __cplusplus
extern "C"
#endif
])dnl
[
char CPXcreateprob();
], [CPXcreateprob()], eval "cc_cv_cplex_$cc_loc_var=5", eval "cc_cv_cplex_$cc_loc_var=4")])
LIBS=$cc_save_LIBS
])
if eval "test \"`echo '$cc_cv_cplex_'$cc_loc_var`\" = 6"; then
   AC_MSG_RESULT([>= 6])
   cc_cplex_version=6
else
   if eval "test \"`echo '$cc_cv_cplex_'$cc_loc_var`\" = 5"; then
      AC_MSG_RESULT([5])
      cc_cplex_version=5
   else
      AC_MSG_RESULT([<= 4])
      cc_cplex_version=4
   fi
fi
])

dnl concorde_SUFFIXES is not currently complete - it just selects
dnl default values
AC_DEFUN(concorde_SUFFIXES,
[AC_REQUIRE([AC_PROG_CC])
AC_PROVIDE([concorde_SUFFIXES])
AC_MSG_CHECKING([for file suffixes])
OBJ_SUFFIX=o
LIB_SUFFIX=a
EXE_SUFFIX=
AC_SUBST(OBJ_SUFFIX)
AC_SUBST(LIB_SUFFIX)
AC_SUBST(EXE_SUFFIX)
AC_MSG_RESULT($OBJ_SUFFIX[, ]$LIB_SUFFIX[, ]$EXE_SUFFIX)
])

AC_DEFUN(concorde_MISSING_PROTOTYPE,
[AC_REQUIRE([AC_PROG_CC])
AC_PROVIDE([concorde_MISSING_PROTOTYPE])
AC_CACHE_CHECK([for missing ]$1[ prototype],cc_cv_proto_$1,
[AC_TRY_COMPILE([$2],[$3],
eval "cc_cv_proto_$1=yes",
eval "cc_cv_proto_$1=no")])
if eval "test \"`echo '$cc_cv_proto_'$1`\" = yes" ; then
   upname=`echo $1 | tr 'abcdefghijklmnopqrstuvwxyz' 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'`
   AC_DEFINE_UNQUOTED(CC_PROTO_$upname)
fi
])

AC_DEFUN(concorde_HEADER_PTHREAD,
[AC_REQUIRE([AC_PROG_CC])
AC_PROVIDE([concorde_HEADER_PTHREAD])
AC_CACHE_CHECK([whether signal.h needs to be included before pthread.h],
cc_cv_header_pthread,
[AC_TRY_COMPILE([#include <pthread.h>],,cc_cv_header_pthread=no,cc_cv_header_pthread=yes)
if test $cc_cv_header_pthread = yes; then
   AC_TRY_COMPILE([#include <signal.h>
#include <pthread.h>],,cc_cv_header_pthread=yes,cc_cv_header_pthread=fail)
fi
])
if test $cc_cv_header_pthread = fail; then
   AC_MSG_ERROR([could not successfully include pthread.h])
fi
if test $cc_cv_header_pthread = yes; then
   AC_DEFINE(CC_SIGNAL_BEFORE_PTHREAD)
fi
])

AC_DEFUN(concorde_SIGNAL_METHOD,
[AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([AC_TYPE_SIGNAL])
AC_PROVIDE([concorde_SIGNAL_METHOD])
AC_CACHE_CHECK([usable sigaction()],cc_cv_signal_sigaction,[AC_TRY_COMPILE(
[#ifdef HAVE_SIGNAL_H
#include <signal.h>
#endif
],
[struct sigaction n;
 extern RETSIGTYPE (h)(int signum);
 int x;
 x = sigemptyset(&n.sa_mask);
 n.sa_handler=h;
 n.sa_flags=4;
 x = sigaction(0,&n,(struct sigaction*)0);
],
[cc_cv_signal_sigaction=yes],
[cc_cv_signal_sigaction=no])])
if test "x$cc_cv_signal_sigaction" = "xyes" ; then
   AC_DEFINE(CCSIGNAL_SIGACTION)
else
AC_CACHE_CHECK([usable signal()],cc_cv_signal_signal,[AC_TRY_COMPILE(
[#ifdef HAVE_SIGNAL_H
#include <signal.h>
#endif
],
[
 extern RETSIGTYPE (*h)(int signum);
 int x;
 h=signal(0,h);
 if (h == SIG_ERR) x=1;
],
[cc_cv_signal_signal=yes],
[cc_cv_signal_signal=no])])
if test "x$cc_cv_signal_signal" = "xyes" ; then
   AC_DEFINE(CCSIGNAL_SIGNAL)
else
   AC_DEFINE(CCSIGNAL_NONE)
fi
fi
])

AC_DEFUN(concorde_ATTRIBUTE,
[AC_REQUIRE([AC_PROG_CC])
AC_PROVIDE([concorde_ATTRIBUTE])
AC_CACHE_CHECK([for __attribute__],cc_cv_attributes,
[AC_TRY_COMPILE(,[int __attribute__ ((unused)) x (__attribute__ ((unused)) int y){return 0;}],cc_cv_attributes=yes,cc_cv_attributes=no)])
if test $cc_cv_attributes = yes; then
   AC_DEFINE(CC_ATTRIBUTE)
fi
])
