AC_DEFUN([AC_COVERAGE],
[
  AC_ARG_ENABLE(coverage,
                [  --enable-coverage       turn on -fprofile-arcs -ftest-coverage ],
                [case "${enableval}" in
                  yes) ENABLE_COVERAGE=1 ;;
                  no) ENABLE_COVERAGE=0 ;;
                  *) AC_MSG_ERROR([bad value ${enableval} for --enable-coverage]) ;;
                esac],
                [ENABLE_COVERAGE=2])

  AC_SUBST([ENABLE_COVERAGE])

  if test "x[$]ENABLE_COVERAGE" = "x1"
    then
      CFLAGS="`echo \"[$]CFLAGS\" | perl -pe 's/-O\d+//g;'` -fprofile-arcs -ftest-coverage"
    fi
])
