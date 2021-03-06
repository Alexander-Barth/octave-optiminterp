#! /bin/sh

## Generate ./configure
rm -f configure.ac
echo "dnl --- DO NOT EDIT --- Automatically generated by autogen.sh" > configure.ac
cat configure.base >> configure.ac
cat <<EOF >> configure.ac
  AC_OUTPUT(\$CONFIGURE_OUTPUTS)
  dnl XXX FIXME XXX chmod is not in autoconf's list of portable functions

  echo " "
  echo "  \"\\\$prefix\" is \$prefix"
  echo "  \"\\\$exec_prefix\" is \$exec_prefix"
  AC_MSG_RESULT([\$STATUS_MSG

find . -name NOINSTALL -print    # shows which toolboxes won't be installed
])
EOF

autoconf

rm -f configure.ac
