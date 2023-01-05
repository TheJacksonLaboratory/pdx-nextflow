#! /bin/bash

# set_permissions.sh
# A component of cga, to set the permissions of any files / directories
# that need to be tweaked to create a properly configured cga installation.

# Must be invoked with a git version number.
if [ "X${1}" = "X" ]; then
    echo "USAGE: ${0} cga-version-number"
    exit 1
fi

BASE_DIR=/opt/compsci/cga/${1}/vcs_connector/setup

if [ ! -e ${BASE_DIR} ] ; then
    echo "ERROR: $1 does not appear to be a cga version."
    echo "    Can't find test file ${BASE_DIR}"
    exit 2
fi

# Set up the vcs_connector config file permissions.
chgrp -R exported-sequence-data ${BASE_DIR}
chmod 775 ${BASE_DIR}
chmod 664 ${BASE_DIR}/config.properties

echo Done.
