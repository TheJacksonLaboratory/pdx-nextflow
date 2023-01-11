#!/bin/bash

export PERL5LIB=$PERL5LIB:/opt/compsci/SOAPfuse/1.27/source/bin/perl_module

TMP_STDERR=${TMPDIR}/soapfuse_stderr.txt

# Capture stderr in case we need to grep it for a string.
perl /opt/compsci/SOAPfuse/1.27/SOAPfuse-RUN.pl -tp $TMPDIR $@ 2> ${TMP_STDERR}
STS=$?

# Now put the captured error text to stderr.
cat ${TMP_STDERR} >&2

# If we completed successfully, we're done.
[[ "X${STS}" == "X0" ]] && rm ${TMP_STDERR} && exit 0

# The perl job didn't complete cleanly, but one failure that we ignore is
# when no fusions are found.

grep -q "Reason: FD2_fusion_result_is_empty" ${TMP_STDERR}
GSTS=$?

# Don't need the captured error text anymore.
rm ${TMP_STDERR}
[[ "X${GSTS}" == "X0" ]] && exit 0

# We failed with other than the expected failure.  Report it.
exit ${STS}
