"""
Purpose :      Run external commands.

Copyright (C): 2016 - Gert Hulselmans
"""

from __future__ import print_function

import os
import sys

if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess


def run_cmd(cmd, stdin=None):
    """
    Execute the passed command line and exit when the command failed.

    :param cmd: list of command line parameters.
    :param stdin: Standard input string or None.
    :return: stdout_data, stderr_data (on error print stdout_data and stderr_data to standard error)
    """

    try:
        pid = subprocess.Popen(args=cmd,
                               bufsize=1,
                               executable=None,
                               stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               preexec_fn=None,
                               close_fds=False,
                               shell=False,
                               cwd=None,
                               env=None,
                               universal_newlines=False,
                               startupinfo=None,
                               creationflags=0)
        stdout_data, stderr_data = pid.communicate(stdin)
    except OSError as msg:
        print("\nERROR during execution of '" + ' '.join(cmd) + "': " + str(msg), file=sys.stderr)
        sys.exit(1)

    if pid.returncode != 0:
        print("\nERROR during execution of '" + ' '.join(cmd) + "':", file=sys.stderr)
        print("\nStandard output:\n----------------", file=sys.stderr)
        print(stdout_data, file=sys.stderr)
        print("\nStandard error:\n---------------", file=sys.stderr)
        print(stderr_data, file=sys.stderr)

        print("\nERROR during execution of '" + ' '.join(cmd) + "':", file=sys.stderr)
        print("\nStandard output:\n----------------", file=sys.stderr)
        print(stdout_data, file=sys.stderr)
        print("\nStandard error:\n---------------", file=sys.stderr)
        print(stderr_data, file=sys.stderr)
        sys.exit(1)

    return stdout_data, stderr_data
