#!/usr/bin/env python3
import sys, inspect
from pathlib import Path

local_dir = Path(inspect.getfile(inspect.currentframe()))
root_dir = local_dir.absolute().parents[3]
sys.path.append(str(root_dir))

import subprocess, traceback
import os


def run(command= 'ls', cmd_args= [ '-l', ] ):
    try:
        subprocess.run(
            [ command ] + cmd_args,
            stdin=sys.stdin,
            stdout=sys.stdout,
            stderr=sys.stderr,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        print(
            "\n^ ERROR ABOVE ^\nThe above error occured during execution of the subprocess, the wrapper stacktrace is:",
            file=sys.stderr,
        )
        traceback.print_tb(e.__traceback__, file=sys.stderr)
        print("Error:", e, file=sys.stderr)
    finally:
    	pass


if __name__ == "__main__":
    cmd = 'ls'
    cmd_args = ['-l', local_dir, ]
    run(command= cmd, cmd_args= cmd_args)
