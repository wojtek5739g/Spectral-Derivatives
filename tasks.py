import cffi
import invoke
import pathlib
import sys
import os
import shutil
import re
import glob

ffi = cffi.FFI()

this_dir = pathlib.Path().absolute()
h_file_name = this_dir / "cmult.h"
with open(h_file_name) as h_file:
    ffi.cdef(h_file.read())

ffi.set_source(
    "cffi_example",
    '#include "cmult.h',
    libraries=["cmult"],
    library_dirs=[this_dir.as_posix()],
    extra_link_args=["-Wl,-rpath,."],
)

ffi.compile()

