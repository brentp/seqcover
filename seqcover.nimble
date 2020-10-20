import os

version       = "0.0.3"
author        = "Brent Pedersen"
description   = "find depth support for DUP/DEL/CNV calls that use PE/SR"
license       = "MIT"


# Dependencies

requires "https://github.com/brentp/d4-nim >= 0.0.1", "argparse", "binaryheap", "hts"
srcDir = "src"
installExt = @["nim"]

bin = @["seqcover"]


skipDirs = @["tests"]


#task test, "run the tests":
#  exec "nim c --lineDir:on --debuginfo -r --threads:on tests/all"
