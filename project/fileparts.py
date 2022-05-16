"""
## Get the various parts of a file with path string.
 
## @param filename.with.path A string of a filename with a path
 
## @return A list with the following components:
## \item{pathname}{The path name}
## \item{filename}{The file name}
## \item{fileext}{The file extension}
## @note This is a modified version of the same function in the matlab R-package.
 
## @family MATLAB
## @export
## @keywords internal


## Function written to match MATLAB function fileparts()
## Author: Caiya Zhang, Yuchen Zheng
"""

import os
import re

def fileparts(filename_with_path):
    pathname = os.path.dirname(filename_with_path)
    filename = os.path.basename(filename_with_path)
    fileext = re.sub(".*(\\.[^\\.]*)$", "\\1", filename)
    filename = re.sub("(.*)(\\.[^\\.]*)$", "\\1", filename)
    if fileext == filename:
        fileext = ""
    fileparts_list = {"pathname": pathname, "filename": filename, "fileext": fileext}
    return fileparts_list
