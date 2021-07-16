def getFileNameFromAbsPath(filepath):
    # Input: Absolute path of file (May not exist for now)
    # Return: file name without dir path
    return filepath.split("/")[-1]

def splitFileSuffix(InputfileName):
    # Input: only file name. (if multiple ".", use: for strN in InputfileName.split("."): ...)
    # Return: The file name that remove original suffix (ie, .fa, .fq, etc).
    return str(InputfileName.split(".")[0])

def splitFileSuffixThenAddNewSuffix(InputfileName, suffix):
    # Input: only file name.
    # Return: Remove original suffix (ie, .fa, .fq, etc) and add a new suffix
    oriFilePart = str(InputfileName.split(".")[0])
    return  oriFilePart + str(suffix)

def newPathName(path, suffix):
    import re
    returnPath = ''
    strList = path.split(".")
    if len(strList) < 2:
        # is mean no suffix:
        returnPath = strList[0] + suffix
    else:
        # cut off last . part:
        for strN in strList[:-1]:
            returnPath += (strN + ".")
        returnPath = re.sub("\.$","",returnPath)
        returnPath += suffix
    return returnPath

