def isRelease () {
    if (env.BRANCH_NAME == 'master' || env.BRANCH_NAME.contains('release'))
        return true;
    return false;
}

def getEnv(arch = 'x64') {
    if (env.NODE_NAME == 'visserv3')
        qtdir = "C:\\libs\\qt\\4.8\\msvc2013-${arch}"
    if (env.NODE_NAME == 'win1') {
        if (arch == 'x32')
            qtdir = "C:\\libs\\qt-4.8.7-x86-msvc2013\\qt-4.8.7-x86-msvc2013"
        else
            qtdir = "C:\\libs\\qt-4.8.7-${arch}-msvc2013\\qt-4.8.7-${arch}-msvc2013"
    }

    return [
        "QTDIR=${qtdir}",
        'Path=$Path;$QTDIR\\bin',
        'CONAN_CMAKE_GENERATOR=Ninja',
        "ARCH=${arch}"
    ]
}

return this;
