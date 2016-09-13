defaultCMakeOptions = '-DCMAKE_BUILD_TYPE=Release -DOGS_32_BIT=ON -DOGS_LIB_BOOST=System ' +
    '-DOGS_LIB_VTK=System -DOGS_DOWNLOAD_ADDITIONAL_CONTENT=ON'

node('win && conan') {
    stage 'Checkout (Win)'
    dir('ogs') { checkout scm }

    if (env.BRANCH_NAME == 'master' || env.BRANCH_NAME.contains('release')) {
        stage 'Data Explorer 32-bit (Win)'
        withEnv(getEnv()) {
            configure 'build-32', '-DOGS_BUILD_GUI=ON -DOGS_BUILD_UTILS=ON ' +
                '-DOGS_BUILD_TESTS=OFF', 'Ninja', '-u -s build_type=Release -s compiler="Visual ' +
                'Studio" -s compiler.version=12 -s arch=x86'
            build 'build-32', 'package'
        }
        archive 'build-32/*.zip'
    }
}

def getEnv() {
    if (env.NODE_NAME == 'visserv3')
        qtdir = 'C:\\libs\\qt\\4.8\\msvc2013-x32'
    if (env.NODE_NAME == 'win1')
        qtdir = 'C:\\libs\\qt-4.8.7-x86-msvc2013\\qt-4.8.7-x86-msvc2013'

    return [
        "QTDIR=${qtdir}",
        'Path=$Path;$QTDIR\\bin',
        'CONAN_CMAKE_GENERATOR=Ninja'
    ]
}

def configure(buildDir, cmakeOptions, generator, conan_args = null) {
    bat("""rd /S /Q ${buildDir}
           mkdir ${buildDir}""".stripIndent())
    if (conan_args != null)
        bat("""cd ${buildDir}
               conan install ../ogs ${conan_args}""".stripIndent())
    if (generator == 'Ninja')
        bat """set path=%path:\"=%
           call "%vs120comntools%..\\..\\VC\\vcvarsall.bat" x86
           cd ${buildDir}
           cmake ../ogs -G "${generator}" ${defaultCMakeOptions} ${cmakeOptions}"""
}

def build(buildDir, target = null) {
    targetString = ""
    if (target != null)
        targetString = "--target ${target}"
    bat("""set path=%path:\"=%
           call "%vs120comntools%..\\..\\VC\\vcvarsall.bat" x86
           cd ${buildDir}
           cmake --build . --config Release ${targetString}""".stripIndent())
}
// *** End helper functions ***
