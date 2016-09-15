def linux(buildDir, cmakeOptions, generator = 'Unix Makefiles', conan_args = null, keepBuildDir = false) {
    if (keepBuildDir == false)
        sh("""rm -rf ${buildDir}
              mkdir ${buildDir}""".stripIndent())
    if (conan_args != null)
        sh("""cd ${buildDir}
              conan install ../ogs ${conan_args}""".stripIndent())
    sh """cd ${buildDir}
          cmake ../ogs -G "${generator}" ${defaultCMakeOptions} ${cmakeOptions}"""
}

def win(buildDir, cmakeOptions, generator, conan_args = null, keepBuildDir = false) {
    if (keepBuildDir == false)
        bat("""rd /S /Q ${buildDir}
               mkdir ${buildDir}""".stripIndent())

    if (conan_args != null)
        bat("""cd ${buildDir}
               conan install ../ogs ${conan_args}""".stripIndent())

    vcvarsallParam = "amd64"
        if (buildDir.endsWith("32"))
            vcvarsallParam = "x86"

    bat """set path=%path:\"=%
           call "%vs120comntools%..\\..\\VC\\vcvarsall.bat" ${vcvarsallParam}
           cd ${buildDir}
           cmake ../ogs -G "${generator}" ${defaultCMakeOptions} ${cmakeOptions}"""
}

return this;
