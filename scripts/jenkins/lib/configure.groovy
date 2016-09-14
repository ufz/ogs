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

return this;
