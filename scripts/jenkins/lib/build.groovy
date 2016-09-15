def linux(buildDir, target = null, cmd = "make -j \$(nproc)") {
    if (target == null) {
        target = 'all'
        if (helper.isRelease())
            target = 'package'
    }
    sh "cd ${buildDir} && ${cmd} ${target}"
}

def win(buildDir, target = null) {
    targetString = ""
    if (target == null && helper.isRelease())
        targetString = "--target package"
    else
        targetString = "--target ${target}"
    bat("""set path=%path:\"=%
           call "%vs120comntools%..\\..\\VC\\vcvarsall.bat" x86_amd64
           cd ${buildDir}
           cmake --build . --config Release ${targetString}""".stripIndent())
}

this.helper = load 'scripts/jenkins/lib/helper.groovy'
return this;
