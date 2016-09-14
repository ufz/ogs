def linux(buildDir, target = null, cmd = "make -j \$(nproc)") {
    if (target == null) {
        target = 'all'
        if (helper.isRelease())
            target = 'package'
    }
    sh "cd ${buildDir} && ${cmd} ${target}"
}

this.helper = load 'scripts/jenkins/lib/helper.groovy'
return this;
