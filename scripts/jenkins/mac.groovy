defaultCMakeOptions = '-DCMAKE_BUILD_TYPE=Release -DOGS_CPU_ARCHITECTURE=core2 -DOGS_LIB_BOOST=System -DOGS_LIB_VTK=System -DOGS_DOWNLOAD_ADDITIONAL_CONTENT=ON'

node('mac && conan') {
    stage 'Checkout (Mac)'
    dir('ogs') { checkout scm }

    stage 'Configure (Mac)'
    configure 'build', '', 'Ninja', ''

    stage 'CLI (Mac)'
    build 'build'

    stage 'Test (Mac)'
    build 'build', 'tests ctest'

    stage 'Data Explorer (Mac)'
    configure 'build', '-DOGS_BUILD_GUI=ON -DOGS_BUILD_UTILS=ON -DOGS_BUILD_TESTS=OFF',
        'Ninja', '', true
    build 'build'

    if (env.BRANCH_NAME == 'master' || env.BRANCH_NAME.contains('release')) {
        stage 'Release (Mac)'
        archive 'build/*.zip,build/*.dmg'
    }

    stage 'Post (Mac)'
    publishTestReports 'build/Testing/**/*.xml', 'build/Tests/testrunner.xml',
        'ogs/scripts/jenkins/msvc-log-parser.rules'
}

def configure(buildDir, cmakeOptions, generator, conan_args = null, keepBuildDir = false) {
    if (keepBuildDir == false)
        sh("""rm -rf ${buildDir}
              mkdir ${buildDir}""".stripIndent())
    if (conan_args != null)
        sh("""cd ${buildDir}
              conan install ../ogs ${conan_args}""".stripIndent())
    sh """cd ${buildDir}
          cmake ../ogs -G "${generator}" ${defaultCMakeOptions} ${cmakeOptions}"""
}

def build(buildDir, target = null) {
    if (target == null) {
        target = 'all'
        if (env.BRANCH_NAME == 'master' || env.BRANCH_NAME.contains('release'))
            target = 'package'
    }
    sh "cd ${buildDir} && ninja ${target}"
}

// *** Helper functions ***
def publishTestReports(ctestPattern, gtestPattern, parseRulefile) {
    step([$class: 'XUnitPublisher', testTimeMargin: '3000', thresholdMode: 1,
        thresholds: [
            [$class: 'FailedThreshold', failureNewThreshold: '', failureThreshold: '',
                unstableNewThreshold: '', unstableThreshold: ''],
            [$class: 'SkippedThreshold', failureNewThreshold: '', failureThreshold: '',
                unstableNewThreshold: '', unstableThreshold: '']],
        tools: [
            [$class: 'CTestType', deleteOutputFiles: true, failIfNotNew: true, pattern:
                "${ctestPattern}", skipNoTestFiles: false, stopProcessingIfError: true],
            [$class: 'GoogleTestType', deleteOutputFiles: true, failIfNotNew: true, pattern:
                "${gtestPattern}", skipNoTestFiles: false, stopProcessingIfError: true]]
    ])

    step([$class: 'LogParserPublisher', failBuildOnError: true, unstableOnWarning: false,
        projectRulePath: "${parseRulefile}", useProjectRule: true])
}
// *** End helper functions ***
