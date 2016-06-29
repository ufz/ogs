defaultDockerArgs = '-v /home/jenkins/.ccache:/usr/src/.ccache'
defaultCMakeOptions = '-DOGS_LIB_BOOST=System -DOGS_LIB_VTK=System'

node('visserv3')
{
    step([$class: 'GitHubSetCommitStatusBuilder', statusMessage: [content: 'Started Jenkins MSVC build']])

    stage 'Checkout'
    dir('ogs') { checkout scm }

    withEnv(['ARCH=msvc2013-x64', 'CMAKE_LIBRARY_SEARCH_PATH=C:\\libs\\$ARCH', 'QTDIR=C:\\libs\\qt\\4.8\\$ARCH', 'Path=$Path;$QTDIR\\bin;$CMAKE_LIBRARY_SEARCH_PATH\\bin']) {

        stage 'Data Explorer'
        configure 'build-de', '-DOGS_BUILD_GUI=ON -DOGS_BUILD_CLI=OFF -DOGS_DOWNLOAD_ADDITIONAL_CONTENT=ON', 'Visual Studio 12 Win64'
        build 'build-de'

        stage 'CLI'
        configure 'build', '-DOGS_DOWNLOAD_ADDITIONAL_CONTENT=ON', 'Visual Studio 12 Win64'
        build 'build', 'tests'
        build 'build', 'ctest'

        if (env.BRANCH_NAME == 'master') {
            stage 'Package'
            build 'build', 'package'
            build 'build-de', 'package'
        }
    }

    stage 'Post'
    publishTestReports 'build/Testing/**/*.xml', 'build/Tests/testrunner.xml',
        'ogs/scripts/jenkins/msvc-log-parser.rules'

    archive 'build*/*.zip'
}

// *** Helper functions ***
def configure(buildDir, cmakeOptions, generator) {
    bat("""rd /S /Q ${buildDir}
           mkdir ${buildDir}
           cd ${buildDir}
           cmake ../ogs -G "${generator}" ${defaultCMakeOptions} ${cmakeOptions}""".stripIndent())
}

def build(buildDir, target=null) {
    targetString = ""
    if (target != null)
        targetString = "--target ${target}"
    bat("""cd ${buildDir}
           cmake --build . --config Release ${targetString}""".stripIndent())
}

def publishTestReports(ctestPattern, gtestPattern, parseRulefile) {
    step([$class: 'XUnitPublisher', testTimeMargin: '3000', thresholdMode: 1,
        thresholds: [
            [$class: 'FailedThreshold', failureNewThreshold: '', failureThreshold: '', unstableNewThreshold: '', unstableThreshold: ''],
            [$class: 'SkippedThreshold', failureNewThreshold: '', failureThreshold: '', unstableNewThreshold: '', unstableThreshold: '']],
        tools: [
            [$class: 'CTestType', deleteOutputFiles: true, failIfNotNew: true, pattern: "${ctestPattern}", skipNoTestFiles: false, stopProcessingIfError: true],
            [$class: 'GoogleTestType', deleteOutputFiles: true, failIfNotNew: true, pattern: "${gtestPattern}", skipNoTestFiles: false, stopProcessingIfError: true]]
    ])

    step([$class: 'LogParserPublisher', failBuildOnError: true, unstableOnWarning: false,
            projectRulePath: "${parseRulefile}", useProjectRule: true])

    step([$class: 'GitHubCommitNotifier', resultOnFailure: 'FAILURE', statusMessage: [content: 'Finished Jenkins gcc build']])

}
// *** End helper functions ***
