defaultDockerArgs = '-v /home/jenkins/.ccache:/usr/src/.ccache'
defaultCMakeOptions = '-DOGS_LIB_BOOST=System -DOGS_LIB_VTK=System -DOGS_DOWNLOAD_ADDITIONAL_CONTENT=ON'

node('visserv3')
{
    step([$class: 'GitHubSetCommitStatusBuilder', statusMessage: [content: 'Started Jenkins MSVC build']])

    stage 'Checkout'
    dir('ogs') { checkout scm }

    withEnv(['ARCH=msvc2013-x64', 'CMAKE_LIBRARY_SEARCH_PATH=C:\\libs\\$ARCH', 'QTDIR=C:\\libs\\qt\\4.8\\$ARCH', 'Path=$Path;$QTDIR\\bin;$CMAKE_LIBRARY_SEARCH_PATH\\bin']) {

        stage 'CLI'
        configure 'build', '', 'Visual Studio 12 Win64'
        build 'build', 'tests'
        build 'build', 'ctest'

        stage 'Data Explorer'
        configure 'build-de', '-DOGS_BUILD_GUI=ON -DOGS_BUILD_UTILS=ON -DOGS_BUILD_TESTS=OFF', 'Visual Studio 12 Win64'
        build 'build-de'

        if (env.BRANCH_NAME == 'master') {
            stage 'Package'
            build 'build', 'package'
            build 'build-de', 'package'
        }
    }
    if (env.BRANCH_NAME == 'master') {
        withEnv(['ARCH=msvc2013-x32', 'QTDIR=C:\\libs\\qt\\4.8\\$ARCH', 'Path=$Path;$QTDIR\\bin;C:\\Tools\\Conan\\conan']) {
            stage 'x32'
            configure 'build-32', '-DOGS_BUILD_GUI=ON -DOGS_BUILD_UTILS=ON -DOGS_BUILD_TESTS=OFF', 'Visual Studio 12', '-u -s build_type=Release -s compiler="Visual Studio" -s compiler.version=12 -s arch=x86'
            build 'build-32', 'package'
        }
    }

    stage 'Post'
    publishTestReports 'build/Testing/**/*.xml', 'build/Tests/testrunner.xml',
        'ogs/scripts/jenkins/msvc-log-parser.rules'

    archive 'build*/*.zip'

    step([$class: 'S3BucketPublisher', dontWaitForConcurrentBuildCompletion: false, entries: [[bucket: 'opengeosys', excludedFile: '', flatten: true, gzipFiles: false, managedArtifacts: true, noUploadOnFailure: true, selectedRegion: 'eu-central-1', sourceFile: 'build*/*.zip', storageClass: 'STANDARD', uploadFromSlave: false, useServerSideEncryption: false]], profileName: 'S3 UFZ', userMetadata: []])

}

// *** Helper functions ***
def configure(buildDir, cmakeOptions, generator, conan_args=null) {
    bat("""rd /S /Q ${buildDir}
           mkdir ${buildDir}""".stripIndent())
    if (conan_args != null)
        bat("""cd ${buildDir}
               conan install ../ogs ${conan_args}""".stripIndent())
    bat """cd ${buildDir}
           cmake ../ogs -G "${generator}" ${defaultCMakeOptions} ${cmakeOptions}"""
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
