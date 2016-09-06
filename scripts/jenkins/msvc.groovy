defaultDockerArgs = '-v /home/jenkins/.ccache:/usr/src/.ccache'
defaultCMakeOptions = '-DOGS_LIB_BOOST=System -DOGS_LIB_VTK=System -DOGS_DOWNLOAD_ADDITIONAL_CONTENT=ON'
env32 = ['ARCH=msvc2013-x32', 'QTDIR=C:\\libs\\qt\\4.8\\$ARCH', 'Path=$Path;$QTDIR\\bin;C:\\Tools\\Conan\\conan']
env64 = ['ARCH=msvc2013-x64', 'CMAKE_LIBRARY_SEARCH_PATH=C:\\libs\\$ARCH', 'QTDIR=C:\\libs\\qt\\4.8\\$ARCH', 'Path=$Path;$QTDIR\\bin;$CMAKE_LIBRARY_SEARCH_PATH\\bin']

node('visserv3')
{
    step([$class: 'GitHubSetCommitStatusBuilder', statusMessage: [content: 'Started Jenkins MSVC build']])

    stage 'Checkout (Win)'
    dir('ogs') { checkout scm }

    withEnv(env64) {
        stage 'Configure (Win)'
        configure 'build', '', 'Visual Studio 12 Win64'

        stage 'CLI (Win)'
        build 'build'

        stage 'Test (Win)'
        build 'build', 'tests'

        stage 'Data Explorer (Win)'
        configure 'build-de', '-DOGS_BUILD_GUI=ON -DOGS_BUILD_UTILS=ON -DOGS_BUILD_TESTS=OFF', 'Visual Studio 12 Win64'
        build 'build-de'
    }

    if (env.BRANCH_NAME == 'master' || env.BRANCH_NAME.contains('release') ) {
        stage 'Release (Win)'
        withEnv(env64) {
            build 'build', 'package'
            build 'build-de', 'package'
        }
        withEnv(env32) {
            configure 'build-32', '-DOGS_32_BIT=ON -DOGS_BUILD_GUI=ON -DOGS_BUILD_UTILS=ON -DOGS_BUILD_TESTS=OFF', 'Visual Studio 12', '-u -s build_type=Release -s compiler="Visual Studio" -s compiler.version=12 -s arch=x86'
            build 'build-32', 'package'
        }
        deploy 'build/*.zip,build-de/*.zip,build-32/*.zip'
    }

    stage 'Post (Win)'
    publishTestReports 'build/Testing/**/*.xml', 'build/Tests/testrunner.xml',
        'ogs/scripts/jenkins/msvc-log-parser.rules'

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

def deploy(files) {
    archive "${files}"
    step([$class: 'S3BucketPublisher', dontWaitForConcurrentBuildCompletion: true, entries: [[bucket: 'opengeosys', excludedFile: '', flatten: true, gzipFiles: false, managedArtifacts: true, noUploadOnFailure: true, selectedRegion: 'eu-central-1', sourceFile: "${files}", storageClass: 'STANDARD', uploadFromSlave: true, useServerSideEncryption: false]], profileName: 'S3 UFZ', userMetadata: []])
}

def publishTestReports(ctestPattern, gtestPattern, parseRulefile) {
    step([$class: 'XUnitPublisher', testTimeMargin: '3000', thresholdMode: 1,
        thresholds: [
            [$class: 'FailedThreshold', failureNewThreshold: '', failureThreshold: '', unstableNewThreshold: '', unstableThreshold: ''],
            [$class: 'SkippedThreshold', failureNewThreshold: '', failureThreshold: '', unstableNewThreshold: '', unstableThreshold: '']],
        tools: [
            [$class: 'GoogleTestType', deleteOutputFiles: true, failIfNotNew: true, pattern: "${gtestPattern}", skipNoTestFiles: false, stopProcessingIfError: true]]
    ])

    step([$class: 'LogParserPublisher', failBuildOnError: true, unstableOnWarning: false,
            projectRulePath: "${parseRulefile}", useProjectRule: true])

    step([$class: 'GitHubCommitNotifier', resultOnFailure: 'FAILURE', statusMessage: [content: 'Finished Jenkins gcc build']])

}
// *** End helper functions ***
