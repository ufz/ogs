defaultDockerArgs = '-v /home/jenkins/.ccache:/usr/src/.ccache'
defaultCMakeOptions = '-DOGS_LIB_BOOST=System -DOGS_LIB_VTK=System'

node('visserv3')
{
    step([$class: 'GitHubSetCommitStatusBuilder', statusMessage: [content: 'Started Jenkins MSVC build']])

    stage 'Checkout'
    dir('ogs') { checkout scm }

    withEnv(['ARCH=msvc2013-x64', 'CMAKE_LIBRARY_SEARCH_PATH=C:\\libs\\$ARCH', 'QTDIR=C:\\libs\\qt\\4.8\\$ARCH', 'Path=$Path;$QTDIR\\bin;$CMAKE_LIBRARY_SEARCH_PATH\\bin']) {
        stage 'Configure & Build & Test'
        bat('''rd /S /Q build
               mkdir build
               cd build
               cmake -G"Visual Studio 12 Win64" -DOGS_BUILD_GUI=ON -DOGS_BUILD_CLI=OFF -DOGS_DOWNLOAD_ADDITIONAL_CONTENT=ON ..\\ogs
               cmake --build . --config Release --target tests
               cmake --build . --config Release --target ctest'''.stripIndent())
        if (env.BRANCH_NAME == 'master')
            bat('''cd build
                   cmake --build . --config Release --target package
                   rename ogs*Windows*.zip ogs-data_explorer-x64.zip'''.stripIndent())
    }

    stage 'Post'
    publishTestReports 'build/Testing/**/*.xml', 'build/Tests/testrunner.xml',
        'ogs/scripts/jenkins/msvc-log-parser.rules'

    archive 'build*/*.zip'
}

// *** Helper functions ***
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
