#!/usr/bin/env groovy
node {
    stage('Deploy S3 Info') {
        sh "curl 'https://svn.ufz.de:8443/job/OGS-6/job/Deploy/lastSuccessfulBuild/api/json?pretty=true&tree=actions" +
            "[artifacts[*]],url' -g --insecure -o lastSuccessfulArtifacts.json"
            archiveArtifacts 'lastSuccessfulArtifacts.json'
        step([$class: 'S3BucketPublisher', dontWaitForConcurrentBuildCompletion: true, entries:
            [[bucket: 'opengeosys', excludedFile: '', flatten: true, gzipFiles: false,
                managedArtifacts: false, noUploadOnFailure: true, selectedRegion: 'eu-central-1',
                sourceFile: "lastSuccessfulArtifacts.json", storageClass: 'STANDARD',
                uploadFromSlave: true, useServerSideEncryption: false]], profileName: 'S3 UFZ',
                userMetadata: [[key: 'Content-Type', value: 'Application/json']]])
}
