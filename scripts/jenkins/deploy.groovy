#!/usr/bin/env groovy
node {
    stage('Deploy to S3') {
        deleteDir()
        step([$class: 'CopyArtifact',
            fingerprintArtifacts: true, flatten: true,
            projectName: 'OGS-6/ufz/master',
            selector: [$class: 'LastCompletedBuildSelector']])
            if (gitTag == "")
                s3upload('*')
            else
                s3upload('*', "opengeosys/ogs6-releases/${gitTag}")
        build job: 'OGS-6/Deploy-Post', wait: false
    }
}

def s3upload(files, bucket = null) {
    def managed = false
    if (bucket == null) {
        managed = true
        bucket = 'opengeosys'
    }

    step([$class: 'S3BucketPublisher',
        dontWaitForConcurrentBuildCompletion: true, entries:
        [[bucket: "${bucket}", excludedFile: '', flatten: true, gzipFiles: false,
            managedArtifacts: managed, noUploadOnFailure: true, selectedRegion: 'eu-central-1',
            sourceFile: "${files}", storageClass: 'STANDARD', uploadFromSlave: true,
            useServerSideEncryption: false]], profileName: 'S3 UFZ', userMetadata: []])
}
