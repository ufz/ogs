def isRelease () {
    if (env.BRANCH_NAME == 'master' || env.BRANCH_NAME.contains('release'))
        return true;
    return false;
}

return this;
