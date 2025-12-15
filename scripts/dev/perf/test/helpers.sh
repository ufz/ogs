# Prints a command on stderr and runs it.
RUN()
{
    echo "+ $*" >&2
    "$@"
}

# Runs a command and prints it on stderr if it fails.
CHECK()
{
    if "$@"
    then
        :
    else
        echo "check $* failed with status $?" >&2
        false
    fi
}

# Runs a command that is expected to fail and prints it on stderr if it succeeds.
CHECK_NOT()
{
    if "$@"
    then
        echo "check $* succeeded but was expected to fail" >&2
        false
    fi
}

# Same as RUN() but with a headline for a new test case
TEST()
{
    echo "== TEST $*" >&2
    "$@"
    echo >&2
}
