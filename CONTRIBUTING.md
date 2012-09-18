# How to contribute



## Getting Started

* Make sure you have a [GitHub account](https://github.com/signup/free)
* Fork the repository on GitHub

## Making Changes

* Create a topic branch from where you want to base your work.
  * This is usually the master branch.
  * Only target other branches if you are certain your fix must be on that
    branch.
  * To quickly create a topic branch based on master; `git branch
    my_topic_branch master` then checkout the new branch with `git
    checkout my_topic_branch`.  Please avoid working directly on the
    `master` branch.
* Make commits of logical units.
* Make sure your code conforms to the [styleguide][styleguide].
* Check for unnecessary whitespace with `git diff --check` before committing.
* Make sure your commit messages are in the proper format.

````
    One sentence summary.

    Detailed description of the commit ...
    ....
````

* TODO: Test your changes

## Submitting Changes

* Push your changes to a topic branch in your fork of the repository.
* Submit a pull request to the main repository.

# Additional Resources

* [General GitHub documentation](http://help.github.com/)
* [GitHub pull request documentation](http://help.github.com/send-pull-requests/)
* [OGS Jenkins-CI server](https://svn.ufz.de/hudson/job/OGS-6/)
* [OGS Styleguide][styleguide]

[styleguide]: http://ufz.github.com/styleguide/cppguide.xml