+++
date = "2020-11-24T09:45:13+01`:00"
title = "Get support / fill a bug report"
author = "Lars Bilke"
weight = 21

[menu]
  [menu.userguide]
    parent = "troubleshooting"
+++

If you encounter issues using or developing OGS let us know. To help you in the best possible way we need detailed information about the problem from you. In the best case we can reproduce the issue by ourselves and possibly fix it for you. Please provide the following information:

### Report contents

- OGS version (e.g. version number or git commit hash)
- Where did you obtained the program (e.g. self-compiled or as a binary download from the web site)
- Operating system / hardware information (e.g. Ubuntu 20.04, Intel i5-10600K, 16 GB Ram)
- Description of the issue:
  - Steps to reproduce the issue
  - **Full log output** as plain text (-file) (e.g. full compiler output including CMake run or output of a simulation run of OGS)! Please **do not** send screenshots of error messages! Instead copy and paste from your terminal application into the report!
  - Information about your environment (e.g. did you used Conan for third-party dependencies or did you installed them by yourself, which library versions)
  - Information about the OGS configuration (e.g. debug or release, serial or parallel, mfront?, ...)

Feel free to skip things from the list if you think they are not related to the issue but in general: the more information the better.


### Submit the report

Submit the report to one of the following channels (beforehand you may want to check these channels if somebody else already described your issue):

| Channel                                                                                       | Visibility                  | Requirement         |
| --------------------------------------------------------------------------------------------- | --------------------------- | ------------------- |
| Send an email to the [Service Desk](mailto:gitlab+ogs-ogs-120-issue-@opengeosys.org)[^desk]   | Private                     | -                   |
| Post on [Discourse](https://discourse.opengeosys.org)                                         | Public                      | A Discourse account |
| Create an issue on [OGS repository](https://gitlab.opengeosys.org/ogs/ogs/-/issues/new?issue) | Public or Private[^private] | A GitLab account    |

[^desk]: This creates a confidential (private) issue on our development repository. An OGS developer will get back to you soon. All communication happens via email. This is the *preferred* method *for external users*.
[^private]: Check the box *This issue is confidential ...* on issue creation.
