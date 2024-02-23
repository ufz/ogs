+++
date = "2018-02-26T11:00:13+01:00"
title = "Code Reviews"
author = "Lars Bilke"
weight = 1015

[menu]
  [menu.devguide]
    parent = "development-workflows"
+++

Once you have submitted a [merge request]({{< ref "setup-fork.md/#optional-working-on-a-new-feature" >}}) the code review process kicks in which can be summarized as:

- A core developer picks your MR for review.
- The core developer checks your code in respect to e.g. style and correctness.
- The core developer may also give hints on how to improve the code regarding to e.g. readability or performance.
- You iterate (make modifications) on the code.
- The core developer again checks the code and finally approves the changes.
- If also the [CI]({{< ref "continuous-integration.md" >}}) is happy (all checks and tests did pass) the MR is merged.

For more information on merge requests see the [GitLab documentation](https://docs.gitlab.com/ee/user/project/merge_requests).

### How to request a code review

The default label upon creating a merge request is `workflow::in development`.

To request a code review set the merge request label to `workflow::please review`!

A merge request reviewer will now have look. If changes are requested the reviewer may change the label to `workflow::needs update` indicating that the merge request author is now back in charge again.

To indicate that you are currently not working on a merge request set it to `workflow::paused` (disables also CI – there is also a dedicated label for disabling CI: `ci skip`).

The following diagram summarises the workflow:

```mermaid
graph TB
  subgraph cr[Code Review]
  review([please review])
  review-->|reviewer requests update|update([needs update])
  update-->|author fixed issues and requests review again|review
  end

  mr(Create a merge request)-->dev
  subgraph ds[Development State]
  dev([in development])-->|author requests review|review
  dev-->blocked([blocked])
  blocked-->dev
  dev<-->paused([paused])
  paused-->dev
  end

  update-->dev

  style dev stroke:#418bc9,stroke-width:2px
  style paused stroke:#418bc9,stroke-width:2px
  style blocked stroke:#ad4363,stroke-width:2px
  style review stroke:#5bb85c,stroke-width:2px
  style update stroke:#7f8b8d,stroke-width:2px
  style cr stroke:#5bb85c,stroke-width:4px
```

### How to indicate that the MR is ready for merge

- Set the merge request label to `workflow::please review`.
- Remove the [`Draft`-flag](https://docs.gitlab.com/ee/user/project/merge_requests/drafts.html#mark-merge-requests-as-ready) (if you have set it previously).

### How to checkout a MR from another developer locally

On the merge request page in the first box which contains information on the MR author and branch name there is button labelled `Check out branch` which will show you instructions on how to locally checkout this MR:

![The checkout branch button in the GitLab web interface](checkout-branch.png)
