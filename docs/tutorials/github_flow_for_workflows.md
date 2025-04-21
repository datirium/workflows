# Github Flow for contributing to workflows

<!-- - [Pre-reqeuisites](#prerequisites) -->
- [Dev-flow](#development-flow)

In order to properly contribute to the development of workflows, contributors should have their GitHub account linked with their  [SciDAP account](https://platform.scidap.com/login). This will allow the contributor to test their custom workflows within the SciDAP platform. 

This isn't required, and contributing to workflows can also be accomplished without scidap. This will make testing workflows that require upstream analyses more difficult though. It will also exclude SciDAP visualizations from the analyses

<!-- ---

## Prerequisites
- [Have a github account](https://docs.github.com/en/github/getting-started-with-github/signing-up-for-github)
- [Setup SciDAP account for contributing](../tutorials/scidap_for_contributors.md)
- [Sync custom workflow repo with your SciDAP account](../scidap-tutorials/sync_github_and_scidap.md)
- [Run workflows locally](./running_workflows_locally.md) -->

---

## Development flow

There are 2 options for developers when utilizing custom workflows on scidap. 

- [**push**](#push-based-flow) based updates
- and [**PR**](#pr-based-flow) based updates.

> **By default**, newly sync'd repo's will be configured to update by **push**

[Differences explained](//TODO)

---

### push based flow

[setup for push based dev-flow](./setup_push_flow.md)

Push based repo's will update workflow versions  in SciDAP automatically whenever changes are pushed to your remote repo.

After forking from datiriums workflows repo and syncing it with your SciDAP account, you're ready to push any workflow changes.

New workflow versions will need to be added to a project in SciDAP in order for them to be selectable when adding a new sample.

 

#### Things to remember about push based flow
- DON'T create or push any branches to the (synced) repo on github. 
- Every push to every branch will update modified workflows. This may cause unintended versioning effects and changes to workflows in development if you use multiple branches

---

### PR based flow

Instead of version updates occurring every push, updates are only triggered when a PR is accepted and merged into the main ("master") branch. 

This allows the owner (and approved github users) to review the updates made, and for those updates to trigger a single new version for modified workflows/tools.

#### Reconfiguring your custom repo to update versions on PR merges
Once you've forked datiriums workflows repo and sync'd it with your SciDAP account, you'll need to change some setting on the repo through GitHub.

STEPS:
1. open ```settings -> webhooks```
2. deselect the "push" event
3. select the "pull_request" event
4. save changes

// TODO: images