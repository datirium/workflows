# Push-based Developer flow



Steps for modifying or developing new workflows 

- [Push-based Developer flow](#push-based-developer-flow)
  - [Local Development](#local-development)
  - [Local Testing](#local-testing)
  - [Custom testing](#custom-testing)
  - [Create PR](#create-pr)

--

## Local Development

In order to start developing locally, you should stage git for the changes you're about to make


1. make sure your workflows are up-to-date
```
git checkout master
git pull origin master
```
2. create a new branch for your development. if there is a jira ticket associated, use that ticket name in the branch  
```git checkout -b SDAP-####-small-summary```


---

## Local Testing

Testing will change based on what was modified.

<!-- TODO -->

<!-- - [docker changes](../testing/docker_changes.md)
- [tool changes](../testing/tool_changes.md)
- [workflow changes](../testing/workflow_changes.md) -->

---

## Custom testing

Once local testing is successful, testing within scidap is appropriate.

Make sure all local changes are added and committed.

If we were developing locally on a branch called "SDAP-5000-new-cwl", then we would squash merge by running:  
```
git checkout main
git merge --squash SDAP-500-new-cwl
```

If there are any merge conflicts, run ```git status``` and then ctrl/cmd click the files with conflicts.
Any conflicts here would likely stem from the git commit history changes, so its likely that all conflicts can be handled by "accepting incoming changes" through vscodes conflict manager.

Once there are no conflicts, make sure the merge is commited to the "main" branch, and run:  
``` git push CUSTOM_REPO main```, where "CUSTOM_REPO" is the name of your reference for the remote repo synced with scidap.

Pushes will trigger automatic version updated, and samples can be updated or created through the platform.


---

## Create PR

If all is good in the custom testing, then you can create a PR from the branch you did local development on, from your forked repo (via github UI).


If custom testing requires changes and you commit directly to the "main" branch to make these changes (which is valid), then you will need to include these changes in the PR branch.

Lets say the branch for local dev/testing was called "SDAP-5000-new-cwl". Run:  
```
git checkout SDAP-5000-new-cwl
git merge --squash main 
```
And then you will likely need to fix merge conflicts.
Its safe to assume all conflicts can be resolved by accepting incoming changes, so (while in the conflict) run: ```git checkout --theirs .```

Now that there is no conflict, you can push this branch to your forked repo:
```git push FORKED_REPO SDAP-5000-new-cwl```

When you make the PR, if you find any files changed that shouldn't be changed, you can reset that file to whats on global by running: ```git checkout master -- path/to/file.txt```
> NOTE: you'll have to add, commit, and push those changes for the PR to get that update



