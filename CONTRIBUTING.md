# Contribution Guide

This guide is intended for developers or administrators who want to contribute a new package, feature, or bugfix to JAWS. It assumes that you have at least some familiarity with Git VCS and Gitlab. The guide will show a few examples of contributing workflows and discuss the granularity of pull-requests (PRs). It will also discuss the tests your PR (Pull Request) must pass in order to be accepted into JAWS.

The changes one proposes in a PR should correspond to one completed feature/bugfix/extension/etc. One can create PRs with changes relevant to different ideas, however reviewing such PRs becomes tedious and error prone. If possible, try to follow the one-PR-one-package/feature rule.

JAWS uses a rough approximation of the [Git Flow branching model](http://nvie.com/posts/a-successful-git-branching-model/). The develop branch contains the latest contributions, and master is always tagged and points to the latest stable release. Therefore, when you send your request, make develop the destination branch on the JAWS repository.


## Continuous Integration

JAWS uses Gitlab CI for Continuous Integration testing. This means that every time you submit a pull request, a series of tests will be run. Your PR will not be accepted until it passes all of these tests. While you can certainly wait for the results of these tests after submitting a PR, we recommend that you run them locally to speed up the review process.


### Tests

#### Flake8 Tests

JAWS uses Flake8 to test for PEP 8 conformance. PEP 8 is a series of style guides for Python that provide suggestions for everything from variable naming to indentation. Your PR needs to comply with PEP 8 in order to be accepted.

Testing for PEP 8 compliance is easy:

```
make test

# or partial tests with
make test-client
make test-site
make test-central
```

Most of the error messages are straightforward, but if you don’t understand what they mean, please ask questions about them when you submit your PR.


#### Unit Tests

Unit tests ensure that core features are working as expected. If you make changes to JAWS, you should run the unit tests to make sure you didn’t break anything.


## Git Workflows

### Branching

The easiest way to contribute a pull request is to make all of your changes on new branches. Make sure your develop is up-to-date and create a new branch off of it:

```
git checkout develop
git pull upstream develop
git branch <descriptive_branch_name>
git checkout <descriptive_branch_name>
```

Here we assume that the local develop branch tracks the upstream develop branch of JAWS. This is not a requirement and you could also do the same with remote branches. But for some it is more convenient to have a local branch that tracks upstream.
Normally we prefer that commits pertaining to a package <package-name> have a message <package-name>: descriptive message. It is important to add descriptive messages so that others, who might be looking at your changes later (in a year or maybe two), would understand the rationale behind them.
Now, you can make your changes while keeping the develop branch pure. Edit a few files and commit them by running:

```
git add <files_to_be_part_of_the_commit>
git commit --message <descriptive_message_of_this_particular_commit>
```
 
Next, push it to your remote fork and create a PR:

```
git push origin <descriptive_branch_name> --set-upstream
```

Gitlab provides a tutorial on how to file a pull request. When you send the request, make develop the destination branch and make sure to use the default merge request template.
If you need this change immediately and don’t have time to wait for your PR to be merged, you can always work on this branch. But if you have multiple PRs, another option is to maintain a Frankenstein branch that combines all of your other branches:

```
git co develop
git branch <your_modified_develop_branch>
git checkout <your_modified_develop_branch>
git merge <descriptive_branch_name>
```

This can be done with each new PR you submit. Just make sure to keep this local branch up-to-date with upstream develop too.
Rebasing
Other developers are constantly making contributions to JAWS, possibly on the same files that your PR changed. If their PR is merged before yours, it can create a merge conflict. This means that your PR can no longer be automatically merged without a chance of breaking your changes. In this case, you will be asked to rebase on top of the latest upstream develop.
First, make sure your develop branch is up-to-date:

```
git checkout develop
git pull upstream develop
```

Now, we need to switch to the branch you submitted for your PR and rebase it on top of develop:

```
git checkout <descriptive_branch_name>
git rebase develop
```
 
Git will likely ask you to resolve conflicts. Edit the file that it says can’t be merged automatically and resolve the conflict. Then, run:

```
git add <file_that_could_not_be_merged>
git rebase --continue
```
 
You may have to repeat this process multiple times until all conflicts are resolved. Once this is done, simply force push your rebased branch to your remote fork:

```
git push --force origin <descriptive_branch_name>
```

