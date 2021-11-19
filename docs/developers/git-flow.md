# Git Flow

We use the `develop/master` variation of the [OneFlow](https://www.endoflineblog.com/oneflow-a-git-branching-model-and-workflow) git flow

## Add New Features
We use feature (topic) branches to implement new features

=== "Internal Developer"
    You are an internal developer if you have writing permissions to the repository.
    
    Most feature branches are never pushed to the repo, only do so if you expect that its development will take days (to avoid losing your work if you computer is damaged). Otherwise follow the following instructions to locally rebase your feature branch into `develop` and push those rebased changes online.

    **Starting your feature branch**

    1. Pull the latest develop 
    ```bash
    git checkout develop
    git pull
    ```
    1. Create your feature branch
    ```bash
    git checkout -b feature/feature1
    ```
    1. Add, modify or delete the necessary files to add your new feature
    1. Update the [change log](../../change-log) (`docs/change-log.md`)
    2. Stage and commit your changes using VS Code git GUI or the following commands
    ```bash
    git add modified-file1 modified-file2
    git commit -m "Add my new feature" # use a concise description
    ```

    **Merging back your feature branch**

    If your changes took time to be implemented it is possible that there are new commits in our `develop` branch, so we need to rebase your feature branch.

    1. Fetch the latest changes to develop
    ```bash
    git fetch origin develop
    ```

    1. Rebase your feature branch
    ```bash
    git checkout feature/feature1
    git rebase -i develop
    ```

    1. Integrate your new feature to `develop`
    ```bash
    git checkout develop
    git merge --no-ff feature/feature1 # (use the default merge message)
    git push origin develop
    git branch -d feature/feature1
    ```

=== "External Developer"
    You are an external developer if you do NOT have writing permissions to the repository.

    **Starting your feature branch**

    1. Fork and clone our repository on Github
    1. Switch to the latest develop 
    ```bash
    git checkout develop
    ```
    1. Create your feature branch
    ```bash
    git checkout -b feature/external-test
    ```
    1. Add, modify or delete the necessary files to add your new feature
    2. Stage and commit your changes using VS Code git GUI or the following commands
    ```bash
    git add modified-file1 modified-file2
    git commit -m "Add my new feature" # use a concise description
    ```
    
    **Merging back your feature branch**
    
    If your changes took time to be implemented, it is possible that there are new commits in our `develop` branch, so we need to rebase your feature branch.

    1. Add our repo as another `remote`
    ```bash
    git remote add upstream https://github.com/carissalow/rapids/
    ```

    1. Fetch the latest changes to develop
    ```bash
    git fetch upstream develop 
    ```
    
    1. Rebase your feature branch
    ```bash
    git checkout feature/external-test
    git rebase -i develop
    ```
    
    1. Push your feature branch online
    ```bash
    git push --set-upstream origin feature/external-test
    ```
    
    1. Open a pull request to the `develop` branch using Github's GUI

## Release a New Version

1. Pull the latest develop 
```bash
git checkout develop
git pull
```
1. Create a new release branch
```bash
git describe --abbrev=0 --tags # Bump the release (0.1.0 to 0.2.0 => NEW_HOTFIX)
git checkout -b release/v[NEW_RELEASE] develop
```
1. Add new tag
```bash
git tag v[NEW_RELEASE]
```
1. Merge and push the release branch
```bash
git checkout develop
git merge release/v[NEW_RELEASE]
git push --tags origin develop
git branch -d release/v[NEW_RELEASE]
```
1. Fast-forward master
```
git checkout master
git merge --ff-only develop
git push # Unlock the master branch before merging
```
1. Release happens automatically after passing the tests

## Release a Hotfix
1. Pull the latest master
```bash
git checkout master
git pull
```
1. Start a hotfix branch
```bash
git describe --abbrev=0 --tags # Bump the hotfix (0.1.0 to 0.1.1 => NEW_HOTFIX)
git checkout -b hotfix/v[NEW_HOTFIX] master
```
1. Fix whatever needs to be fixed
1. Update the change log
1. Tag and merge the hotfix
```bash
git tag v[NEW_HOTFIX]
git checkout develop
git merge hotfix/v[NEW_HOTFIX]
git push --tags origin develop
git branch -d hotfix/v[NEW_HOTFIX]
```
1. Fast-forward master
```
git checkout master
git merge --ff-only v[NEW_HOTFIX]
git push # Unlock the master branch before merging
```
1. Release happens automatically after passing the tests
