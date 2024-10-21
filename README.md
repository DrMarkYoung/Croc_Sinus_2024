# Repository `<repository-name>`

<repository-description>

# License

This template is licensed under the [Creative Commons Attribution 4.0
International license](https://creativecommons.org/licenses/by/4.0/).
Please see the [license file](LICENSE.md).

When using this template, please don’t forget to:

- Adapt this license section to your own needs

- License your own content, and remember that [open is
  better](https://choosealicense.com/).

- Link to the [original
  license](https://creativecommons.org/licenses/by/4.0/) and give
  appropriate credit; please do so by including the following in the
  “License” section of the README.md file in your own project:

  > ## Attributions
  >
  > This project makes use of the
  > [rproj-template](https://github.com/DaniMori/rproj-template) Github
  > template created by [Daniel Morillo](https://github.com/DaniMori)
  > and licensed under the [Creative Commons Attribution 4.0
  > International
  > license](https://creativecommons.org/licenses/by/4.0/).

# How to use this template

1.  Click on the Green "Use this template" button on the top right corner.

2.  In the drop-down menu, click on "Create a new repository".

3.  In the "Create a new repository" page, fill in the field "Repository name"
    with a short name that represents the project you want to create.

4.  Optionally, add a free text description in the "Description" field, and
    change the visibility to "Private" if you do not want to make it freely
    accessible on the internet.

5.  Click on "Create repository".

6.  Follow the instructions in the next section to install (i.e. clone) your
    repository locally.

## Cloning your repository locally

To clone your repository locally, follow [these
instructions](https://book.cds101.com/using-rstudio-server-to-clone-a-github-repo-as-a-new-project.html#step---2)
(steps 2 through 7). This will create a local copy of ('clone') the GitHub
repository as an Rstudio project in the folder specified. The URL that must be
entered into the `Repository URL` text box will be of the form

```         
https://github.com/<username>/<project-name>
```

where `<username>` is your GitHub user name, and `<project-name>` is the value
you input before in the "Repository name" field in step 3.

**IMPORTANT:** It is totally unrecommended to clone a git repository inside a
cloud storage folder (e.g., Dropbox, OneDrive). Please note that GitHub serves
the purpose of backing up the repository, so no cloud storage is necessary.
Similarly, cloning the repository in a network folder may cause problems with
the `renv` environment (see below); do it at your own risk!

After cloning the repository, the Rstudio project will open automatically in the
Rstudio IDE. If it doesn't, or you want to return later to the project in
Rstudio, you can do so by double clicking on the file `rstudio_project.Rproj`
that has been created in the project folder when cloning the repository.

## Restoring the environment

The reproducible environment created by `renv` must be restored to install all
the packages this project needs to be built properly. If `renv` does not
initialize automatically (check the console for messages about this), you will
need to manually install the package first:

``` r
install.packages("renv")
```

Once it is successfully installed, use the "renv" -\> "Restore library..."
button in Rstudio's "Packages" tab to restore the environment. Alternatively,
you can type in the console:

``` r
renv::restore(prompt = FALSE)
```

## Creating your own "README.md" file

This template provides the "README.Rmd" Rmarkdown file to facilitate creating
your own "README.md" file. To use it, follow these instructions:

1.  Open file "README.Rmd" in the Rstudio editor

2.  Substitute the values of the objects `REPOSITORY_NAME` in line 12 and
    `GITHUB_USERNAME` in line 16 by your repository name and your GitHub user
    name, respectively.

3.  Add a free text repository description (line 22), to help other people
    understand the purpose of your project. This can be the same text you input
    previously in the "Description" field when creating your repository (see
    step 4 in the "How to use this template" main section), or a more detailed
    one.

4.  If necessary, [choose a license](https://choosealicense.com/) and update the
    "License" section. If you do so, remember to substitute [the license
    file](LICENSE.md) by your own license.

5.  In section "Attributions", add any additional attributions to software
    components you will use in your project.

6.  Knit the "README.Rmd" file.(\*)

7.  Commit (i.e., `git commit`) the changes to "README.md" and push (`git push`)
    them to your GitHub repository.

(\*) Please take into account that knitting "README.Rmd" will output a new
"README.md" file, overwriting this one.
