# Step-by-step record: creating the `pasture-metagenomics` GitHub repository

This document records the steps used to prepare, initialize, and push the metagenomics project repository to GitHub.

Repository name:

```text
pasture-metagenomics
```

GitHub URL:

```text
https://github.com/Sereyboth94/pasture-metagenomics
```

Local repository folder:

```text
C:\Users\soths\OneDrive - Lincoln University\Writing\NatCom\Github
```

---

## 1. Prepare the local project folder

The final repository folder contained:

```text
Github/
├── bioinformatics/
├── docs/
├── function/
├── location/
├── manuscript_figures/
├── taxonomy/
├── .gitignore
└── README.md
```

The three main standalone R workflows were:

```text
location/run_location_figure1_pipeline.R
taxonomy/scripts/run_taxonomy_pipeline.R
function/scripts/run_function_pipeline.R
```

---

## 2. Check that the workflows run locally

From RStudio:

```r
setwd("C:/Users/soths/OneDrive - Lincoln University/Writing/NatCom/Github")

source("taxonomy/scripts/run_taxonomy_pipeline.R")
source("function/scripts/run_function_pipeline.R")
source("location/run_location_figure1_pipeline.R")
```

Expected output folders:

```text
manuscript_figures/Main_Figures/
manuscript_figures/Supplementary_Taxonomy_Figures/
manuscript_figures/Supplementary_Function_Figures/
```

---

## 3. Check for large files before GitHub upload

From PowerShell:

```powershell
cd "C:\Users\soths\OneDrive - Lincoln University\Writing\NatCom\Github"

Get-ChildItem -Recurse |
Sort-Object Length -Descending |
Select-Object -First 20 Name,@{N="MB";E={[math]::Round($_.Length/1MB,2)}}
```

Large raw sequencing files should not be committed to GitHub.

Avoid committing:

```text
*.fastq.gz
*.fq.gz
*.tar
*.tar.gz
large databases
large intermediate files
```

---

## 4. Create a new private GitHub repository

On GitHub:

```text
https://github.com/new
```

Settings used:

```text
Repository name: pasture-metagenomics
Visibility: Private
README: optional, but ideally OFF if a local README already exists
.gitignore: No .gitignore
License: optional
```

Repository created at:

```text
https://github.com/Sereyboth94/pasture-metagenomics
```

---

## 5. Initialize Git locally

From PowerShell:

```powershell
cd "C:\Users\soths\OneDrive - Lincoln University\Writing\NatCom\Github"

git init
git branch -M main
```

---

## 6. Stage and commit files

```powershell
git add .
git commit -m "Initial release of pasture metagenomics workflows"
```

---

## 7. Connect local repository to GitHub

```powershell
git remote add origin https://github.com/Sereyboth94/pasture-metagenomics.git
```

Check the remote:

```powershell
git remote -v
```

Expected:

```text
origin  https://github.com/Sereyboth94/pasture-metagenomics.git (fetch)
origin  https://github.com/Sereyboth94/pasture-metagenomics.git (push)
```

---

## 8. Push to GitHub

```powershell
git push -u origin main
```

---

## 9. Fix remote conflict if GitHub already has a README

If push fails with:

```text
! [rejected] main -> main (fetch first)
```

then the GitHub repository already contains a commit, usually a README.

Option A: merge remote history

```powershell
git pull origin main --allow-unrelated-histories
git push -u origin main
```

If there is a README conflict and the local README should be kept:

```powershell
git checkout --ours README.md
git add README.md
git commit -m "Merge remote repository setup"
git push -u origin main
```

Option B: overwrite the fresh remote repository

Use this only if the GitHub repository has no important files:

```powershell
git push -u origin main --force
```

---

## 10. Check repository status

```powershell
git status
```

Expected:

```text
nothing to commit, working tree clean
```

Check repository size:

```powershell
git count-objects -vH
```

---

## 11. Verify GitHub contents

After pushing, check that GitHub contains:

```text
README.md
.gitignore
bioinformatics/
docs/
function/
location/
manuscript_figures/
taxonomy/
```

---

## 12. Recommended repository release workflow

Before publication:

1. Keep the repository private while checking scripts and documentation.
2. Confirm that all three workflows run.
3. Confirm final manuscript figures are present.
4. Add SRA/BioProject accession once available.
5. Add Zenodo DOI after archiving.
6. Make repository public.
7. Create a GitHub release.

Suggested release tag:

```text
v1.0.0-natcom-submission
```

---

## 13. Zenodo archiving

After the repository is public and ready:

1. Log in to Zenodo.
2. Connect GitHub to Zenodo.
3. Enable the `pasture-metagenomics` repository.
4. Create a GitHub release.
5. Zenodo will archive the release and generate a DOI.
6. Add the DOI to `README.md`.

---

## 14. Useful update commands after future edits

After editing files:

```powershell
cd "C:\Users\soths\OneDrive - Lincoln University\Writing\NatCom\Github"

git status
git add .
git commit -m "Update documentation and workflows"
git push
```

---

## 15. Suggested citation placeholder

```text
Soth S. et al. Pasture metagenomics workflows and manuscript figures.
GitHub repository: https://github.com/Sereyboth94/pasture-metagenomics
Zenodo DOI: TO BE ADDED
```
