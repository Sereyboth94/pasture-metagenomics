# Repository setup and biomass update for `pasture-metagenomics`

This document records the original repository setup and the current steps for adding the plant biomass workflow to the existing private GitHub repository. Sections 4–9 describe the original setup; use Section 14 for this update.

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

After the biomass update, the repository folder should contain:

```text
Github/
├── Biomass/
│   ├── Plant_Biomass.xlsx
│   ├── Plant_Biomass_ANOVA.R
│   ├── run_location_figure1_pipeline.R
│   ├── README.md
│   └── .gitignore
├── bioinformatics/
├── docs/
├── function/
├── location/
├── manuscript_figures/
├── taxonomy/
├── .gitignore
└── README.md
```

The analysis has three main execution steps. The biomass directory also contains a standalone analysis script and the plot-level input workbook:

```text
taxonomy/scripts/run_taxonomy_pipeline.R
function/scripts/run_function_pipeline.R
Biomass/run_location_figure1_pipeline.R
Biomass/Plant_Biomass_ANOVA.R
Biomass/Plant_Biomass.xlsx
```

`Biomass/run_location_figure1_pipeline.R` reruns `Plant_Biomass_ANOVA.R` before combining biomass panel A with the taxonomy and pathway RDS panels. The old `location/run_location_figure1_pipeline.R` path should not be used for the new Figure 1 workflow.

---

## 2. Check that the workflows run locally

From RStudio:

```r
setwd("C:/Users/soths/OneDrive - Lincoln University/Writing/NatCom/Github")

source("taxonomy/scripts/run_taxonomy_pipeline.R")
source("function/scripts/run_function_pipeline.R")
source("Biomass/run_location_figure1_pipeline.R")
```

Expected output folders:

```text
Biomass/Plant_Biomass_ANOVA_results/
location/results/
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
Biomass/Plant_Biomass_ANOVA_results/
*.backup
large databases
large intermediate files
```

The small `Biomass/Plant_Biomass.xlsx` workbook contains plot-level measurements needed to reproduce panel A. Check that the project team permits its inclusion before making the repository public.

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

This conflict section applies to initial repository creation. For the biomass update, preserve remote history and use the targeted commands in Section 14; do not force-push an established repository.

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
Biomass/
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
2. Confirm the taxonomy, function and biomass/Figure 1 workflows run in the order shown in Section 2.
3. Confirm final manuscript figures are present and the biomass workbook is approved for sharing.
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

## 14. Add Biomass to the existing repository

1. Before copying files, update the local checkout from PowerShell while the working tree is clean:

```powershell
cd "C:\Users\soths\OneDrive - Lincoln University\Writing\NatCom\Github"
git status
git pull --ff-only origin main
```

2. Place the supplied `Biomass/` directory directly inside the local `Github` repository. Keep the capital `B` in the folder name. It contains `Plant_Biomass.xlsx`, `Plant_Biomass_ANOVA.R`, `run_location_figure1_pipeline.R`, a folder README and a local `.gitignore` for generated results.
3. Add a short section to the root `README.md` linking to `Biomass/README.md` and showing the run order in Section 2. Do not replace unrelated README content.
4. Review and push only these changes to the existing `main` branch:

```powershell
git add Biomass README.md
git diff --cached --stat
git diff --cached -- Biomass/README.md README.md
git commit -m "Add plant biomass analysis and Figure 1 assembly"
git push origin main
git status
```

If the initial `git status` shows other local changes, review or commit those before pulling. The commit commands stage only the biomass folder and the root README. Generated `Biomass/Plant_Biomass_ANOVA_results/` files are ignored by `Biomass/.gitignore`.

5. Verify the GitHub repository displays the five tracked files in `Biomass/` and the root README link. Run the three scripts in Section 2 locally before treating the figure as reproduced. This guide updates the repository procedure; it does not itself push files to GitHub.

---

## 15. Suggested citation placeholder

```text
Soth S. et al. Pasture metagenomics workflows and manuscript figures.
GitHub repository: https://github.com/Sereyboth94/pasture-metagenomics
Zenodo DOI: TO BE ADDED
```
