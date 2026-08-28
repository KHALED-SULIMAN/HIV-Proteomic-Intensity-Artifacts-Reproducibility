# Upload this repository to GitHub

## Recommended repository name

`HIV-Proteomic-Intensity-Artifacts-Reproducibility`

## Browser method

1. Create a new **Public** repository on GitHub.
2. Do not initialize it with a README, license, or `.gitignore`, because those files already exist here.
3. Extract the ZIP package locally.
4. Open the extracted repository folder and select **all files and folders inside it**.
5. On the empty GitHub repository page choose **Add file > Upload files** and drag the selected contents into the upload area.
6. Commit with the message: `Initial reproducibility release for Proteomes manuscript 4517853`.
7. Confirm that the repository root displays `README.md`, `scripts/`, `data/`, `results/`, `matlab/`, `REPRODUCIBILITY.md`, and `CITATION.cff`.

Do not upload the ZIP itself as the repository content, and do not create a nested `github_repo/` directory.

## Command-line method

From the extracted repository folder:

```bash
git init
git branch -M main
git add .
git commit -m "Initial reproducibility release for Proteomes manuscript 4517853"
git remote add origin https://github.com/YOUR_GITHUB_USERNAME/HIV-Proteomic-Intensity-Artifacts-Reproducibility.git
git push -u origin main
```

Replace `YOUR_GITHUB_USERNAME` with your actual GitHub username.
