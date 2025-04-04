# Seurat Shiny Browser Package Cleanup Guide

## File Duplication Issue

Currently, there are duplicate files in the repository structure:

1. Original files in the top-level directories:
   - `server/`
   - `tabs/`
   - `global.R`
   - `server.R`
   - `ui.R`

2. Copies of these files in the package structure:
   - `inst/shiny-app/server/`
   - `inst/shiny-app/tabs/`
   - `inst/shiny-app/global.R`
   - `inst/shiny-app/server.R`
   - `inst/shiny-app/ui.R`

## Recommended Cleanup

For a clean package structure, you should keep only one copy of each file. Here's the recommended approach:

### Option 1: Remove Top-Level Files (Recommended for Package Distribution)

```bash
# Remove the original files and directories that are now duplicated
rm -rf server/ tabs/ global.R server.R ui.R
```

### Option 2: Keep Top-Level Files for Development (Alternative)

If you prefer to work directly with files at the top level during development:

1. Remove the duplicate files in the package structure:
   ```bash
   rm -rf inst/shiny-app/server/ inst/shiny-app/tabs/ inst/shiny-app/global.R inst/shiny-app/server.R inst/shiny-app/ui.R
   ```

2. Create symlinks from the package directory to the top-level files:
   ```bash
   mkdir -p inst/shiny-app
   ln -s "$(pwd)/server" "$(pwd)/inst/shiny-app/server"
   ln -s "$(pwd)/tabs" "$(pwd)/inst/shiny-app/tabs"
   ln -s "$(pwd)/global.R" "$(pwd)/inst/shiny-app/global.R"
   ln -s "$(pwd)/server.R" "$(pwd)/inst/shiny-app/server.R"
   ln -s "$(pwd)/ui.R" "$(pwd)/inst/shiny-app/ui.R"
   ```

## Final Structure

After cleanup, your package should have one of these structures:

### Package Distribution Structure
```
seurat-shiny-browser/
├── DESCRIPTION
├── LICENSE
├── NAMESPACE
├── R/                  # Package functions
├── inst/
│   └── shiny-app/      # All app files live here
│       ├── server/     # Server components
│       ├── tabs/       # UI tab components
│       ├── app.R       # Entry point
│       ├── global.R    # Global settings
│       ├── server.R    # Server logic
│       └── ui.R        # UI definition
├── man/                # Documentation
└── README.md
```

### Development Structure with Symlinks
```
seurat-shiny-browser/
├── DESCRIPTION
├── LICENSE
├── NAMESPACE
├── R/                  # Package functions
├── server/             # Server components (original files)
├── tabs/               # UI tab components (original files)
├── global.R            # Global settings (original file)
├── server.R            # Server logic (original file)
├── ui.R                # UI definition (original file)
├── inst/
│   └── shiny-app/      # Contains symlinks to original files
├── man/                # Documentation
└── README.md
```

## Notes

- The `app.R` file in `inst/shiny-app/` is needed for the package to run properly
- The files in `R/` directory (package functions) are separate from the Shiny app files
- Only one copy of each file should exist to avoid confusion and maintenance issues
- R package conventions expect the actual app files to be in `inst/shiny-app/`