# RABDAM GUI Technology Plan

## Goal

Build a local desktop GUI for RABDAM3 that can run on Linux and Windows, display protein structures in 3D, and connect to the existing Python scientific backend.

## Proposed Stack

```text
Electron
+ React
+ TypeScript
+ Mol*
+ Python backend
```

## Role of Each Part

### Python

Python remains the scientific engine.

It handles:

- RABDAM3 calculations
- structure parsing and analysis
- BDamage/Bnet scoring
- database or CSV generation
- validation and benchmarking
- file input/output

### Electron

Electron provides the local desktop application shell.

It handles:

- opening the app window
- packaging the app for Linux and Windows
- providing a consistent Chromium/WebGL environment
- running the GUI locally without requiring the user to open a browser

### React

React builds the user interface.

It handles:

- buttons
- panels
- menus
- settings
- progress displays
- results tables
- interaction between the user, Python backend, and Mol* viewer

### TypeScript

TypeScript is used for safer frontend code.

It helps define structured data such as:

- structures
- chains
- residues
- atoms
- BDamage/Bnet scores
- validation results
- output tables

### Mol*

Mol* provides the built-in 3D protein structure viewer.

It handles:

- loading PDB/mmCIF structures
- rendering cartoon, sticks, atoms, ligands, and surfaces
- rotating, zooming, and selecting structures
- highlighting residues or atoms
- coloring structures by RABDAM/Bnet scores

## Basic Architecture

```text
User opens desktop app
        ↓
Electron launches local GUI
        ↓
React renders controls, tables, and Mol* viewer
        ↓
React sends commands to Python backend
        ↓
Python runs RABDAM3 calculations
        ↓
Results return to React
        ↓
React updates tables and tells Mol* what to highlight
```

## Example GUI Layout

```text
Left panel:
  - choose PDB/mmCIF file
  - choose output folder
  - analysis settings
  - run button

Center:
  - Mol* 3D protein viewer

Right panel:
  - selected residue/atom information
  - BDamage/Bnet scores
  - warnings or rejection reasons

Bottom:
  - progress log
  - results table
  - export buttons
```

## Packaging Plan

Use PyInstaller to package the Python backend as a local executable.

Use Electron Builder or Electron Forge to package the full desktop app.

Expected outputs:

```text
Linux:
  AppImage, .deb, or similar package

Windows:
  .exe installer or portable build
```

## Guiding Principle

Keep the scientific code separate from the GUI.

```text
Python backend:
  scientific correctness

React frontend:
  user interaction and display

Mol*:
  molecular visualization

Electron:
  local desktop packaging
```

This keeps RABDAM3 usable from the command line while allowing a polished GUI to be built around it.
