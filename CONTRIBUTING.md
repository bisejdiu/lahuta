# CONTRIBUTING

Thank you for considering contributing to `Lahuta`! We welcome contributions from anyone, whether you're new to the project or you've been around a long time. This document provides some guidelines to help ensure a smooth collaboration process.

## Table of Contents

- [Getting Started](#getting-started)
- [Code of Conduct](#code-of-conduct)
- [How Can I Contribute?](#how-can-i-contribute)
  - [Reporting Bugs](#reporting-bugs)
  - [Suggesting Enhancements](#suggesting-enhancements)
  - [Your First Code Contribution](#your-first-code-contribution)
- [Pull Request Process](#pull-request-process)
- [Styleguides](#styleguides)
  - [Git Commit Messages](#git-commit-messages)
  - [Code Styleguide](#code-styleguide)
- [Additional Notes](#additional-notes)

## Getting Started

1. **Fork the repository**: Click on the 'Fork' button at the top right of this page and clone your forked repository to your local machine.

2. **Understanding the branches**:
   - `main`: This is the stable branch that end-users use. This is what gets packaged and distributed.
   - `dev`: This is where most development work happens. New features, bug fixes, and most other changes are merged here.
   - `docs`: This branch is specifically for documentation-related changes. Though the distinction between `docs` and `dev` is minimal, they are kept separate to organize the project better.
   - `feature/...`: When working on new features or major changes, branch off from `dev` using the `feature/` prefix. For instance, if you're working on a feature named "msa_indexing", your branch might be called `feature/msa_indexing`.

3. **Creating your branch**:
   - For features: Branch off from `dev` and name it `feature/your-feature-name`.
   - For documentation: You might prefer branching off from `docs`, or `dev` if the documentation changes are minor.
   - For bug fixes, general enhancements, or other: Typically, branch off from `dev`.

4. **Commit your changes**: Write a clear and concise commit message describing what you did.

5. **Push to your fork**: Once you're ready, push your branch to your fork on GitHub.

6. **Create a Pull Request**: 
   - If you're pushing a new feature or a major change, create a Pull Request to merge your `feature/...` branch into the `dev` branch of the main repository.
   - For documentation, create a Pull Request either to the `docs` or `dev` branch, depending on the nature of the changes.
   - In general, most contributions will be merged into the `dev` branch first before they eventually make their way to `main` during a release cycle.

## Code of Conduct

By participating in this project, you are expected to uphold our [Code of Conduct](CODE_OF_CONDUCT.md).

## How Can I Contribute?

### Reporting Bugs

Please check the issue tracker before creating a new bug report. If you're unable to find an open issue addressing the problem, open a new one.

### Suggesting Enhancements

If you have ideas or improvements, please send them our way! Just open an issue describing your suggestion.

## Pull Request Process

1. Make sure any install or build dependencies are removed before the end of the layer when doing a build.

2. Update the README.md with details of changes to the interface, this includes new environment variables, exposed ports, useful file locations, and container parameters.

3. Your Pull Request will be reviewed by the maintainers. They might ask for some changes or improvements before it gets merged.

## Styleguides

### Git Commit Messages

- Use the present tense ("Add feature" not "Added feature").
- Use the imperative mood ("Move cursor to..." not "Moves cursor to...").
- Limit lines to 120 characters or less.

### Code Styleguide

Specify linters ... 

## Additional Notes

Thank you for investing your time in Lahuta's development! We appreciate your effort and are always open to feedback and suggestions.
