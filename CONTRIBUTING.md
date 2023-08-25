# CONTRIBUTING

Thank you for considering contributing to `Lahuta`! We welcome contributions from anyone, whether you're new to the project or you've been around a long time. This document provides some guidelines to help ensure a smooth collaboration process.

## Table of Contents

- [Getting Started](#getting-started)
- [Code of Conduct](#code-of-conduct)
- [How Can I Contribute?](#how-can-i-contribute)
- [Styleguides](#styleguides)

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

By participating in this project, you are expected to respect common codes of conduct associated with Python projects. 
As a specific example, consider adhering to the principles outlined in the [NumPy Code of Conduct](https://numpy.org/code-of-conduct/).

## How Can I Contribute?

### Reporting Bugs

Please check the issue tracker before creating a new bug report. If you're unable to find an open issue addressing the problem, open a new one.

### Suggesting Enhancements

If you have ideas or improvements, please send them our way! Just open an issue describing your suggestion.

## Styleguides

### Python Static Typing with `mypy`

`Lahuta` strives to maintain a modern and robust codebase, and as such, it leverages static typing. We use `mypy` in strict mode. When contributing:

- Ensure all your Python code is statically typed.
- Run `mypy` in strict mode to verify your types.
- Our Continuous Integration (CI) pipeline checks for any typing issues, and pull requests will fail if any discrepancies are detected.

### Linting with `ruff`

We use the `ruff` linter to ensure our code remains consistent and adheres to best practices.

- Before submitting your pull request, ensure you've run `ruff` to check for any linting issues.
- Address any warnings or errors `ruff` raises before submitting your contribution.

### Formatting with `black`

To maintain a consistent code format across the entire project, we employ `black`. 

- Before committing, make sure to format your code with `black`.
- This ensures all code adheres to the same formatting standards.
- We do not strictly enforce this. 

### Configuration with `pyproject.toml`

All project configuration, including for tools like `mypy`, `ruff`, and `black`, is centralized in `pyproject.toml`. This file serves as the single source of truth for project settings.

- If you need to alter configurations, make the changes in `pyproject.toml`.
- When adding dependencies or changing project metadata, update `pyproject.toml` accordingly.

Following these guidelines will ensure a smoother review process and will expedite the merging of your contributions.

## Additional Notes

Thank you for investing your time in Lahuta's development! We appreciate your effort and are always open to feedback and suggestions.
