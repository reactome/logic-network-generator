# Contributing to Logic Network Generator

Thank you for your interest in contributing! This document provides guidelines for contributing to the project.

## Getting Started

### Prerequisites

- Python 3.9+
- Poetry
- Docker (for Neo4j database)
- Git

### Development Setup

1. **Fork and clone the repository**
   ```bash
   git clone https://github.com/YOUR_USERNAME/logic-network-generator.git
   cd logic-network-generator
   ```

2. **Install dependencies**
   ```bash
   poetry install
   ```

3. **Start Neo4j database** (for integration tests)
   ```bash
   docker run -p 7474:7474 -p 7687:7687 \
     -e NEO4J_dbms_memory_heap_maxSize=8g \
     public.ecr.aws/reactome/graphdb:Release94
   ```

4. **Install pre-commit hooks**
   ```bash
   poetry run pre-commit install
   ```

## Development Workflow

### 1. Create a Branch

Create a feature branch from `main`:
```bash
git checkout -b feature/your-feature-name
# or
git checkout -b fix/your-bug-fix
```

Branch naming conventions:
- `feature/` - New features
- `fix/` - Bug fixes
- `docs/` - Documentation updates
- `refactor/` - Code refactoring
- `test/` - Test improvements

### 2. Make Changes

- Write clean, readable code
- Follow existing code style and patterns
- Add type hints to all functions
- Write docstrings for public functions and classes
- Keep commits atomic and focused

### 3. Write Tests

- **Unit tests** are required for all new features and bug fixes
- Add tests to the appropriate file in `tests/`
- Ensure tests pass locally before pushing

Run unit tests (fast, no database required):
```bash
poetry run pytest tests/ -v -m "not database"
```

Run all tests including integration tests (requires Neo4j):
```bash
poetry run pytest tests/ -v
```

### 4. Code Quality

Before committing, ensure your code passes all quality checks:

**Run linter:**
```bash
poetry run ruff check src/
poetry run ruff format src/
```

**Run type checker (optional but recommended):**
```bash
poetry run mypy --ignore-missing-imports src/
```

**Or use pre-commit to run all checks:**
```bash
poetry run pre-commit run --all-files
```

### 5. Commit Changes

Write clear, descriptive commit messages:
```bash
git add .
git commit -m "Add feature: brief description

Longer explanation of what changed and why (if needed).

Fixes #123"
```

Commit message guidelines:
- Use present tense ("Add feature" not "Added feature")
- First line should be 50 characters or less
- Reference issue numbers when applicable

### 6. Push and Create Pull Request

```bash
git push origin feature/your-feature-name
```

Then create a pull request on GitHub:
- Fill out the PR template completely
- Link related issues
- Describe what was changed and why
- Include screenshots or output if relevant

## Code Style Guidelines

### Python Style

We use Ruff for linting and formatting:
- Maximum line length: 100 characters
- Use type hints for function signatures
- Follow PEP 8 naming conventions
- Use descriptive variable names

### Documentation Style

- Use Google-style docstrings
- Document all public functions, classes, and modules
- Include examples in docstrings when helpful
- Keep README and documentation up to date

Example docstring:
```python
def generate_logic_network(pathway_id: str) -> pd.DataFrame:
    """Generate a logic network for a Reactome pathway.

    Args:
        pathway_id: Reactome pathway database identifier

    Returns:
        DataFrame containing the logic network edges

    Raises:
        ValueError: If pathway_id is invalid
        ConnectionError: If cannot connect to Neo4j

    Example:
        >>> network = generate_logic_network("69620")
        >>> print(len(network))
        1234
    """
```

### Test Style

- Test file names: `test_*.py`
- Test function names: `test_description_of_what_is_tested`
- Use descriptive test names that explain the scenario
- Use arrange-act-assert pattern
- One assertion per test when possible

## Testing Guidelines

### Unit Tests

- Test individual functions in isolation
- Mock external dependencies (database, file I/O)
- Fast to run (milliseconds per test)
- No database required
- Mark with default pytest markers

### Integration Tests

- Test end-to-end functionality
- Require Neo4j database
- Slower to run (seconds per test)
- Mark with `@pytest.mark.database`

Example:
```python
import pytest

@pytest.mark.database
class TestPathwayValidation:
    """Integration tests requiring Neo4j."""

    def test_validates_against_database(self):
        # Test implementation
        pass
```

## Pull Request Process

1. **Ensure all tests pass**
   - Unit tests must pass
   - Integration tests should pass (if you can run them)

2. **Update documentation**
   - Update README.md if adding features
   - Add entry to CHANGELOG.md
   - Update docstrings

3. **Request review**
   - Tag relevant maintainers
   - Respond to feedback promptly
   - Make requested changes

4. **Merge requirements**
   - All CI checks must pass
   - At least one approval from maintainer
   - No merge conflicts with main branch

## Reporting Bugs

Use the [Bug Report](https://github.com/reactome/logic-network-generator/issues/new?template=bug_report.md) template and include:
- Clear description of the bug
- Steps to reproduce
- Expected vs actual behavior
- Environment details (OS, Python version, etc.)
- Error messages or logs

## Suggesting Features

Use the [Feature Request](https://github.com/reactome/logic-network-generator/issues/new?template=feature_request.md) template and include:
- Clear description of the feature
- Problem it solves
- Proposed solution
- Use cases and examples

## Questions?

- Open a [GitHub Discussion](https://github.com/reactome/logic-network-generator/discussions)
- Check existing issues and documentation
- Contact the maintainers

## Code of Conduct

- Be respectful and inclusive
- Welcome newcomers
- Focus on constructive feedback
- Assume good intentions

## License

By contributing, you agree that your contributions will be licensed under the Apache 2.0 License.
