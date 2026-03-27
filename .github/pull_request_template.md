## Description

Brief description of what this PR does.

## Type of Change

- [ ] Bug fix (non-breaking change that fixes an issue)
- [ ] New feature (non-breaking change that adds functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to not work as expected)
- [ ] Documentation update
- [ ] Code quality improvement (refactoring, performance, etc.)

## Related Issue

Fixes #(issue number)

## Changes Made

- Change 1
- Change 2
- Change 3

## Testing

### Unit Tests
- [ ] All existing unit tests pass locally (`poetry run pytest tests/ -v -m "not database"`)
- [ ] Added new unit tests for changes (if applicable)

### Integration Tests (Optional - requires Neo4j)
- [ ] All integration tests pass locally (`poetry run pytest tests/ -v`)

### Manual Testing
Describe any manual testing performed:
- Tested with pathway ID(s):
- Verified output files:

## Code Quality

- [ ] Code follows project style guidelines (ruff)
- [ ] Ran `poetry run ruff check src/` with no errors
- [ ] Ran `poetry run ruff format src/`
- [ ] Type hints added/updated where applicable
- [ ] Ran `poetry run mypy --ignore-missing-imports src/` (optional)

## Documentation

- [ ] Updated README.md (if needed)
- [ ] Updated CHANGELOG.md
- [ ] Added/updated docstrings
- [ ] Updated relevant documentation in `docs/`

## Checklist

- [ ] Self-review completed
- [ ] Code is well-commented, particularly in complex areas
- [ ] No debugging code left in (print statements, breakpoints, etc.)
- [ ] No credentials or sensitive information in code
- [ ] Git commit messages are clear and descriptive

## Screenshots (if applicable)

Add screenshots or terminal output if it helps explain the changes.

## Additional Notes

Any additional information that reviewers should know.
