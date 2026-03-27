# Changelog

All notable changes to this project.

## [0.2.0] - 2025-11-11

### Added
- **Position-Aware UUIDs**: Same entity at different pathway positions now receives unique UUIDs, eliminating unwanted self-loops
- **UUID Mapping Export**: Maps UUIDs back to Reactome IDs with position context (`uuid_mapping_{pathway_id}.csv`)
- **Comprehensive Validation System**: 11 tests validate logic networks against source database
  - Loop/cycle analysis
  - Regulator matching
  - Identifier resolution (UniProt, gene symbols, Ensembl)
  - Root input identification
  - Topological equivalence
  - Information loss checking
- **Ultra-Comprehensive Validation**: 8 additional tests for production confidence
  - Find root inputs by UniProt (e.g., TP53)
  - Trace entities through all positions
  - Verify no spurious loops introduced
- **Output Folder Organization**: All generated files now saved to `output/` directory

### Fixed
- Self-loop bug where same entity at different positions incorrectly merged into single node
- Test portability - removed hardcoded local paths

### Changed
- Output files relocated from root to `output/` folder for better organization
- Test suite expanded from 52 to 73+ tests (including position-aware UUID tests)
- Enhanced logging for UUID registry statistics and union-find operations

## [0.1.0] - 2025-01-29

### Added
- **Database ID Mapping Tool**: Convert Reactome IDs to human-readable names with full CLI options
- **Regulator Tests**: 9 comprehensive tests for negative regulators, positive regulators, and catalysts
- **Usage Examples**: Working examples in `examples/` directory with documentation
- **Architecture Documentation**: Complete system architecture and design decisions in `docs/ARCHITECTURE.md`
- **Error Handling**: Comprehensive error messages with troubleshooting guidance
- **Type Hints**: Added type annotations across codebase (~95% coverage)
- **Input Validation**: Validate DataFrame inputs with helpful error messages
- **CI/CD**: GitHub Actions workflow for automated testing
- **Coverage Reporting**: pytest-cov integration with HTML reports

### Changed
- Terminology alignment: "molecule" → "physical entity" to match Reactome schema
- Enhanced logging throughout codebase
- Improved function documentation with detailed docstrings

### Removed
- Debug print statements and verbose logging
- Temporary instrumentation code

### Testing
- Test suite: 52 tests with 100% pass rate
- Coverage configuration in `pyproject.toml`
- Pytest configuration for consistent test execution

## Initial Release

### Core Features
- Generate logic networks from Reactome pathways
- Decompose complexes and entity sets into components
- AND/OR logic determination based on pathway structure
- Support for negative regulators, positive regulators, and catalysts
- Neo4j database integration
- Batch processing with pathway lists
- Caching for improved performance
