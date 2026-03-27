#!/usr/bin/env poetry run python
"""Run comprehensive pathway validation.

Usage:
    poetry run python validate_pathway.py [pathway_id]

Example:
    poetry run python validate_pathway.py 69620
"""

import sys
import subprocess
from pathlib import Path

def main():
    # Get pathway ID from command line or use default
    pathway_id = sys.argv[1] if len(sys.argv) > 1 else "69620"

    print(f"Running comprehensive validation for pathway {pathway_id}...")
    print("=" * 80)

    # Run the validation tests
    result = subprocess.run(
        ["poetry", "run", "pytest", "tests/test_pathway_validation.py", "-v", "-s"],
        cwd=Path(__file__).parent
    )

    sys.exit(result.returncode)

if __name__ == "__main__":
    main()
