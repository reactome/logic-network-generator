# Examples

This directory contains example scripts demonstrating how to use the Logic Network Generator.

## Available Examples

### 1. `generate_pathway_example.py`

**Purpose**: Complete example showing how to generate and analyze a pathway logic network.

**What it demonstrates**:
- Generating a logic network for a specific Reactome pathway
- Analyzing network properties (edges, nodes, logic relationships)
- Finding root inputs and terminal outputs
- Handling common errors (connection failures, invalid pathways)

**Usage**:
```bash
# Ensure Neo4j is running
docker run -p 7474:7474 -p 7687:7687 \
  -e NEO4J_dbms_memory_heap_maxSize=8g \
  public.ecr.aws/reactome/graphdb:Release94

# Run the example
poetry run python examples/generate_pathway_example.py
```

**Expected Output**:
```
Logic Network Generator - Example Usage
======================================================================

Generating logic network for pathway: Cell Cycle, Mitotic
Pathway ID: 69620

Step 1: Fetching reactions from Neo4j...
Step 2: Decomposing complexes and entity sets...
Step 3: Creating logic network...

======================================================================
Generation Complete!
======================================================================

Network Analysis:
  Total edges: 4995

  Edge types:
    - input: 3200
    - output: 1200
    - catalyst: 350
    - regulator: 245

  Logic relationships:
    - AND edges (required): 4100
    - OR edges (alternatives): 895

  Network structure:
    - Root inputs (starting points): 9
    - Terminal outputs (endpoints): 11
    - Unique physical entities: 458
```

## Example Pathways

Here are some good pathways to try:

| Pathway ID | Pathway Name | Complexity | Description |
|------------|-------------|------------|-------------|
| 69620 | Cell Cycle, Mitotic | Medium | Well-studied cell cycle pathway |
| 68875 | Apoptosis | Medium | Programmed cell death pathway |
| 1640170 | Cell Cycle | Large | Complete cell cycle regulation |
| 112316 | Neuronal System | Large | Neural signaling pathways |
| 382551 | Transport of small molecules | Large | Molecular transport mechanisms |

## Common Usage Patterns

### Pattern 1: Generate Multiple Pathways

```python
pathway_ids = ["69620", "68875", "112316"]

for pathway_id in pathway_ids:
    generate_pathway_file(
        pathway_id=pathway_id,
        taxon_id="9606",
        pathway_name=f"Pathway_{pathway_id}",
        decompose=False
    )
```

### Pattern 2: Load and Analyze Existing Network

```python
import pandas as pd
from src.logic_network_generator import find_root_inputs, find_terminal_outputs

# Load previously generated network
network = pd.read_csv("pathway_logic_network_69620.csv")

# Find starting and ending points
roots = find_root_inputs(network)
terminals = find_terminal_outputs(network)

# Analyze specific subsets
and_edges = network[network['and_or'] == 'and']
or_edges = network[network['and_or'] == 'or']

print(f"Network has {len(roots)} entry points and {len(terminals)} exit points")
print(f"AND edges: {len(and_edges)}, OR edges: {len(or_edges)}")
```

### Pattern 3: Export for Cytoscape

```python
import pandas as pd

# Load network
network = pd.read_csv("pathway_logic_network_69620.csv")

# Create Cytoscape-compatible format
cytoscape_edges = network[['source_id', 'target_id', 'and_or', 'edge_type']].copy()
cytoscape_edges.columns = ['Source', 'Target', 'Logic', 'EdgeType']

# Save for Cytoscape import
cytoscape_edges.to_csv("network_for_cytoscape.csv", index=False)
print("Exported to network_for_cytoscape.csv")
print("Import in Cytoscape: File → Import → Network from File")
```

## Troubleshooting

### Neo4j Connection Issues

**Error**: `ConnectionError: Failed to connect to Neo4j database`

**Solution**:
```bash
# Check if Neo4j is running
docker ps | grep reactome

# Start Neo4j if not running
docker run -p 7474:7474 -p 7687:7687 \
  -e NEO4J_dbms_memory_heap_maxSize=8g \
  public.ecr.aws/reactome/graphdb:Release94

# Wait 30 seconds for Neo4j to start, then try again
```

### Invalid Pathway ID

**Error**: `ValueError: No reactions found for pathway ID: 12345`

**Solution**:
- Verify the pathway ID exists at https://reactome.org/PathwayBrowser/
- Check that you're using the numeric database ID (not the stable identifier)
- Try a known working pathway like 69620

### Out of Memory

**Error**: `MemoryError` or very slow performance

**Solution**:
- Start with smaller pathways (< 500 reactions)
- Increase Neo4j memory: `-e NEO4J_dbms_memory_heap_maxSize=16g`
- Run on a machine with more RAM

## Additional Resources

- **Architecture Documentation**: `docs/ARCHITECTURE.md`
- **Test Suite**: `tests/` directory with 43 tests
- **Improvement Ideas**: `IMPROVEMENT_RECOMMENDATIONS.md`
- **Reactome Database**: https://reactome.org/
