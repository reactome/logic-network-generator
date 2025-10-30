"""
Example showing improved code structure with:
- Type hints
- Input validation
- Clear variable names
- Good docstrings
- Error handling
- No global state

Compare this to the current implementation to see the improvements.
"""

from typing import Dict, List, Any, Tuple
import pandas as pd
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)


@dataclass
class TransformationEdge:
    """Represents a single transformation edge in the network."""
    reactant_uuid: str  # Molecule consumed (input)
    product_uuid: str   # Molecule produced (output)
    logic_type: str     # 'and' or 'or'
    edge_category: str  # 'input' or 'output'
    regulation: str = 'pos'  # 'pos' or 'neg'


class LogicNetworkGenerator:
    """
    Generates logic networks from Reactome pathway data.

    This class transforms biological pathway data into directed graphs where:
    - Nodes are molecules (identified by UUIDs)
    - Edges are transformations within reactions (reactant → product)
    - AND/OR logic indicates whether multiple sources are alternatives

    Example:
        >>> from py2neo import Graph
        >>> graph = Graph("bolt://localhost:7687", auth=("neo4j", "test"))
        >>> generator = LogicNetworkGenerator(graph)
        >>> network = generator.generate(
        ...     decomposed_mapping=pd.read_csv('mapping.csv'),
        ...     reaction_connections=pd.read_csv('connections.csv'),
        ...     best_matches=pd.read_csv('matches.csv')
        ... )
    """

    def __init__(self, neo4j_graph):
        """
        Initialize the generator.

        Args:
            neo4j_graph: Connected py2neo Graph instance
        """
        self.graph = neo4j_graph
        self._molecule_uuid_cache: Dict[int, str] = {}

    def generate(
        self,
        decomposed_mapping: pd.DataFrame,
        reaction_connections: pd.DataFrame,
        best_matches: pd.DataFrame,
    ) -> pd.DataFrame:
        """
        Generate a logic network from pathway data.

        Args:
            decomposed_mapping: DataFrame with columns:
                - uid: Hash of molecule combination
                - reactome_id: Biological reaction ID
                - input_or_output_reactome_id: Terminal molecule ID
            reaction_connections: DataFrame with columns:
                - preceding_reaction_id: Upstream reaction
                - following_reaction_id: Downstream reaction
            best_matches: DataFrame with columns:
                - incomming: Input hash (within reaction)
                - outgoing: Output hash (within reaction)

        Returns:
            DataFrame representing the logic network with columns:
                - source_id: UUID of input molecule (reactant)
                - target_id: UUID of output molecule (product)
                - and_or: Logic type ('and' or 'or')
                - edge_type: Edge category ('input', 'output', etc.)
                - pos_neg: Regulation type ('pos' or 'neg')

        Raises:
            ValueError: If input DataFrames are invalid
            RuntimeError: If network generation fails
        """
        # Validate inputs
        self._validate_inputs(decomposed_mapping, reaction_connections, best_matches)

        try:
            # Create virtual reactions from best matches
            virtual_reactions = self._create_virtual_reactions(
                decomposed_mapping, best_matches
            )

            # Generate transformation edges
            edges = self._generate_transformation_edges(
                virtual_reactions, decomposed_mapping
            )

            # Add catalyst and regulator edges
            edges.extend(
                self._generate_catalyst_edges(virtual_reactions)
            )

            # Convert to DataFrame
            return self._edges_to_dataframe(edges)

        except Exception as e:
            logger.error(f"Failed to generate network: {e}")
            raise RuntimeError(f"Network generation failed: {e}") from e

    def _validate_inputs(
        self,
        decomposed_mapping: pd.DataFrame,
        reaction_connections: pd.DataFrame,
        best_matches: pd.DataFrame,
    ) -> None:
        """
        Validate input DataFrames have required structure.

        Raises:
            ValueError: If validation fails
        """
        # Check not empty
        if decomposed_mapping.empty:
            raise ValueError("decomposed_mapping cannot be empty")
        if best_matches.empty:
            raise ValueError("best_matches cannot be empty")

        # Check required columns
        required_mapping_cols = {'uid', 'reactome_id', 'input_or_output_reactome_id'}
        missing = required_mapping_cols - set(decomposed_mapping.columns)
        if missing:
            raise ValueError(
                f"decomposed_mapping missing columns: {missing}"
            )

        required_matches_cols = {'incomming', 'outgoing'}
        missing = required_matches_cols - set(best_matches.columns)
        if missing:
            raise ValueError(
                f"best_matches missing columns: {missing}"
            )

        logger.info("Input validation passed")

    def _generate_transformation_edges(
        self,
        virtual_reactions: List[Dict[str, Any]],
        decomposed_mapping: pd.DataFrame,
    ) -> List[TransformationEdge]:
        """
        Generate edges representing biochemical transformations.

        Each virtual reaction's inputs are connected to its outputs,
        representing the transformation that occurs.

        Args:
            virtual_reactions: List of reaction dictionaries
            decomposed_mapping: Mapping from hashes to molecules

        Returns:
            List of TransformationEdge objects
        """
        edges = []

        for reaction in virtual_reactions:
            # Extract terminal molecules
            reactant_ids = self._extract_terminal_molecules(
                decomposed_mapping, reaction['input_hash']
            )
            product_ids = self._extract_terminal_molecules(
                decomposed_mapping, reaction['output_hash']
            )

            # Skip if no terminal molecules
            if not reactant_ids or not product_ids:
                continue

            # Assign UUIDs to molecules
            reactant_uuids = [
                self._get_or_create_uuid(mol_id) for mol_id in reactant_ids
            ]
            product_uuids = [
                self._get_or_create_uuid(mol_id) for mol_id in product_ids
            ]

            # Determine AND/OR logic based on number of preceding reactions
            num_preceding = reaction['num_preceding_reactions']
            logic_type, edge_category = self._determine_logic(num_preceding)

            # Create cartesian product of reactants × products
            for reactant_uuid in reactant_uuids:
                for product_uuid in product_uuids:
                    edges.append(TransformationEdge(
                        reactant_uuid=reactant_uuid,
                        product_uuid=product_uuid,
                        logic_type=logic_type,
                        edge_category=edge_category,
                    ))

        logger.info(f"Generated {len(edges)} transformation edges")
        return edges

    def _determine_logic(self, num_preceding: int) -> Tuple[str, str]:
        """
        Determine AND/OR logic based on number of preceding reactions.

        Logic:
        - Single source (num_preceding == 1) → AND (required)
        - Multiple sources (num_preceding > 1) → OR (alternatives)

        Args:
            num_preceding: Number of reactions feeding into this one

        Returns:
            Tuple of (logic_type, edge_category)
        """
        if num_preceding > 1:
            return ('or', 'output')
        else:
            return ('and', 'input')

    def _extract_terminal_molecules(
        self,
        decomposed_mapping: pd.DataFrame,
        hash_value: str
    ) -> List[int]:
        """
        Extract terminal molecule IDs for a given hash.

        Terminal molecules are those that weren't further decomposed
        (e.g., individual proteins, not complexes).

        Args:
            decomposed_mapping: DataFrame containing mappings
            hash_value: Hash to look up

        Returns:
            List of Reactome IDs for terminal molecules
        """
        rows = decomposed_mapping[decomposed_mapping['uid'] == hash_value]
        terminal_ids = rows['input_or_output_reactome_id'].dropna().unique()
        return [int(id) for id in terminal_ids]

    def _get_or_create_uuid(self, reactome_id: int) -> str:
        """
        Get or create a UUID for a Reactome ID.

        Uses caching to ensure the same Reactome ID always gets
        the same UUID.

        Args:
            reactome_id: Reactome database ID

        Returns:
            UUID string for this molecule
        """
        if reactome_id not in self._molecule_uuid_cache:
            import uuid
            self._molecule_uuid_cache[reactome_id] = str(uuid.uuid4())

        return self._molecule_uuid_cache[reactome_id]

    def _create_virtual_reactions(
        self,
        decomposed_mapping: pd.DataFrame,
        best_matches: pd.DataFrame,
    ) -> List[Dict[str, Any]]:
        """
        Create virtual reactions from best matches.

        Each best match represents a pairing of input/output molecule
        combinations that forms a virtual reaction.

        Args:
            decomposed_mapping: Mapping from hashes to reactions
            best_matches: Pairings of input and output hashes

        Returns:
            List of virtual reaction dictionaries
        """
        virtual_reactions = []

        for _, match in best_matches.iterrows():
            incoming_hash = match['incomming']
            outgoing_hash = match['outgoing']

            # Get the biological reaction ID
            reactome_id = self._get_reactome_id_from_hash(
                decomposed_mapping, incoming_hash
            )

            virtual_reactions.append({
                'reactome_id': reactome_id,
                'input_hash': incoming_hash,
                'output_hash': outgoing_hash,
                'num_preceding_reactions': 1,  # Simplified for example
            })

        return virtual_reactions

    def _get_reactome_id_from_hash(
        self,
        decomposed_mapping: pd.DataFrame,
        hash_value: str
    ) -> int:
        """
        Extract Reactome ID for a given hash.

        Args:
            decomposed_mapping: Mapping DataFrame
            hash_value: Hash to look up

        Returns:
            Reactome ID as integer

        Raises:
            ValueError: If hash not found
        """
        result = decomposed_mapping.loc[
            decomposed_mapping['uid'] == hash_value, 'reactome_id'
        ].values

        if len(result) == 0:
            raise ValueError(f"Hash not found: {hash_value}")

        return int(result[0])

    def _generate_catalyst_edges(
        self,
        virtual_reactions: List[Dict[str, Any]]
    ) -> List[TransformationEdge]:
        """
        Generate edges for catalysts.

        (Simplified placeholder - real implementation would query Neo4j)
        """
        # TODO: Implement catalyst edge generation
        return []

    def _edges_to_dataframe(
        self,
        edges: List[TransformationEdge]
    ) -> pd.DataFrame:
        """
        Convert TransformationEdge objects to DataFrame.

        Args:
            edges: List of edge objects

        Returns:
            DataFrame with standard column names
        """
        return pd.DataFrame([
            {
                'source_id': edge.reactant_uuid,
                'target_id': edge.product_uuid,
                'and_or': edge.logic_type,
                'edge_type': edge.edge_category,
                'pos_neg': edge.regulation,
            }
            for edge in edges
        ])


# Example usage
if __name__ == '__main__':
    # This is a usage example - requires actual data files
    print("""
    Example usage:

    from py2neo import Graph

    # Connect to database
    graph = Graph("bolt://localhost:7687", auth=("neo4j", "test"))

    # Create generator
    generator = LogicNetworkGenerator(graph)

    # Load data
    mapping = pd.read_csv('decomposed_uid_mapping_69620.csv')
    connections = pd.read_csv('reaction_connections_69620.csv')
    matches = pd.read_csv('best_matches_69620.csv')

    # Generate network
    network = generator.generate(mapping, connections, matches)

    # Save result
    network.to_csv('pathway_logic_network_69620.csv', index=False)
    print(f"Generated network with {len(network)} edges")
    """)
