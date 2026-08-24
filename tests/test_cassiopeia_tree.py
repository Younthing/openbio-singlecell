from __future__ import annotations

import copy
from dataclasses import FrozenInstanceError
from types import SimpleNamespace

import pytest

from openbio_singlecell import cassiopeia_tree
from openbio_singlecell.cassiopeia_tree import CassiopeiaTree
from openbio_singlecell.contracts import TableResult
from openbio_singlecell.node_types import AnnDataType, CassiopeiaTreeType, TableResultType
from openbio_singlecell.nodes_lineage import (
    LINEAGE_NODE_CLASSES,
    OpenBioSingleCellCassiopeiaExpansionTest,
    OpenBioSingleCellCassiopeiaPlasticity,
    OpenBioSingleCellReconstructCassiopeiaTree,
)


def output_values(node_output):
    return node_output.result


class FakeSolvedTree:
    def __init__(self, character_matrix, priors):
        self.n_cell = int(character_matrix.shape[0])
        self.root = "root"
        self.nodes = ["root", "clade", "cell_1", "cell_2"]
        self.leaves = ["cell_1", "cell_2"]
        self.character_matrix_received = character_matrix
        self.priors_received = priors
        self.attributes = {}
        self.cell_meta = None
        self.solved = False

    def copy(self):
        return copy.deepcopy(self)

    def depth_first_traverse_nodes(self):
        return list(self.nodes)

    def get_attribute(self, node, name):
        return self.attributes[(node, name)]

    def leaves_in_subtree(self, node):
        if node in self.leaves:
            return [node]
        return list(self.leaves)

    def parent(self, node):
        return {"clade": "root", "cell_1": "clade", "cell_2": "clade"}[node]

    def set_attribute(self, node, name, value):
        self.attributes[(node, name)] = value


@pytest.fixture
def fake_cassiopeia(science, monkeypatch):
    calls = {
        "construct": 0,
        "convert": 0,
        "expansion": 0,
        "priors": 0,
        "read": 0,
        "solve": 0,
    }
    allele_table = science.pd.DataFrame(
        {
            "sample": ["tumor-a", "tumor-a", "tumor-b"],
            "barcode": ["cell_1", "cell_2", "other_cell"],
            "integration": ["int_1", "int_1", "int_2"],
            "family": ["mut_1", "mut_2", "mut_3"],
        }
    )

    def read_allele_table(path, first_column_as_index):
        calls["read"] += 1
        calls["read_arguments"] = (path, first_column_as_index)
        return allele_table.copy(deep=True)

    def compute_empirical_indel_priors(table, grouping_variables):
        calls["priors"] += 1
        calls["prior_columns"] = tuple(table.columns)
        calls["prior_rows"] = int(table.shape[0])
        calls["grouping_variables"] = tuple(grouping_variables)
        return {"indel": 0.5}

    def convert_alleletable_to_character_matrix(table, *, allele_rep_thresh, mutation_priors):
        calls["convert"] += 1
        calls["converted_tumors"] = tuple(table["Tumor"])
        calls["allele_rep_thresh"] = allele_rep_thresh
        calls["mutation_priors"] = mutation_priors
        matrix = science.pd.DataFrame(
            [[0, 1], [1, 0]],
            index=["cell_1", "cell_2"],
            columns=["character_1", "character_2"],
        )
        return matrix, {"character_1": {1: 0.5}}, None

    def construct_tree(*, character_matrix, priors):
        calls["construct"] += 1
        tree = FakeSolvedTree(character_matrix, priors)
        calls["constructed_tree"] = tree
        return tree

    class FakeSolver:
        def solve(self, tree):
            calls["solve"] += 1
            tree.solved = True

    def compute_expansion_pvalues(tree, *, min_clade_size, min_depth):
        calls["expansion"] += 1
        calls["expansion_arguments"] = (min_clade_size, min_depth)
        for node, pvalue in {
            "root": 0.5,
            "clade": 0.005,
            "cell_1": 0.8,
            "cell_2": 0.8,
        }.items():
            tree.set_attribute(node, "expansion_pvalue", pvalue)

    def score_small_parsimony(tree, *, meta_item, root=None):
        calls.setdefault("parsimony", []).append((meta_item, root))
        if root is None:
            return 4
        return len(tree.leaves_in_subtree(root))

    fake = SimpleNamespace(
        data=SimpleNamespace(CassiopeiaTree=construct_tree),
        pp=SimpleNamespace(
            compute_empirical_indel_priors=compute_empirical_indel_priors,
            convert_alleletable_to_character_matrix=convert_alleletable_to_character_matrix,
        ),
        solver=SimpleNamespace(VanillaGreedySolver=FakeSolver),
        tl=SimpleNamespace(
            compute_expansion_pvalues=compute_expansion_pvalues,
            score_small_parsimony=score_small_parsimony,
        ),
    )
    monkeypatch.setattr(cassiopeia_tree, "_require_cassiopeia", lambda: (fake, object()))
    monkeypatch.setattr(cassiopeia_tree, "_read_allele_table", read_allele_table)
    return calls


@pytest.fixture
def reconstructed_tree(fake_cassiopeia):
    (tree,) = output_values(
        OpenBioSingleCellReconstructCassiopeiaTree.execute(
            "alleles.tsv",
            tumor=" tumor-a ",
            first_column_as_index=False,
            tumor_column="sample",
            cell_barcode_column="barcode",
            integration_barcode_column="integration",
            mutation_family_column="family",
            allele_representation_threshold=0.75,
        )
    )
    return tree


def test_lineage_schemas_use_one_reusable_tree_contract():
    reconstruct = OpenBioSingleCellReconstructCassiopeiaTree.GET_SCHEMA()
    expansion = OpenBioSingleCellCassiopeiaExpansionTest.GET_SCHEMA()
    plasticity = OpenBioSingleCellCassiopeiaPlasticity.GET_SCHEMA()

    assert [input_.id for input_ in reconstruct.inputs] == [
        "allele_table_file",
        "tumor",
        "first_column_as_index",
        "tumor_column",
        "cell_barcode_column",
        "integration_barcode_column",
        "mutation_family_column",
        "allele_representation_threshold",
    ]
    assert [(output.display_name, output.io_type) for output in reconstruct.outputs] == [
        ("tree", CassiopeiaTreeType.io_type)
    ]
    assert [input_.id for input_ in expansion.inputs] == [
        "tree",
        "minimum_clade_fraction",
        "minimum_depth",
        "expansion_pvalue_threshold",
    ]
    assert expansion.inputs[0].io_type == CassiopeiaTreeType.io_type
    assert expansion.outputs[0].io_type == TableResultType.io_type
    assert [input_.id for input_ in plasticity.inputs] == [
        "adata",
        "tree",
        "annotation_key",
        "output_key",
        "summary_key",
    ]
    assert plasticity.inputs[0].io_type == AnnDataType.io_type
    assert plasticity.inputs[1].io_type == CassiopeiaTreeType.io_type
    assert "validate_inputs" in OpenBioSingleCellReconstructCassiopeiaTree.__dict__
    assert "fingerprint_inputs" in OpenBioSingleCellReconstructCassiopeiaTree.__dict__
    assert "validate_inputs" not in OpenBioSingleCellCassiopeiaExpansionTest.__dict__
    assert "fingerprint_inputs" not in OpenBioSingleCellCassiopeiaExpansionTest.__dict__
    assert "validate_inputs" not in OpenBioSingleCellCassiopeiaPlasticity.__dict__
    assert "fingerprint_inputs" not in OpenBioSingleCellCassiopeiaPlasticity.__dict__
    assert LINEAGE_NODE_CLASSES.index(OpenBioSingleCellReconstructCassiopeiaTree) < LINEAGE_NODE_CLASSES.index(
        OpenBioSingleCellCassiopeiaExpansionTest
    )
    assert LINEAGE_NODE_CLASSES.index(OpenBioSingleCellReconstructCassiopeiaTree) < LINEAGE_NODE_CLASSES.index(
        OpenBioSingleCellCassiopeiaPlasticity
    )


def test_reconstruction_builds_and_solves_the_tree_once(reconstructed_tree, fake_cassiopeia):
    assert isinstance(reconstructed_tree, CassiopeiaTree)
    assert reconstructed_tree.tumor == "tumor-a"
    assert reconstructed_tree.input_cells == 2
    assert reconstructed_tree.character_count == 2
    assert reconstructed_tree.provenance["operation"] == "reconstruct_cassiopeia_tree"
    assert reconstructed_tree.provenance["parameters"]["allele_representation_threshold"] == 0.75
    assert not hasattr(reconstructed_tree, "solver")
    assert not hasattr(reconstructed_tree, "character_matrix")
    assert not hasattr(reconstructed_tree, "priors")

    assert fake_cassiopeia["read"] == 1
    assert fake_cassiopeia["read_arguments"] == ("alleles.tsv", False)
    assert fake_cassiopeia["priors"] == 1
    assert fake_cassiopeia["prior_rows"] == 3
    assert fake_cassiopeia["prior_columns"] == ("Tumor", "cellBC", "intBC", "MetFamily")
    assert fake_cassiopeia["grouping_variables"] == ("intBC", "MetFamily")
    assert fake_cassiopeia["convert"] == 1
    assert fake_cassiopeia["converted_tumors"] == ("tumor-a", "tumor-a")
    assert fake_cassiopeia["mutation_priors"] == {"indel": 0.5}
    assert fake_cassiopeia["construct"] == 1
    assert fake_cassiopeia["solve"] == 1
    assert fake_cassiopeia["constructed_tree"].solved is True

    with pytest.raises(FrozenInstanceError):
        reconstructed_tree.tumor = "other"
    with pytest.raises(TypeError):
        reconstructed_tree.provenance["operation"] = "changed"


def test_consumers_reuse_solved_tree_without_reconstruction(reconstructed_tree, fake_cassiopeia, science):
    (expansion_result,) = output_values(
        OpenBioSingleCellCassiopeiaExpansionTest.execute(
            reconstructed_tree,
            minimum_clade_fraction=0.25,
            minimum_depth=1,
            expansion_pvalue_threshold=0.01,
        )
    )
    assert isinstance(expansion_result, TableResult)
    assert expansion_result.table.loc[expansion_result.table["node"] == "clade", "is_expansion"].item() is True
    assert expansion_result.parameters["tree_cells"] == 2
    assert expansion_result.parameters["tree_characters"] == 2
    assert fake_cassiopeia["expansion_arguments"] == (0.5, 1)

    obs = science.pd.DataFrame({"cell_type": ["A", "B"]}, index=["cell_1", "cell_2"])
    var = science.pd.DataFrame(index=["gene_1"])
    adata = science.ad.AnnData(science.np.asarray([[1.0], [2.0]]), obs=obs, var=var)
    (output,) = output_values(
        OpenBioSingleCellCassiopeiaPlasticity.execute(
            adata,
            reconstructed_tree,
            annotation_key="cell_type",
            output_key="plasticity",
            summary_key="lineage_plasticity",
        )
    )

    assert "plasticity" not in adata.obs
    assert "lineage_plasticity" not in adata.uns
    assert output.obs["plasticity"].tolist() == [1.0, 1.0]
    assert output.uns["lineage_plasticity"] == {
        "tumor": "tumor-a",
        "annotation_key": "cell_type",
        "parsimony": 4,
        "effective_plasticity_score": 1.0,
        "tree_nodes": 4,
        "tree_leaves": 2,
        "tree_characters": 2,
    }
    history = output.uns["openbio_singlecell"]["analysis_history"]
    assert history["000000"]["operation"] == "cassiopeia_plasticity"

    assert fake_cassiopeia["read"] == 1
    assert fake_cassiopeia["priors"] == 1
    assert fake_cassiopeia["convert"] == 1
    assert fake_cassiopeia["construct"] == 1
    assert fake_cassiopeia["solve"] == 1
    assert fake_cassiopeia["constructed_tree"].attributes == {}
    assert fake_cassiopeia["constructed_tree"].cell_meta is None


def test_reconstruction_rejects_invalid_threshold_before_file_access(fake_cassiopeia):
    with pytest.raises(ValueError, match="allele_representation_threshold must be between 0 and 1"):
        OpenBioSingleCellReconstructCassiopeiaTree.execute(
            "alleles.tsv",
            tumor="tumor-a",
            allele_representation_threshold=1.1,
        )

    assert fake_cassiopeia["read"] == 0


def test_tree_wrapper_rejects_values_without_solved_tree_capabilities():
    with pytest.raises(TypeError, match="missing required solved-tree capabilities"):
        CassiopeiaTree(
            tumor="tumor-a",
            input_cells=2,
            character_count=1,
            provenance={"operation": "test"},
            _solved_tree=object(),
        )

    with pytest.raises(ValueError, match="provenance cannot be empty"):
        CassiopeiaTree(
            tumor="tumor-a",
            input_cells=2,
            character_count=1,
            provenance={},
            _solved_tree=object(),
        )
