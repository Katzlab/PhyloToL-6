# Unit tests for Scripts/contamination.py

# Leaf names follow the EukPhylo convention: the first 10 characters are the taxon code
# (e.g. Sr_ci_AAAA), and the first 2 are the major clade (Ba, Za, Op, Pl, Am, Ex, Sr).
# Trees are kept tiny so the expected result of each test can be worked out by hand.

from types import SimpleNamespace
from unittest import mock

import ete3
import pytest

import contamination


def make_args(**overrides):
	defaults = dict(
		clade_grabbing_exceptions=None,
		clade_grabbing_rules_file=None,
		target_taxa='Sr_ci',
		num_contams=0,
		min_target_presence=2,
		required_taxa=None,
		required_taxa_num=0,
		cocontaminants=None,
		subsister_rules=None,
		cl_exclude_taxa=None,
		output=None,
	)
	defaults.update(overrides)
	return SimpleNamespace(**defaults)


def write_file(tmp_path, name, content):
	path = tmp_path / name
	path.write_text(content)
	return str(path)


def leaf_names(tree):
	return sorted(leaf.name for leaf in tree)


# --------------------------------------------------------------------------------------
# test get_newick
# --------------------------------------------------------------------------------------

def test_get_newick_plain_newick(tmp_path):
	f = write_file(tmp_path, 'OG1.tre', '((A:1,B:1):1,C:1);')
	assert contamination.get_newick(f) == '((A:1,B:1):1,C:1);'


def test_get_newick_nexus_tree1_line(tmp_path):
	nexus = '#NEXUS\nbegin trees;\ntree tree1=((A:1,B:1):1,C:1);\nend;\n'
	f = write_file(tmp_path, 'OG1.nex', nexus)
	assert contamination.get_newick(f).strip() == '((A:1,B:1):1,C:1);'


def test_get_newick_strips_quotes_and_backslashes(tmp_path):
	f = write_file(tmp_path, 'OG1.tre', "(('A':1,'B\\':1):1,C:1);")
	assert contamination.get_newick(f) == '((A:1,B:1):1,C:1);'


def test_get_newick_no_tree_returns_empty_string(tmp_path):
	f = write_file(tmp_path, 'OG1.tre', 'not a tree\n')
	assert contamination.get_newick(f) == ''


# --------------------------------------------------------------------------------------
# test reroot
# --------------------------------------------------------------------------------------

def root_child_leaf_sets(tree):
	return [set(leaf.name for leaf in child) for child in tree.children]


def test_reroot_on_largest_prokaryote_clade():
	# The sister group (Am, Pl, Ex) has 3 eukaryotes, so the prokaryote clade can't absorb it
	# (reroot tolerates up to 2 non-target leaves in the chosen clade)
	prok = {'Ba_bb_Bsub_1', 'Ba_bb_Ecol_1', 'Za_aa_Mjan_1', 'Ba_cc_Paer_1'}
	tree = ete3.Tree('(Op_me_Hsap_1,(((Ba_bb_Bsub_1,Ba_bb_Ecol_1),(Za_aa_Mjan_1,Ba_cc_Paer_1)),(Am_tu_Xxxx_1,Pl_gr_Atha_1,Ex_eu_Yyyy_1)));')

	tree = contamination.reroot(tree)

	assert prok in root_child_leaf_sets(tree)


def test_reroot_falls_back_to_opisthokonts_without_prokaryotes():
	op = {'Op_me_Hsap_1', 'Op_me_Mmus_1', 'Op_fu_Scer_1', 'Op_fu_Ncra_1'}
	tree = ete3.Tree('((((Op_me_Hsap_1,Op_me_Mmus_1),(Op_fu_Scer_1,Op_fu_Ncra_1)),(Pl_gr_Atha_1,Pl_gr_Osat_1,Pl_gr_Zmay_1,Pl_gr_Ppat_1)),(Am_tu_Xxxx_1,Am_tu_Yyyy_1));')

	tree = contamination.reroot(tree)

	assert op in root_child_leaf_sets(tree)


def test_reroot_leaves_tree_unchanged_without_large_enough_clade():
	newick = '((Op_me_Hsap_1,Op_me_Mmus_1,Op_fu_Scer_1),(Am_tu_Xxxx_1,Pl_gr_Atha_1));'
	tree = ete3.Tree(newick)

	tree = contamination.reroot(tree)

	assert tree.write(format=9) == ete3.Tree(newick).write(format=9)


# --------------------------------------------------------------------------------------
# test get_subtrees (clade-based contamination removal)
# --------------------------------------------------------------------------------------

# Three ciliates (A, B, C) form a clean clade; a fourth (D) sits inside a clade of other taxa
TREE_ONE_NESTED_CILIATE = '(((Sr_ci_AAAA_1,Sr_ci_BBBB_1),Sr_ci_CCCC_1),((Op_me_Hsap_1,Sr_ci_DDDD_1),(Am_tu_Xxxx_1,Pl_gr_Atha_1)));'
NON_TARGETS = ['Am_tu_Xxxx_1', 'Op_me_Hsap_1', 'Pl_gr_Atha_1']

# Ciliates A, B and C form a clade together with one contaminant (Op)
TREE_ONE_CONTAMINANT = '(((Sr_ci_AAAA_1,Sr_ci_BBBB_1),(Sr_ci_CCCC_1,Op_me_Hsap_1)),(Am_tu_Xxxx_1,Pl_gr_Atha_1));'


def test_get_subtrees_removes_target_outside_robust_clade(tmp_path):
	f = write_file(tmp_path, 'OG1.tre', TREE_ONE_NESTED_CILIATE)

	kept = contamination.get_subtrees(make_args(), f)

	assert sorted(kept) == sorted(['Sr_ci_AAAA_1', 'Sr_ci_BBBB_1', 'Sr_ci_CCCC_1'] + NON_TARGETS)


def test_get_subtrees_always_keeps_non_target_sequences(tmp_path):
	f = write_file(tmp_path, 'OG1.tre', TREE_ONE_NESTED_CILIATE)

	kept = contamination.get_subtrees(make_args(), f)

	assert set(NON_TARGETS) <= set(kept)


@pytest.mark.parametrize('num_contams, ciliate_c_kept', [
	(0, False),    # no contaminants allowed: the A/B/C/Op clade is rejected
	(1, True),     # one contaminant allowed (absolute count)
	(0.5, True),   # 1 contaminant <= 0.5 * 3 targets (proportion)
	(0.3, False),  # 1 contaminant >  0.3 * 3 targets
])
def test_get_subtrees_num_contams_threshold(tmp_path, num_contams, ciliate_c_kept):
	f = write_file(tmp_path, 'OG1.tre', TREE_ONE_CONTAMINANT)

	kept = contamination.get_subtrees(make_args(num_contams=num_contams), f)

	assert {'Sr_ci_AAAA_1', 'Sr_ci_BBBB_1', 'Op_me_Hsap_1'} <= set(kept)
	assert ('Sr_ci_CCCC_1' in kept) == ciliate_c_kept


def test_get_subtrees_min_target_presence_not_met_removes_all_targets(tmp_path):
	f = write_file(tmp_path, 'OG1.tre', TREE_ONE_NESTED_CILIATE)

	kept = contamination.get_subtrees(make_args(min_target_presence=4), f)

	assert sorted(kept) == NON_TARGETS


def test_get_subtrees_target_taxa_from_file(tmp_path):
	f = write_file(tmp_path, 'OG1.tre', TREE_ONE_NESTED_CILIATE)
	targets = write_file(tmp_path, 'targets.txt', 'Sr_ci\n\n')

	kept = contamination.get_subtrees(make_args(target_taxa=targets), f)

	assert 'Sr_ci_DDDD_1' not in kept
	assert {'Sr_ci_AAAA_1', 'Sr_ci_BBBB_1', 'Sr_ci_CCCC_1'} <= set(kept)


def test_get_subtrees_required_taxa_present(tmp_path):
	f = write_file(tmp_path, 'OG1.tre', TREE_ONE_NESTED_CILIATE)

	kept = contamination.get_subtrees(make_args(required_taxa='Sr_ci_AAAA', required_taxa_num=1), f)

	assert {'Sr_ci_AAAA_1', 'Sr_ci_BBBB_1', 'Sr_ci_CCCC_1'} <= set(kept)
	assert 'Sr_ci_DDDD_1' not in kept


def test_get_subtrees_required_taxa_missing_removes_all_targets(tmp_path):
	f = write_file(tmp_path, 'OG1.tre', TREE_ONE_NESTED_CILIATE)

	kept = contamination.get_subtrees(make_args(required_taxa='Sr_ci_ZZZZ', required_taxa_num=1), f)

	assert sorted(kept) == NON_TARGETS


def test_get_subtrees_exceptions_are_always_kept(tmp_path):
	f = write_file(tmp_path, 'OG1.tre', TREE_ONE_NESTED_CILIATE)
	exceptions = write_file(tmp_path, 'exceptions.txt', 'Sr_ci_DDDD_1\n')

	kept = contamination.get_subtrees(make_args(clade_grabbing_exceptions=exceptions), f)

	assert 'Sr_ci_DDDD_1' in kept


def test_get_subtrees_rules_file(tmp_path):
	f = write_file(tmp_path, 'OG1.tre', TREE_ONE_NESTED_CILIATE)
	# Columns: target_taxa, num_contams, min_target_presence, required_taxa, required_taxa_num.
	# Lines without exactly 5 columns are ignored.
	rules = write_file(tmp_path, 'rules.tsv', 'Sr_ci\t0\t2\tna\tna\nbad\tline\n')

	kept = contamination.get_subtrees(make_args(clade_grabbing_rules_file=rules, target_taxa=None), f)

	assert sorted(kept) == sorted(['Sr_ci_AAAA_1', 'Sr_ci_BBBB_1', 'Sr_ci_CCCC_1'] + NON_TARGETS)


@pytest.mark.parametrize('arg', ['clade_grabbing_exceptions', 'clade_grabbing_rules_file'])
def test_get_subtrees_missing_input_file_exits(tmp_path, arg):
	f = write_file(tmp_path, 'OG1.tre', TREE_ONE_NESTED_CILIATE)

	with pytest.raises(SystemExit):
		contamination.get_subtrees(make_args(**{arg: str(tmp_path / 'does_not_exist.txt')}), f)


# --------------------------------------------------------------------------------------
# test get_sisters (sister-based contamination removal)
# --------------------------------------------------------------------------------------

# Ciliate A is sister to a single bacterium; all branch lengths 0.1
TREE_CILIATE_SISTER_TO_BACTERIUM = '((Sr_ci_AAAA_1:0.1,Ba_bb_Ecol_1:0.1):0.1,(Op_me_Hsap_1:0.1,Am_tu_Xxxx_1:0.1):0.1,(Pl_gr_Atha_1:0.1,Ex_eu_Yyyy_1:0.1):0.1);'

# Rule: remove Sr_ci sequences that are sister to Ba, if branch length <= 2 * mean branch length
SISTER_RULES = {'Sr_ci': [('Ba', 2.0)]}


def test_get_sisters_removes_sequence_sister_to_contaminant(tmp_path):
	f = write_file(tmp_path, 'OG1.tre', TREE_CILIATE_SISTER_TO_BACTERIUM)

	kept = contamination.get_sisters(make_args(), f, SISTER_RULES, {})

	assert 'Sr_ci_AAAA_1' not in kept
	assert len(kept) == 5


def test_get_sisters_keeps_sequence_with_long_branch(tmp_path):
	# Mean leaf branch length = (1.0 + 5 * 0.1) / 6 = 0.25, so the cutoff is 0.5 and A (1.0) is kept
	newick = TREE_CILIATE_SISTER_TO_BACTERIUM.replace('Sr_ci_AAAA_1:0.1', 'Sr_ci_AAAA_1:1.0')
	f = write_file(tmp_path, 'OG1.tre', newick)

	kept = contamination.get_sisters(make_args(), f, SISTER_RULES, {})

	assert 'Sr_ci_AAAA_1' in kept


def test_get_sisters_keeps_sequence_with_mixed_sisters(tmp_path):
	# A is sister to (Ba, Op); only Ba is a listed contaminant, so A is kept
	newick = '((Sr_ci_AAAA_1:0.1,(Ba_bb_Ecol_1:0.1,Op_me_Hsap_1:0.1):0.1):0.1,(Am_tu_Xxxx_1:0.1,Pl_gr_Atha_1:0.1):0.1,(Ex_eu_Yyyy_1:0.1,Ex_eu_Zzzz_1:0.1):0.1);'
	f = write_file(tmp_path, 'OG1.tre', newick)

	kept = contamination.get_sisters(make_args(), f, SISTER_RULES, {})

	assert 'Sr_ci_AAAA_1' in kept


def test_get_sisters_no_rules_keeps_everything(tmp_path):
	f = write_file(tmp_path, 'OG1.tre', TREE_CILIATE_SISTER_TO_BACTERIUM)

	kept = contamination.get_sisters(make_args(), f, {}, {})

	assert sorted(kept) == leaf_names(ete3.Tree(TREE_CILIATE_SISTER_TO_BACTERIUM))


def test_get_sisters_cocontaminants_are_not_counted_as_sisters(tmp_path):
	# Declaring A and the bacterium as co-contaminants (the same "cell") means the bacterium
	# is no longer A's sister, so A is kept
	f = write_file(tmp_path, 'OG1.tre', TREE_CILIATE_SISTER_TO_BACTERIUM)
	coconts = write_file(tmp_path, 'coconts.tsv', 'Sr_ci_AAAA\tBa_bb_Ecol\n')

	kept = contamination.get_sisters(make_args(cocontaminants=coconts), f, SISTER_RULES, {})

	assert 'Sr_ci_AAAA_1' in kept


def test_get_sisters_missing_cocontaminants_file_exits(tmp_path):
	f = write_file(tmp_path, 'OG1.tre', TREE_CILIATE_SISTER_TO_BACTERIUM)

	with pytest.raises(SystemExit):
		contamination.get_sisters(make_args(cocontaminants=str(tmp_path / 'missing.tsv')), f, SISTER_RULES, {})


def test_get_sisters_subsister_rule(tmp_path):
	# A's sister is Op (not a contaminant), but its sub-sister is Ba, which the sub-sister rule removes
	newick = '(((Sr_ci_AAAA_1:0.1,Op_me_Hsap_1:0.1):0.1,Ba_bb_Ecol_1:0.1):0.1,(Am_tu_Xxxx_1:0.1,Pl_gr_Atha_1:0.1):0.1,(Ex_eu_Yyyy_1:0.1,Ex_eu_Zzzz_1:0.1):0.1);'
	f = write_file(tmp_path, 'OG1.tre', newick)

	kept = contamination.get_sisters(make_args(subsister_rules='rules.tsv'), f, {}, {'Sr_ci': ['Ba']})

	assert 'Sr_ci_AAAA_1' not in kept
	assert len(kept) == 6


@pytest.mark.xfail(strict=True, reason='Known bug: the sister search walks past the root when every leaf belongs to one taxon')
def test_get_sisters_single_taxon_tree(tmp_path):
	newick = '((Sr_ci_AAAA_1:0.1,Sr_ci_AAAA_2:0.1):0.1,Sr_ci_AAAA_3:0.1);'
	f = write_file(tmp_path, 'OG1.tre', newick)

	kept = contamination.get_sisters(make_args(), f, SISTER_RULES, {})

	assert sorted(kept) == ['Sr_ci_AAAA_1', 'Sr_ci_AAAA_2', 'Sr_ci_AAAA_3']


# --------------------------------------------------------------------------------------
# test write_new_preguidance
# --------------------------------------------------------------------------------------

@pytest.fixture
def output_dir(tmp_path):
	(tmp_path / 'Output' / 'Pre-Guidance').mkdir(parents=True)
	return tmp_path


SEQS_PER_OG = {
	'OG1_preguidance.fasta': {
		'Sr_ci_AAAA_1': 'MKVL',
		'Op_me_Hsap_1': 'MKIL',
	}
}


def read_fasta_ids(path):
	return [line[1:].strip() for line in open(path) if line.startswith('>')]


def test_write_new_preguidance_writes_only_kept_sequences(output_dir):
	params = make_args(output=str(output_dir))

	seq_file, removed = contamination.write_new_preguidance(params, ['Sr_ci_AAAA_1'], SEQS_PER_OG, 'OG1_preguidance.FastTree.tre')

	assert seq_file == 'OG1_preguidance.fasta'
	assert removed == ['Op_me_Hsap_1']
	assert read_fasta_ids(output_dir / 'Output' / 'Pre-Guidance' / seq_file) == ['Sr_ci_AAAA_1']


@pytest.mark.parametrize('exclude_code', ['Op', 'Op_me', 'Op_me_Hsap'])
def test_write_new_preguidance_excludes_taxa(output_dir, tmp_path, exclude_code):
	# Codes can be a major clade (2 chars), minor clade (5 chars) or full taxon code (10 chars)
	exclude = write_file(tmp_path, 'exclude.txt', exclude_code + '\n')
	params = make_args(output=str(output_dir), cl_exclude_taxa=exclude)

	seq_file, removed = contamination.write_new_preguidance(params, ['Sr_ci_AAAA_1', 'Op_me_Hsap_1'], SEQS_PER_OG, 'OG1_preguidance.FastTree.tre')

	assert read_fasta_ids(output_dir / 'Output' / 'Pre-Guidance' / seq_file) == ['Sr_ci_AAAA_1']


def test_write_new_preguidance_no_matching_sequence_file(output_dir):
	params = make_args(output=str(output_dir))

	assert contamination.write_new_preguidance(params, [], SEQS_PER_OG, 'OG2.FastTree.tre') == (None, [])


def test_write_new_preguidance_multiple_matching_sequence_files(output_dir):
	params = make_args(output=str(output_dir))
	seqs_per_og = {'OG1_a.fasta': {}, 'OG1_b.fasta': {}}

	assert contamination.write_new_preguidance(params, [], seqs_per_og, 'OG1.FastTree.tre') == (None, [])


def test_write_new_preguidance_missing_exclude_file_exits(output_dir):
	params = make_args(output=str(output_dir), cl_exclude_taxa=str(output_dir / 'missing.txt'))

	with pytest.raises(SystemExit):
		contamination.write_new_preguidance(params, [], SEQS_PER_OG, 'OG1_preguidance.FastTree.tre')


# --------------------------------------------------------------------------------------
# test external tool wrappers (the tools themselves are mocked, never run)
# --------------------------------------------------------------------------------------

def test_cl_fasttree_runs_on_fasta_files_only(tmp_path):
	guidance_dir = tmp_path / 'Output' / 'Guidance'
	guidance_dir.mkdir(parents=True)
	(guidance_dir / 'OG1.70gapTrimmed.fasta').write_text('>a\nMK\n')
	(guidance_dir / 'notes.txt').write_text('ignore me')
	params = make_args(output=str(tmp_path))

	with mock.patch('contamination.os.system') as system:
		contamination.cl_fasttree(params)

	system.assert_called_once_with('FastTree ' + str(tmp_path) + '/Output/Guidance/OG1.70gapTrimmed.fasta > ' + str(tmp_path) + '/Output/Trees/OG1.FastTree.tre')
